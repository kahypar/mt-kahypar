/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/node_swapper/node_swapper.h"

#include <queue>

#include "tbb/concurrent_vector.h"
#include "tbb/parallel_sort.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"


namespace mt_kahypar {

#define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }


namespace {
struct VertexMove {
  HypernodeID hn;
  PartitionID from;
  Gain gain;
};

class VertexMoveComparator {
public:
  bool operator()(const VertexMove& lhs, const VertexMove& rhs) {
    return lhs.gain < rhs.gain;
  }
};

// Simple concurrent priority queue.
// Our concurrent priority queue contains p (number of threads) priority queues.
// push and pop selects a random priority queue and protect writes to it via a spin lock.
class ConcurrentPQ {

using PQ = std::priority_queue<VertexMove, std::vector<VertexMove>, VertexMoveComparator>;

public:
  ConcurrentPQ(const size_t num_threads) :
    _num_threads(num_threads),
    _pq_lock(num_threads),
    _pqs(num_threads) { }

  void push(VertexMove&& move) {
    const size_t i = selectPQ();
    _pq_lock[i].lock();
    _pqs[i].emplace(std::move(move));
    _pq_lock[i].unlock();
  }

  VertexMove pop() {
    const size_t i = selectPQ();
    _pq_lock[i].lock();
    if ( !_pqs[i].empty() ) {
      VertexMove move = _pqs[i].top();
      _pqs[i].pop();
      _pq_lock[i].unlock();
      return move;
    }
    _pq_lock[i].unlock();
    return VertexMove { kInvalidHypernode, kInvalidPartition, kInvalidGain };
  }

private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t selectPQ() {
    return utils::Randomize::instance().getRandomInt(
      0, static_cast<int>(_num_threads) - 1, sched_getcpu());
  }

  const size_t _num_threads;
  vec<SpinLock> _pq_lock;
  vec<PQ> _pqs;
};

} // namespace

HyperedgeWeight NodeSwapper::refine() {
  if ( !_hg.isGainCacheInitialized() ) {
    _hg.initializeGainCache();
  }
  ASSERT(_hg.checkTrackedPartitionInformation());

  CAtomic<HyperedgeWeight> delta(0);
  const HyperedgeWeight initial_km1 = metrics::km1(_hg);
  HyperedgeWeight last_round_km1 = initial_km1;
  double last_round_reduction = 2.0;

  // We perform several rounds until the overall improvement is below 0.25%
  while ( last_round_reduction > 0.0025 ) {
    _in_queue.reset();

    // First, we compute for each node its preferred block that improves the objective function
    // and insert it into the corresponding priority queue.
    vec<ConcurrentPQ> to_pq(_context.partition.k, ConcurrentPQ(_context.shared_memory.num_threads));
    _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID from = _hg.partID(hn);
      PartitionID best_to = from;
      Gain best_gain = 0;
      for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
        if ( from != to ) {
          const Gain gain = _hg.km1Gain(hn, from, to);
          if ( gain > best_gain || ( gain > 0 && gain == best_gain &&
              utils::Randomize::instance().flipCoin(sched_getcpu()) )) {
            best_gain = gain;
            best_to = to;
          }
        }
      }

      if ( from != best_to ) {
        to_pq[best_to].push(VertexMove { hn, from, best_gain });
        _in_queue.set(hn);
      }
    });

    // Second, for each node u not contained in one of the priority queue,
    // we try to swap it with the highest gain node that want to move into
    // the block of node u.
    utils::Randomize::instance().parallelShuffleVector(_nodes, 0UL, _nodes.size());
    tbb::parallel_for(0UL, _nodes.size(), [&](const size_t i) {
      const HypernodeID u = _nodes[i];
      if ( _hg.nodeIsEnabled(u) && !_in_queue[u] ) {
        const PartitionID u_from = _hg.partID(u);
        VertexMove move = to_pq[u_from].pop();
        while ( move.hn != kInvalidHypernode &&
                move.gain != _hg.km1Gain(move.hn, move.from, u_from) ) {
          move.gain = _hg.km1Gain(move.hn, move.from, u_from);
          if ( move.gain > 0 ) {
            to_pq[u_from].push(std::move(move));
          }
          move = to_pq[u_from].pop();
        }

        if ( move.hn != kInvalidHypernode ) {
          const HypernodeID v = move.hn;
          const PartitionID v_from = move.from;
          const Gain u_gain = _hg.km1Gain(u, u_from, v_from);

          // Note that the actual swap gain can differ, if u and v share
          // some common edges. However, we double check the estimated gain
          // with attributed gain and if the swap worsen the solution quality
          // we revert it immediatly.
          const Gain estimated_swap_gain = move.gain + u_gain;
          if ( estimated_swap_gain > 0 ) {
            Gain real_swap_gain = 0;
            auto delta_gain_func = [&](const HyperedgeID he,
                                      const HyperedgeWeight edge_weight,
                                      const HypernodeID edge_size,
                                      const HypernodeID pin_count_in_from_part_after,
                                      const HypernodeID pin_count_in_to_part_after) {
              real_swap_gain -= km1Delta(he, edge_weight, edge_size,
                pin_count_in_from_part_after, pin_count_in_to_part_after);
            };

            // Perform swap
            _hg.changeNodePartWithGainCacheUpdate(u, u_from, v_from,
              std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
            _hg.changeNodePartWithGainCacheUpdate(v, v_from, u_from,
              std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);

            if ( real_swap_gain < 0 ) {
              // Revert swap, if it has worsen the solution quality
              _hg.changeNodePartWithGainCacheUpdate(u, v_from, u_from,
                std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
              _hg.changeNodePartWithGainCacheUpdate(v, u_from, v_from,
                std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
              move.gain = _hg.km1Gain(move.hn, v_from, u_from);
              to_pq[u_from].push(std::move(move));
            }

            delta += real_swap_gain;
          } else {
            to_pq[u_from].push(std::move(move));
          }
        }
      }
    });

    _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      _hg.recomputeMoveFromBenefit(hn);
    });

    const HyperedgeWeight current_km1 = initial_km1 - delta.load(std::memory_order_relaxed);
    last_round_reduction = static_cast<double>(last_round_km1) / current_km1 - 1.0;
    last_round_km1 = current_km1;
  }

  return -delta.load(std::memory_order_relaxed);
}

} // namespace mt_kahypar