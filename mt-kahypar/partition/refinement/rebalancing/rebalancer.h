/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include <queue>
#include <boost/dynamic_bitset.hpp>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"

namespace mt_kahypar {
template <template <typename> class GainPolicy>
class Rebalancer {
 private:
  using GainCalculator = GainPolicy<PartitionedHypergraph>;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static constexpr Gain MIN_PQ_GAIN_THRESHOLD = 5;

  struct MoveGainComparator {
    bool operator()(const Move& lhs, const Move& rhs) {
      return lhs.gain > rhs.gain || (lhs.gain == rhs.gain && lhs.node < rhs.node);
    }
  };

  using MovePQ = std::priority_queue<Move, vec<Move>, MoveGainComparator>;

  struct IndexedMovePQ {
    explicit IndexedMovePQ(const size_t idx) :
      idx(idx),
      pq() { }

    size_t idx;
    MovePQ pq;
  };

 public:
  explicit Rebalancer(PartitionedHypergraph& hypergraph,
                       const Context& context,
                       const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _gain(context),
    _part_weights(_context.partition.k) { }

  Rebalancer(const Rebalancer&) = delete;
  Rebalancer(Rebalancer&&) = delete;

  Rebalancer & operator= (const Rebalancer &) = delete;
  Rebalancer & operator= (Rebalancer &&) = delete;

  void rebalance(kahypar::Metrics& best_metrics) {
    // If partition is imbalanced, rebalancer is activated
    if ( !metrics::isBalanced(_hg, _context) ) {
      _gain.reset();

      // This function is passed as lambda to the changeNodePart function and used
      // to calculate the "real" delta of a move (in terms of the used objective function).
      auto objective_delta = [&](const HyperedgeID he,
                                const HyperedgeWeight edge_weight,
                                const HypernodeID edge_size,
                                const HypernodeID pin_count_in_from_part_after,
                                const HypernodeID pin_count_in_to_part_after) {
                              _gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                                              pin_count_in_from_part_after,
                                                              pin_count_in_to_part_after);
                            };




      vec<Move> moves_to_empty_blocks = repairEmptyBlocks();
      for (Move& m : moves_to_empty_blocks) {
        moveVertex(m.node, m, objective_delta);
      }

      // We first try to perform moves that does not worsen solution quality of the partition
      // Moves that would worsen the solution quality are gathered in a thread local priority queue
      // and processed afterwards if partition is still imbalanced
      std::atomic<size_t> idx(0);
      tbb::enumerable_thread_specific<IndexedMovePQ> move_pqs([&] {
        return IndexedMovePQ(idx++);
      });
      _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID from = _hg.partID(hn);
        if ( _hg.isBorderNode(hn) && _hg.partWeight(from) > _context.partition.max_part_weights[from] ) {
          Move rebalance_move = _gain.computeMaxGainMove(_hg, hn, true /* rebalance move */);
          if ( rebalance_move.gain <= 0 ) {
            moveVertex(hn, rebalance_move, objective_delta);
          } else if ( rebalance_move.gain != std::numeric_limits<Gain>::max() ) {
            move_pqs.local().pq.emplace(std::move(rebalance_move));
          }
        }
      });

      ASSERT([&] {
        for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
          if ( _part_weights[block] != _hg.partWeight(block) ) {
            return false;
          }
        }
        return true;
      }(), "Rebalancer part weights are wrong");

      // If partition is still imbalanced, we try execute moves stored into
      // the thread local priority queue which could possibly worsen solution quality
      if ( !metrics::isBalanced(_hg, _context) ) {

        // Initialize minimum gain value of each priority queue
        parallel::scalable_vector<uint8_t> active_pqs(idx.load(), false);
        parallel::scalable_vector<Gain> min_pq_gain(idx.load(),
          std::numeric_limits<Gain>::max() - MIN_PQ_GAIN_THRESHOLD);
        for ( const IndexedMovePQ& idx_pq : move_pqs ) {
          if ( !idx_pq.pq.empty() ) {
            min_pq_gain[idx_pq.idx] = idx_pq.pq.top().gain;
          }
        }

        // Function returns minimum gain value of all priority queues
        auto global_pq_min_gain = [&](const bool only_active_pqs) {
          Gain min_gain = std::numeric_limits<Gain>::max() - MIN_PQ_GAIN_THRESHOLD;
          for ( size_t i = 0; i < min_pq_gain.size(); ++i ) {
            if ( (!only_active_pqs || active_pqs[i]) && min_pq_gain[i] < min_gain ) {
              min_gain = min_pq_gain[i];
            }
          }
          return min_gain;
        };

        // We process each priority queue in parallel. When we perform
        // a move we make sure that the current minimum gain value of the local
        // PQ is within a certain threshold of the global minimum gain value.
        // Otherwise, we perform busy waiting until all moves with a better gain
        // are processed.
        tbb::parallel_for_each(move_pqs, [&](IndexedMovePQ& idx_pq) {
          const size_t idx = idx_pq.idx;
          MovePQ& pq = idx_pq.pq;
          active_pqs[idx] = true;
          Gain current_global_min_pq_gain = global_pq_min_gain(false);
          while ( !pq.empty() ) {
            Move move = pq.top();
            min_pq_gain[idx] = move.gain;
            pq.pop();

            // If the minimum gain value of the local priority queue is not within
            // a certain threshold of the global priority queue, we perform busy waiting
            // until all moves with a better gain of other pqs are performed.
            while ( move.gain > current_global_min_pq_gain + MIN_PQ_GAIN_THRESHOLD ) {
              current_global_min_pq_gain = global_pq_min_gain(true);
            }

            const PartitionID from = move.from;
            if ( _hg.partWeight(from) > _context.partition.max_part_weights[from] ) {
              Move real_move = _gain.computeMaxGainMove(_hg, move.node, true /* rebalance move */);
              if ( real_move.gain <= move.gain ) {
                moveVertex(real_move.node, real_move, objective_delta);
              } else if ( real_move.gain != std::numeric_limits<Gain>::max() ) {
                pq.emplace(std::move(real_move));
              }
            }
          }
          active_pqs[idx] = false;
          min_pq_gain[idx] = std::numeric_limits<Gain>::max() - MIN_PQ_GAIN_THRESHOLD;
        });
      }

      // Update metrics statistics
      HyperedgeWeight current_metric = best_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      Gain delta = _gain.delta();
      HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
        metrics::objective(_hg, _context.partition.objective),
        V(current_metric) << V(delta) <<
        V(metrics::objective(_hg, _context.partition.objective)));
      best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
    }
  }

 private:

  vec<Move> repairEmptyBlocks() {
    // First detect if there are any empty blocks.
    boost::dynamic_bitset<> is_empty(_context.partition.k);
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      _part_weights[block] = _hg.partWeight(block);
      if (_hg.partWeight(block) == 0) {
        is_empty.set(block, true);
      }
    }

    vec<Move> moves_to_empty_blocks;

    // If so, find the best vertices to move to that block
    while (is_empty.any()) {
      const PartitionID k = _context.partition.k;
      tbb::enumerable_thread_specific< vec<Gain> > ets_scores(k, 0);

      // positive gain values correspond to "good" improvement. MovePQ uses std::greater (MinHeap)
      // --> stores worst gains at the top where we can eject them
      tbb::enumerable_thread_specific< vec< vec<Move> > > ets_best_move(k);

      _hg.doParallelForAllNodes([&](const HypernodeID u) {
        vec<Gain>& scores = ets_scores.local();
        vec< vec<Move> >& move_proposals = ets_best_move.local();

        const PartitionID pu = _hg.partID(u);
        Gain unremovable = 0;
        for (HyperedgeID e : _hg.incidentEdges(u)) {
          const HyperedgeWeight edge_weight = _hg.edgeWeight(e);
          if (_hg.pinCountInPart(e, pu) > 1) {
            unremovable += edge_weight;
          }
          for (PartitionID i : _hg.connectivitySet(e)) {
            scores[i] += edge_weight;
          }
        }

        // maintain thread local priority queues of up to k best gains
        for (PartitionID i = 0; i < k; ++i) {
          if (i != pu && is_empty[i] && _hg.partWeight(pu) > _hg.nodeWeight(u)) {
            const Gain gain = unremovable - scores[i];
            vec<Move>& c = move_proposals[i];
            if (c.size() < k) {
              c.push_back(Move { pu, i, u, gain });
              std::push_heap(c.begin(), c.end(), MoveGainComparator());
            } else if (c.front().gain < gain) {
              std::pop_heap(c.begin(), c.end(), MoveGainComparator());
              c.back() = { pu, i, u, gain };
              std::push_heap(c.begin(), c.end(), MoveGainComparator());
            }
          }
          scores[i] = 0;
        }
      });

      vec< vec<Move> > best_moves_per_part(k);

      for (vec<vec<Move>>& tlpq : ets_best_move) {
        size_t i = is_empty.find_first();
        while (i != is_empty.npos) {
          std::copy(tlpq[i].begin(), tlpq[i].end(), std::back_inserter(best_moves_per_part[i]));
          i = is_empty.find_next(i);
        }
      }

      auto prefer_highest_gain = [&](const Move& lhs, const Move& rhs) {
        return lhs.gain > rhs.gain || (lhs.gain == rhs.gain &&
                                       _hg.partWeight(_hg.partID(lhs.node) > _hg.partWeight(_hg.partID(rhs.node))));

      };

      size_t i = is_empty.find_first();
      while (i != is_empty.npos) {
        vec<Move>& c = best_moves_per_part[i];
        std::nth_element(c.begin(), c.begin() + k, c.end(), prefer_highest_gain);
        c.erase(c.begin() + k, c.end());
        i = is_empty.find_next(i);
      }

      // do greedy matching



    }
    return moves_to_empty_blocks;
  }

  template<typename F>
  bool moveVertex(const HypernodeID hn, const Move& move, const F& objective_delta) {
    ASSERT(_hg.partID(hn) == move.from);
    const PartitionID from = move.from;
    const PartitionID to = move.to;
    const HypernodeWeight node_weight = _hg.nodeWeight(hn);
    if ( from != to ) {
      // Before moving, we ensure that the block we move the vertex to does
      // not become overloaded
      _part_weights[to] += node_weight;
      if ( _part_weights[to] <= _context.partition.max_part_weights[to] ) {
        if ( _hg.changeNodePart(hn, from, to, objective_delta) ) {
          DBG << "Moved vertex" << hn << "from block" << from << "to block" << to
              << "with gain" << move.gain;
          _part_weights[from] -= node_weight;
          return true;
        }
      }
      _part_weights[to] -= node_weight;
    }
    return false;
  }

  PartitionedHypergraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  GainCalculator _gain;
  parallel::scalable_vector<AtomicWeight> _part_weights;
};

using Km1Rebalancer = Rebalancer<Km1Policy>;
using CutRebalancer = Rebalancer<CutPolicy>;
}  // namespace kahypar
