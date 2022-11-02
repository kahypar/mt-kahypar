/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/rebalancing//rebalancer.h"


#include <boost/dynamic_bitset.hpp>

#include <tbb/parallel_for_each.h>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {

  template <template <typename> class GainPolicy>
  void Rebalancer<GainPolicy>::rebalance(Metrics& best_metrics) {
    // If partition is imbalanced, rebalancer is activated
    if ( !metrics::isBalanced(_hg, _context) ) {
      _gain.reset();

      for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
        _part_weights[block] = _hg.partWeight(block);
      }

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
              Mode::direct, _context.partition.objective);
      Gain delta = _gain.delta();
      HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
                              metrics::objective(_hg, _context.partition.objective),
                              V(current_metric) << V(delta) <<
                                                V(metrics::objective(_hg, _context.partition.objective)));
      best_metrics.updateMetric(current_metric + delta, Mode::direct, _context.partition.objective);
    }
  }

  template <template <typename> class GainPolicy>
  vec<Move> Rebalancer<GainPolicy>::repairEmptyBlocks() {
    // First detect if there are any empty blocks.
    const size_t k = size_t(_context.partition.k);
    boost::dynamic_bitset<> is_empty(k);
    for (size_t i = 0; i < k; ++i) {
      if (_hg.partWeight(PartitionID(i)) == 0) {
        is_empty.set(i, true);
      }
    }

    vec<Move> moves_to_empty_blocks;

    // If so, find the best vertices to move to that block
    while (is_empty.any()) {

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
        for (PartitionID i = 0; i < PartitionID(k); ++i) {
          if (i != pu && is_empty[i] && _hg.partWeight(pu) > _hg.nodeWeight(u)
              && _hg.nodeWeight(u) <= _context.partition.max_part_weights[i]) {
            const Gain gain = scores[i] - unremovable;
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
        const HypernodeWeight pwl = _hg.partWeight(_hg.partID(lhs.node));
        const HypernodeWeight pwr = _hg.partWeight(_hg.partID(rhs.node));
        return std::tie(lhs.gain, pwl, lhs.node) > std::tie(rhs.gain, pwr, rhs.node);
      };

      auto node_already_used = [&](HypernodeID node) {
        return std::any_of(moves_to_empty_blocks.begin(),
                           moves_to_empty_blocks.end(),
                           [node](const Move& m) { return m.node == node; }
        );
      };

      size_t i = is_empty.find_first();
      while (i != is_empty.npos) {
        vec<Move>& c = best_moves_per_part[i];
        std::sort(c.begin(), c.end(), prefer_highest_gain);

        size_t j = 0;
        while (j < c.size() && node_already_used(c[j].node)) { ++j; }
        if (j != c.size()) {
          moves_to_empty_blocks.push_back(c[j]);
          is_empty.set(i, false);
        }

        i = is_empty.find_next(i);
      }

    }
    return moves_to_empty_blocks;
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  template class Rebalancer<Km1Policy>;
  template class Rebalancer<CutPolicy>;
}