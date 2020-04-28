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
  using GainCalculator = GainPolicy<PartitionedHypergraph<>>;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit Rebalancer(PartitionedHypergraph<>& hypergraph,
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
    if ( metrics::imbalance(_hg, _context) > _context.partition.epsilon ) {
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

      // We first try to perform moves that does not worsen solution quality of the partition
      _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID from = _hg.partID(hn);
        if ( _hg.isBorderNode(hn) && _hg.partWeight(from) > _context.partition.max_part_weights[from] ) {
          Move rebalance_move = _gain.computeMaxGainMove(_hg, hn, true /* rebalance move */);
          if ( rebalance_move.gain <= 0 ) {
            moveVertex(hn, rebalance_move, objective_delta);
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

      // If partition is still imbalanced, we just move the first vertices we find
      // from an overloaded block to an other block that maximizes the objective function.
      if ( metrics::imbalance(_hg, _context) > _context.partition.epsilon ) {
        _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          const PartitionID from = _hg.partID(hn);
          if ( _hg.partWeight(from) > _context.partition.max_part_weights[from] ) {
            moveVertex(hn, _gain.computeMaxGainMove(_hg, hn, true /* rebalance move */), objective_delta);
          }
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

  PartitionedHypergraph<>& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  GainCalculator _gain;
  parallel::scalable_vector<AtomicWeight> _part_weights;
};

using Km1Rebalancer = Rebalancer<Km1Policy>;
using CutRebalancer = Rebalancer<CutPolicy>;
}  // namespace kahypar
