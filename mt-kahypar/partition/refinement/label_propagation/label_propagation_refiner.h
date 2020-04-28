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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"

#include "kahypar/meta/mandatory.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <typename TypeTraits,
          template <typename> class GainPolicy,
          bool track_border_vertices = TRACK_BORDER_VERTICES>
class LabelPropagationRefinerT final : public IRefinerT<TypeTraits, track_border_vertices> {
 private:
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<track_border_vertices>;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using GainCalculator = GainPolicy<HyperGraph>;
  using ActiveNodes = parallel::scalable_vector<HypernodeID>;
  using NextActiveNodes = ds::StreamingVector<HypernodeID>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefinerT(HyperGraph&,
                                    const Context& context,
                                    const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _current_num_nodes(kInvalidHypernode),
    _current_num_edges(kInvalidHyperedge),
    _gain(context),
    _active_nodes(),
    _next_active(),
    _visited_he() { }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT & operator= (const LabelPropagationRefinerT &) = delete;
  LabelPropagationRefinerT & operator= (LabelPropagationRefinerT &&) = delete;

 private:
  bool refineImpl(HyperGraph& hypergraph,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();
    _next_active.reset();

    // Perform Label Propagation
    labelPropagation(hypergraph);

    // Update global part weight and sizes
    best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

    // Update metrics statistics
    HyperedgeWeight current_metric = best_metrics.getMetric(
      kahypar::Mode::direct_kway, _context.partition.objective);
    Gain delta = _gain.delta();
    ASSERT(delta <= 0, "LP refiner worsen solution quality");
    HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
      metrics::objective(hypergraph, _context.partition.objective),
      V(current_metric) << V(delta) <<
      V(metrics::objective(hypergraph, _context.partition.objective)));
    best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
    utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
    return delta < 0;
  }

  void labelPropagation(HyperGraph& hypergraph) {
    NextActiveNodes next_active_nodes;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i;

      utils::Timer::instance().start_timer(
        "lp_round_" + std::to_string(i), "Label Propagation Round " + std::to_string(i), true);

      if ( _active_nodes.size() > 0 ) {
        labelPropagationRound(hypergraph, next_active_nodes);
      }

      if ( _context.refinement.label_propagation.execute_sequential ) {
        _active_nodes = next_active_nodes.copy_sequential();
        next_active_nodes.clear_sequential();
      } else {
        _active_nodes = next_active_nodes.copy_parallel();
        next_active_nodes.clear_parallel();
      }
      utils::Timer::instance().stop_timer("lp_round_" + std::to_string(i));

      if ( _active_nodes.size() == 0 ) {
        break;
      }
    }
  }

  bool labelPropagationRound(HyperGraph& hypergraph,
                             NextActiveNodes& next_active_nodes) {
    _visited_he.reset();
    _next_active.reset();
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
                             _gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                                            pin_count_in_from_part_after, pin_count_in_to_part_after);
                           };

    // Shuffle Vector
    bool converged = true;
    if ( _context.refinement.label_propagation.execute_sequential ) {
      utils::Randomize::instance().shuffleVector(
        _active_nodes, 0UL, _active_nodes.size(), sched_getcpu());

      for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
        const HypernodeID hn = _active_nodes[j];
        converged &= !moveVertex(hypergraph, hn,
          next_active_nodes, objective_delta);
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
        _active_nodes, 0UL, _active_nodes.size());

      tbb::parallel_for(0UL, _active_nodes.size(), [&](const size_t& j) {
          const HypernodeID hn = _active_nodes[j];
          converged &= !moveVertex(hypergraph, hn,
            next_active_nodes, objective_delta);
        });
    }
    return converged;
  }

  template<typename F>
  bool moveVertex(HyperGraph& hypergraph,
                  const HypernodeID hn,
                  NextActiveNodes& next_active_nodes,
                  const F& objective_delta) {
    bool is_moved = false;
    ASSERT(hn != kInvalidHypernode);
    if ( hypergraph.isBorderNode(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));

      Move best_move = _gain.computeMaxGainMove(hypergraph, hn);
      // We perform a move if it either improves the solution quality or, in case of a
      // zero gain move, the balance of the solution.
      const bool positive_gain = best_move.gain < 0;
      const bool zero_gain_move = (_context.refinement.label_propagation.rebalancing &&
                                    best_move.gain == 0 &&
                                    hypergraph.partWeight(best_move.from) - 1 >
                                    hypergraph.partWeight(best_move.to) + 1 &&
                                    hypergraph.partWeight(best_move.to) <
                                    _context.partition.perfect_balance_part_weights[best_move.to]);
      const bool perform_move = positive_gain || zero_gain_move;
      if (best_move.from != best_move.to && perform_move) {
        PartitionID from = best_move.from;
        PartitionID to = best_move.to;

        Gain delta_before = _gain.localDelta();
        if (hypergraph.changeNodePart(hn, from, to, objective_delta)) {
          // In case the move to block 'to' was successful, we verify that the "real" gain
          // of the move is either equal to our computed gain or if not, still improves
          // the solution quality.
          Gain move_delta = _gain.localDelta() - delta_before;
          bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
          if (accept_move) {
            DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";

            // Set all neighbors of the vertex to active
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
              if ( hypergraph.edgeSize(he) <=
                    ID(_context.refinement.label_propagation.hyperedge_size_activation_threshold) ) {
                if ( !_visited_he[he] ) {
                  for (const HypernodeID& pin : hypergraph.pins(he)) {
                    if ( (_context.refinement.label_propagation.rebalancing ||
                           hypergraph.isBorderNode(hn)) &&
                         !_next_active[pin] ) {
                      next_active_nodes.stream(pin);
                      _next_active.set(pin, true);
                    }
                  }
                  _visited_he.set(he, true);
                }
              }
            }
            if ( _next_active[hn] ) {
              next_active_nodes.stream(hn);
              _next_active.set(hn, true);
            }
            is_moved = true;
          } else {
            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
            // In case, the real gain is not equal with the computed gain and
            // worsen the solution quality we revert the move.
            ASSERT(hypergraph.partID(hn) == to);
            hypergraph.changeNodePart(hn, to, from, objective_delta);
          }
        }
      }
    }

    return is_moved;
  }

  void initializeImpl(HyperGraph& hypergraph) override final {
    ActiveNodes tmp_active_nodes;
    _active_nodes = std::move(tmp_active_nodes);

    if ( _context.refinement.label_propagation.execute_sequential ) {
      // Setup active nodes sequential
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        if ( _context.refinement.label_propagation.rebalancing ||
             hypergraph.isBorderNode(hn) ) {
          _active_nodes.push_back(hn);
        }
      }
    } else {
      // Setup active nodes in parallel
      // A node is active, if it is a border vertex.
      NextActiveNodes tmp_active_nodes;

      auto add_vertex = [&](const HypernodeID& hn) {
        tmp_active_nodes.stream(hn);
      };

      hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID& hn) {
        if ( _context.refinement.label_propagation.rebalancing ||
             hypergraph.isBorderNode(hn) ) {
          add_vertex(hn);
        }
      });

      _active_nodes = tmp_active_nodes.copy_parallel();
    }

    const HypernodeID new_num_nodes = hypergraph.initialNumNodes();
    if ( new_num_nodes > _current_num_nodes || _current_num_nodes == kInvalidHypernode ) {
      _next_active = kahypar::ds::FastResetFlagArray<>(hypergraph.initialNumNodes());
      _current_num_nodes = new_num_nodes;
    }

    const HypernodeID new_num_edges = hypergraph.initialNumEdges();
    if ( new_num_edges > _current_num_edges || _current_num_edges == kInvalidHyperedge ) {
      _visited_he = kahypar::ds::FastResetFlagArray<>(hypergraph.initialNumEdges());
      _current_num_edges = new_num_edges;
    }
  }

  const Context& _context;
  const TaskGroupID _task_group_id;
  HypernodeID _current_num_nodes;
  HyperedgeID _current_num_edges;
  GainCalculator _gain;
  ActiveNodes _active_nodes;
  kahypar::ds::FastResetFlagArray<> _next_active;
  kahypar::ds::FastResetFlagArray<> _visited_he;
};

using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, Km1Policy>;
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, CutPolicy>;
}  // namespace kahypar
