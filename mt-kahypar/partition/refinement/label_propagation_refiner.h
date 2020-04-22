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
          template <typename> class GainPolicy>
class LabelPropagationRefinerT final : public IRefinerT<TypeTraits> {
 private:
  using HyperGraph = typename TypeTraits::PartitionedHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using GainCalculator = GainPolicy<HyperGraph>;
  using NumaActiveNodes = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using NumaNextActiveNodes = parallel::scalable_vector<ds::StreamingVector<HypernodeID>>;
  using NumaFastResetFlagArray = parallel::scalable_vector<kahypar::ds::FastResetFlagArray<>>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefinerT(HyperGraph&,
                                    const Context& context,
                                    const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _is_numa_aware(false),
    _current_num_nodes(kInvalidHypernode),
    _current_num_edges(kInvalidHyperedge),
    _gain(context),
    _active_nodes(),
    _next_active(),
    _visited_he(),
    _numa_lp_round_synchronization(0) {
    _is_numa_aware = _context.refinement.label_propagation.numa_aware &&
                     !_context.refinement.label_propagation.execute_sequential &&
                     TBB::instance().num_used_numa_nodes() > 1;
  }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT & operator= (const LabelPropagationRefinerT &) = delete;
  LabelPropagationRefinerT & operator= (LabelPropagationRefinerT &&) = delete;

  ~LabelPropagationRefinerT() {
    parallel::parallel_free(_active_nodes, _next_active, _visited_he);
  }

 private:
  bool refineImpl(HyperGraph& hypergraph,
                  const parallel::scalable_vector<HypernodeID>&,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();

    _numa_lp_round_synchronization = 0;
    for ( int node = 0; node < static_cast<int>(_next_active.size()); ++node ) {
      _next_active[node].reset();
    }

    NumaNextActiveNodes next_active_nodes;
    if (_is_numa_aware) {
      // Execute label propagation on all numa nodes
      next_active_nodes.resize(TBB::instance().num_used_numa_nodes());
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
          // In case we execute numa-aware label propagation, we perform
          // label propagation on all nodes of the current numa node.
          ASSERT(node < static_cast<int>(_active_nodes.size()));
          labelPropagation(hypergraph, next_active_nodes, node);
        });
    } else {
      // In case, we execute non-numa-aware label propagation, we
      // perform label propagation on all vertices
      next_active_nodes.resize(1);
      labelPropagation(hypergraph, next_active_nodes, 0);
    }

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

  void labelPropagation(HyperGraph& hypergraph,
                        NumaNextActiveNodes& next_active_nodes,
                        const int node) {
    DBG << "LP refinement" << V(_context.partition.k) << V(_context.partition.objective) << V(hypergraph.initialNumNodes()) << V(TBB::instance().total_number_of_threads());
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i << "on NUMA Node" << node << V(_active_nodes[node].size())
              << V(metrics::km1(hypergraph)) << V(metrics::hyperedgeCut(hypergraph));

      utils::Timer::instance().start_timer(
        "lp_round_" + std::to_string(i) + (node != -1 ? "_" + std::to_string(node) : ""),
        "Label Propagation Round " + std::to_string(i) + (node != -1 ? " - " + std::to_string(node) : ""), true);

      if ( _active_nodes[node].size() > 0 ) {
        labelPropagationRound(hypergraph, next_active_nodes, node);
      }

      if ( _is_numa_aware ) {
        // In case, we execute numa-aware label propagation, we need to synchronize
        // the rounds in order that swapping the active node set is thread-safe
        ++_numa_lp_round_synchronization;
        const size_t synchro_barrier = ( i + 1 ) * TBB::instance().num_used_numa_nodes();
        while ( _numa_lp_round_synchronization.load() < synchro_barrier) { }
      }

      if ( _context.refinement.label_propagation.execute_sequential ) {
        _active_nodes[node] = next_active_nodes[node].copy_sequential();
        next_active_nodes[node].clear_sequential();
      } else {
        _active_nodes[node] = next_active_nodes[node].copy_parallel();
        next_active_nodes[node].clear_parallel();
      }
      utils::Timer::instance().stop_timer(
        "lp_round_" + std::to_string(i) + (node != -1 ? "_" + std::to_string(node) : ""));
    }
  }

  bool labelPropagationRound(HyperGraph& hypergraph,
                             NumaNextActiveNodes& next_active_nodes,
                             const int node) {
    _visited_he[node].reset();
    _next_active[node].reset();
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
        _active_nodes[node], 0UL, _active_nodes[node].size(), sched_getcpu());

      for ( size_t j = 0; j < _active_nodes[node].size(); ++j ) {
        const HypernodeID hn = _active_nodes[node][j];
        converged &= !moveVertex(hypergraph, hn,
          next_active_nodes, node, objective_delta);
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
        _active_nodes[node], 0UL, _active_nodes[node].size());

      tbb::parallel_for(0UL, _active_nodes[node].size(), [&](const size_t& j) {
          const HypernodeID hn = _active_nodes[node][j];
          converged &= !moveVertex(hypergraph, hn,
            next_active_nodes, node, objective_delta);
        });
    }
    return converged;
  }

  template<typename F>
  bool moveVertex(HyperGraph& hypergraph,
                  const HypernodeID hn,
                  NumaNextActiveNodes& next_active_nodes,
                  const int node,
                  const F& objective_delta) {
    bool is_moved = false;
    ASSERT(hn != kInvalidHypernode);
    if ( hypergraph.isBorderNode(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(!_is_numa_aware || common::get_numa_node_of_vertex(hn) == node);

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
        bool changed_part = hypergraph.changeNodePart(hn, from, to, objective_delta);
        if (changed_part) {
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
                const HyperedgeID original_he_id = hypergraph.originalEdgeID(he);
                if ( !_visited_he[node][original_he_id] ) {
                  for (const HypernodeID& pin : hypergraph.pins(he)) {
                    const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
                    int pin_numa_node = _is_numa_aware ?
                      common::get_numa_node_of_vertex(pin) : 0;
                    if ( (_context.refinement.label_propagation.rebalancing ||
                           hypergraph.isBorderNode(hn)) &&
                         !_next_active[pin_numa_node][original_pin_id] ) {
                      next_active_nodes[pin_numa_node].stream(pin);
                      _next_active[pin_numa_node].set(original_pin_id, true);
                    }
                  }
                  _visited_he[node].set(original_he_id, true);
                }
              }
            }
            if ( _next_active[node][original_id] ) {
              next_active_nodes[node].stream(hn);
              _next_active[node].set(original_id, true);
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
    const int num_numa_nodes = _is_numa_aware ? TBB::instance().num_used_numa_nodes() : 1;
    NumaActiveNodes tmp_active_nodes(num_numa_nodes);
    _active_nodes = std::move(tmp_active_nodes);

    if ( _context.refinement.label_propagation.execute_sequential ) {
      // Setup active nodes sequential
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        if ( _context.refinement.label_propagation.rebalancing ||
             hypergraph.isBorderNode(hn) ) {
          _active_nodes[0].push_back(hn);
        }
      }
    } else {
      // Setup active nodes in parallel
      // A node is active, if it is a border vertex.
      NumaNextActiveNodes tmp_active_nodes(num_numa_nodes);

      auto add_vertex = [&](const HypernodeID& hn) {
        const int node =  _is_numa_aware ?
          common::get_numa_node_of_vertex(hn) : 0;
        ASSERT(node < static_cast<int>(tmp_active_nodes.size()));
        tmp_active_nodes[node].stream(hn);
      };

      hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID& hn) {
        if ( _context.refinement.label_propagation.rebalancing ||
             hypergraph.isBorderNode(hn) ) {
          add_vertex(hn);
        }
      });

      if ( _is_numa_aware ) {
        TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
          ASSERT(node < static_cast<int>(_active_nodes.size()));
          _active_nodes[node] = tmp_active_nodes[node].copy_parallel();
        });
      } else {
        _active_nodes[0] = tmp_active_nodes[0].copy_parallel();
      }
    }

    const HypernodeID new_num_nodes = hypergraph.initialNumNodes();
    if ( new_num_nodes > _current_num_nodes || _current_num_nodes == kInvalidHypernode ) {
      _next_active.clear();
      for ( int node = 0; node < num_numa_nodes; ++node ) {
        _next_active.emplace_back(hypergraph.initialNumNodes());
      }
      _current_num_nodes = new_num_nodes;
    }

    const HypernodeID new_num_edges = hypergraph.initialNumEdges();
    if ( new_num_edges > _current_num_edges || _current_num_edges == kInvalidHyperedge ) {
      _visited_he.clear();
      for ( int node = 0; node < num_numa_nodes; ++node ) {
        _visited_he.emplace_back(hypergraph.initialNumEdges());
      }
      _current_num_edges = new_num_edges;
    }
  }

  const Context& _context;
  const TaskGroupID _task_group_id;
  bool _is_numa_aware;
  HypernodeID _current_num_nodes;
  HyperedgeID _current_num_edges;
  GainCalculator _gain;
  NumaActiveNodes _active_nodes;
  NumaFastResetFlagArray _next_active;
  NumaFastResetFlagArray _visited_he;
  std::atomic<size_t> _numa_lp_round_synchronization;
};

using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, Km1Policy>;
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, CutPolicy>;
}  // namespace kahypar
