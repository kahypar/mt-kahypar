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
  using NumaFastResetFlagArray = parallel::scalable_vector<kahypar::ds::FastResetFlagArray<>>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefinerT(HyperGraph&,
                                    const Context& context,
                                    const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _numa_nodes_indices(),
    _nodes(),
    _gain(context),
    _active(),
    _next_active(),
    _visited_he(),
    _numa_lp_round_synchronization(0) { }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT & operator= (const LabelPropagationRefinerT &) = delete;
  LabelPropagationRefinerT & operator= (LabelPropagationRefinerT &&) = delete;

 private:
  bool refineImpl(HyperGraph& hypergraph,
                  const parallel::scalable_vector<HypernodeID>&,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();

    utils::Timer::instance().start_timer("label_propagation", "Label Propagation");

    HEAVY_REFINEMENT_ASSERT([&] {
        // Assertion verifies, that all enabled nodes are contained in _nodes
        parallel::scalable_vector<HypernodeID> tmp_nodes = _nodes;
        std::sort(tmp_nodes.begin(), tmp_nodes.end());
        size_t pos = 0;
        for (const HypernodeID& hn : hypergraph.nodes()) {
          HypernodeID current_hn = tmp_nodes[pos++];
          while ( !hypergraph.nodeIsEnabled(current_hn) ) {
            current_hn = tmp_nodes[pos++];
          }
          if (current_hn != hn) {
            LOG << "Hypernode" << hn << "not contained in vertex set";
            return false;
          }
        }
        return true;
      } (), "LP Refiner does not contain all vertices hypergraph");

    _numa_lp_round_synchronization = 0;
    for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
      _active[node].reset();
      _next_active[node].reset();
    }

    if (_context.refinement.label_propagation.numa_aware &&
        TBB::instance().num_used_numa_nodes() > 1) {
      // Execute label propagation on all numa nodes
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
          // In case we execute numa-aware label propagation, we perform
          // label propagation on all nodes of the current numa node.
          size_t start = _numa_nodes_indices[node];
          size_t end = _numa_nodes_indices[node + 1];
          labelPropagation(hypergraph, _nodes, start, end, node);
        });
    } else {
      // In case, we execute non-numa-aware label propagation, we
      // perform label propagation on all vertices
      labelPropagation(hypergraph, _nodes, 0UL, _nodes.size());
    }

    // Update global part weight and sizes
    best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

    // Update metrics statistics
    HyperedgeWeight current_metric = best_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective);
    Gain delta = _gain.delta();
    ASSERT(delta <= 0, "LP refiner worsen solution quality");
    HEAVY_REFINEMENT_ASSERT(current_metric + delta == metrics::objective(hypergraph, _context.partition.objective),
                            V(current_metric) << V(delta) << V(metrics::objective(hypergraph, _context.partition.objective)));
    best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
    utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
    utils::Timer::instance().stop_timer("label_propagation");
    return delta < 0;
  }

  void labelPropagation(HyperGraph& hypergraph,
                        parallel::scalable_vector<HypernodeID>& refinement_nodes,
                        const size_t start,
                        const size_t end,
                        int node = -1) {
    int numa_node = std::max(node, 0);
    bool converged = false;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i << "on NUMA Node" << node;

      utils::Timer::instance().start_timer(
        "lp_round_" + std::to_string(i) + (node != -1 ? "_" + std::to_string(node) : ""),
        "Label Propagation Round " + std::to_string(i) + (node != -1 ? " - " + std::to_string(node) : ""), true);
      if ( !converged ) {
        converged = labelPropagationRound(hypergraph, refinement_nodes, start, end, node, i == 0);
      }

      if ( _context.refinement.label_propagation.numa_aware ) {
        // In case, we execute numa-aware label propagation, we need to synchronize
        // the rounds in order that swapping the active node set is thread-safe
        ++_numa_lp_round_synchronization;
        const size_t synchro_barrier = ( i + 1 ) * TBB::instance().num_used_numa_nodes();
        while ( _numa_lp_round_synchronization.load() < synchro_barrier) { }
      }

      ASSERT(static_cast<size_t>(numa_node) < _active.size());
      _active[numa_node].swap(_next_active[numa_node]);
      _next_active[numa_node].reset();
      utils::Timer::instance().stop_timer(
        "lp_round_" + std::to_string(i) + (node != -1 ? "_" + std::to_string(node) : ""));
    }
  }

  bool labelPropagationRound(HyperGraph& hypergraph,
                             parallel::scalable_vector<HypernodeID>& refinement_nodes,
                             const size_t start,
                             const size_t end,
                             const int node,
                             const bool is_first_round) {
    int numa_node = node == -1 ? 0 : node;
    _visited_he[numa_node].reset();
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
    if ( _context.refinement.label_propagation.execute_sequential ) {
      utils::Randomize::instance().localizedShuffleVector(
        refinement_nodes, start, end, _context.shared_memory.shuffle_block_size);
    } else {
      utils::Randomize::instance().localizedParallelShuffleVector(
        refinement_nodes, start, end, _context.shared_memory.shuffle_block_size);
    }

    bool converged = true;
    if ( _context.refinement.label_propagation.execute_sequential ) {
      size_t local_iteration_cnt = 0;
      for ( size_t j = start; j < end; ++j ) {
        const HypernodeID hn = refinement_nodes[j];
        converged &= !moveVertex(hypergraph, hn, local_iteration_cnt,
          node, is_first_round, objective_delta);
      }
    } else {
      tbb::enumerable_thread_specific<size_t> iteration_cnt(0);
      tbb::parallel_for(start, end, [&](const size_t& j) {
          size_t& local_iteration_cnt = iteration_cnt.local();
          const HypernodeID hn = refinement_nodes[j];
          converged &= !moveVertex(hypergraph, hn, local_iteration_cnt,
            node, is_first_round, objective_delta);
        });
    }
    return converged;
  }

  template<typename F>
  bool moveVertex(HyperGraph& hypergraph,
                  const HypernodeID hn,
                  size_t& local_iteration_cnt,
                  const int node,
                  const bool is_first_round,
                  const F& objective_delta) {
    bool is_moved = false;
    if ( hn != kInvalidHypernode &&
         hypergraph.isBorderNode(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));
      int numa_node = node == -1 ? 0 : node;
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(node == -1 || common::get_numa_node_of_vertex(hn) == node);

      // We only compute the max gain move for a node if we are either in the first round of label
      // propagation or if the vertex is still active. A vertex is active, if it changed its block
      // in the last round or one of its neighbors.

      if (is_first_round || _active[0][original_id]) {
        Move best_move = _gain.computeMaxGainMove(hypergraph, hn);
        // We perform a move if it either improves the solution quality or, in case of a
        // zero gain move, the balance of the solution.
        bool perform_move = best_move.gain < 0 ||
                            (_context.refinement.label_propagation.rebalancing &&
                              best_move.gain == 0 &&
                              hypergraph.partWeight(best_move.from) - 1 >
                              hypergraph.partWeight(best_move.to) + 1 &&
                              hypergraph.partWeight(best_move.to) <
                              _context.partition.perfect_balance_part_weights[best_move.to]);
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
                const HyperedgeID original_he_id = hypergraph.originalEdgeID(he);
                if ( !_visited_he[numa_node][original_he_id] ) {
                  for (const HypernodeID& pin : hypergraph.pins(he)) {
                    int pin_numa_node = _context.refinement.label_propagation.numa_aware ?
                      common::get_numa_node_of_vertex(pin) : 0;
                    _next_active[pin_numa_node].set(hypergraph.originalNodeID(pin), true);
                  }
                  _visited_he[numa_node].set(original_he_id, true);
                }
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
          _next_active[0].set(original_id, true);
        }
      }

      ++local_iteration_cnt;
    }

    return is_moved;
  }

  void initializeImpl(HyperGraph& hypergraph) override final {
    const int num_numa_nodes = hypergraph.numNumaHypergraphs();
    _numa_nodes_indices.assign(num_numa_nodes + 1, 0);
    for ( int node = 1; node <= num_numa_nodes; ++node ) {
      _numa_nodes_indices[node] = _numa_nodes_indices[node - 1] + hypergraph.initialNumNodes(node - 1);
    }

    _nodes.assign(hypergraph.initialNumNodes(), kInvalidHypernode);
    auto add_vertex = [&](const HypernodeID hn) {
      const int node = common::get_numa_node_of_vertex(hn);
      const HypernodeID local_id = common::get_local_position_of_vertex(hn);
      ASSERT(_numa_nodes_indices[node] + local_id < _nodes.size());
      _nodes[_numa_nodes_indices[node] + local_id] = hn;
    };

    if ( _context.refinement.label_propagation.execute_sequential ) {
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        add_vertex(hn);
      }
    } else {
      tbb::parallel_for(ID(0), hypergraph.initialNumNodes(), [&](const HypernodeID id) {
        const HypernodeID hn = hypergraph.globalNodeID(id);
        if ( hypergraph.nodeIsEnabled(hn) ) {
          add_vertex(hn);
        }
      });
    }

    _active.clear();
    _next_active.clear();
    _visited_he.clear();
    for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
      _active.emplace_back(hypergraph.initialNumNodes());
      _next_active.emplace_back(hypergraph.initialNumNodes());
      _visited_he.emplace_back(hypergraph.initialNumEdges());
    }
  }

  const Context& _context;
  const TaskGroupID _task_group_id;
  // ! Vector that poins to the first vertex of a consecutive range of nodes
  // ! that belong to the same numa node in _nodes.
  std::vector<size_t> _numa_nodes_indices;
  // ! Contains all nodes that are considered during label propagation.
  // ! The nodes are ordered (increasing) by their numa node which they are assigned to.
  // ! Inside a consecutive range of nodes, that belongs to the same numa node, the
  // ! vertices are pseudo-random shuffled.
  parallel::scalable_vector<HypernodeID> _nodes;
  // ! Computes max gain moves
  GainCalculator _gain;
  // ! Indicate which vertices are active and considered for LP
  NumaFastResetFlagArray _active;
  NumaFastResetFlagArray _next_active;
  NumaFastResetFlagArray _visited_he;
  std::atomic<size_t> _numa_lp_round_synchronization;
};

using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, Km1Policy>;
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, CutPolicy>;
}  // namespace kahypar
