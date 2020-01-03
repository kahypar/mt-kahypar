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
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/partition/refinement/policies/execution_policy.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <typename TypeTraits,
          typename ExecutionPolicy,
          template <typename> class GainPolicy>
class LabelPropagationRefinerT final : public IRefiner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using GainCalculator = GainPolicy<HyperGraph>;
  using NumaFasetFlagArray = parallel::scalable_vector<kahypar::ds::FastResetFlagArray<>>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefinerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _numa_nodes_indices(),
    _nodes(),
    _current_level(0),
    _execution_policy(context.refinement.label_propagation.execution_policy_alpha),
    _gain(_hg, _context),
    _active(),
    _next_active(),
    _numa_lp_round_synchronization(0) {
    initialize();
    _nodes.reserve(_hg.initialNumNodes());
    for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
      _active.emplace_back(_hg.initialNumNodes());
      _next_active.emplace_back(_hg.initialNumNodes());
    }
  }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT & operator= (const LabelPropagationRefinerT &) = delete;
  LabelPropagationRefinerT & operator= (LabelPropagationRefinerT &&) = delete;

 private:
  bool refineImpl(const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();

    // We collect all contraction partners of the uncontractions in order
    // to avoid rebuilding and reshuffling the set of enabled nodes each
    // time from scratch if we execute label propagation.
    // NOTE, the assumption is that the LP Refiner is called after
    // each uncontraction and that the contraction partner is on the
    // second position of vector refinement_nodes.
    ASSERT(refinement_nodes.size() % 2 == 0);
    if ( !_context.refinement.label_propagation.localized ) {
      for (size_t i = 0; i < refinement_nodes.size() / 2; ++i) {
        addVertex(refinement_nodes[2 * i + 1]);  // add all contraction partners
        ++_current_level;
      }
    } else {
      _current_level += refinement_nodes.size() / 2;
    }

    // Label propagation is not executed on all levels of the n-level hierarchy.
    // If LP should be executed on the current level is determined by the execution policy.
    if (!_context.refinement.label_propagation.execute_always && !_execution_policy.execute(_current_level)) {
      return false;
    }

    utils::Timer::instance().start_timer("label_propagation", "Label Propagation");

    HEAVY_REFINEMENT_ASSERT([&] {
        // Assertion verifies, that all enabled nodes are contained in _nodes
        if ( !_context.refinement.label_propagation.localized ) {
          std::vector<HypernodeID> tmp_nodes;
          tmp_nodes.insert(tmp_nodes.begin(), _nodes.begin(), _nodes.end());
          for (int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node) {
            size_t pos = _numa_nodes_indices[node];
            size_t end = _numa_nodes_indices[node + 1];
            std::sort(tmp_nodes.begin() + pos, tmp_nodes.begin() + end);
            for (const HypernodeID& hn : _hg.nodes(node)) {
              HypernodeID current_hn = tmp_nodes[pos++];
              if (StreamingHyperGraph::get_numa_node_of_vertex(current_hn) != node) {
                LOG << "Hypernode" << hn << "is not on numa node" << node;
                return false;
              }
              if (current_hn != hn) {
                LOG << "Hypernode" << hn << "not contained in vertex set";
                return false;
              }
            }
          }
        }
        return true;
      } (), "LP Refiner does not contain all vertices hypergraph");

    _numa_lp_round_synchronization = 0;
    for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
      _active[node].reset();
      _next_active[node].reset();
    }

    if ( _context.refinement.label_propagation.localized ) {
      tbb::parallel_for(0UL, refinement_nodes.size(), [&](const size_t i) {
        size_t node = StreamingHyperGraph::get_numa_node_of_vertex(refinement_nodes[i]);
        ASSERT(node < _active.size());
        _active[node].set(_hg.originalNodeID(refinement_nodes[i]), true);
      });
    }

    if (_context.refinement.label_propagation.numa_aware && TBB::instance().num_used_numa_nodes() > 1) {
      // Execute label propagation on all numa nodes
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
          if ( _context.refinement.label_propagation.localized ) {
            // In case, we perform numa-aware localized label propagation, we
            // collect all vertices of the current uncontracted node that belongs
            // to actual numa node and perform label propagation on those nodes.
            parallel::scalable_vector<HypernodeID> current_refinement_nodes;
            for ( const HypernodeID& hn : refinement_nodes ) {
              if ( StreamingHyperGraph::get_numa_node_of_vertex(hn) == node ) {
                current_refinement_nodes.push_back(hn);
              }
            }
            labelPropagation(current_refinement_nodes, 0UL, current_refinement_nodes.size(), node);
          } else {
            // In case we execute numa-aware non-localized label propagation, we perform
            // label propagation on all nodes of the current numa node.
            size_t start = _numa_nodes_indices[node];
            size_t end = _numa_nodes_indices[node + 1];
            labelPropagation(_nodes, start, end, node);
          }
        });
    } else {
      // Execute label propagation
      if ( _context.refinement.label_propagation.localized ) {
        // In case, we execute non-numa-aware localized label propagation, we
        // perform label propagation on the current uncontracted vertices
        labelPropagation(refinement_nodes, 0UL, refinement_nodes.size());
      } else {
        // In case, we execute non-numa-aware non-localized label propagation, we
        // perform label propagation on all vertices
        labelPropagation(_nodes, 0UL, _nodes.size());
      }
    }

    // Update global part weight and sizes
    _hg.updateGlobalPartInfos();
    best_metrics.imbalance = metrics::imbalance(_hg, _context);

    // Update metrics statistics
    HyperedgeWeight current_metric = best_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective);
    Gain delta = _gain.delta();
    ASSERT(delta <= 0, "LP refiner worsen solution quality");
    HEAVY_REFINEMENT_ASSERT(current_metric + delta == metrics::objective(_hg, _context.partition.objective),
                            V(current_metric) << V(delta) << V(metrics::objective(_hg, _context.partition.objective)));
    best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
    utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
    utils::Timer::instance().stop_timer("label_propagation");
    return delta < 0;
  }

  void initialize() {
    HypernodeID current_num_nodes = 0;
    _numa_nodes_indices.assign(TBB::instance().num_used_numa_nodes() + 1, 0);
    for (const HypernodeID& hn : _hg.nodes()) {
      const size_t node = StreamingHyperGraph::get_numa_node_of_vertex(hn);
      ASSERT(node + 1 < _numa_nodes_indices.size());
      ++_numa_nodes_indices[node + 1];
      _nodes.emplace_back(hn);
      ++current_num_nodes;
    }
    _execution_policy.initialize(_hg, current_num_nodes);

    // Sort vertices such that all hypernodes belonging to the same
    // numa node are within a consecutive range in _nodes
    std::sort(_nodes.begin(), _nodes.end(),
              [&](const HypernodeID& lhs, const HypernodeID& rhs) {
        return StreamingHyperGraph::get_numa_node_of_vertex(lhs) <
        StreamingHyperGraph::get_numa_node_of_vertex(rhs);
      });

    // Build up a vector that points to the start position of the vertices
    // belonging to the same numa node in _nodes
    for (size_t i = 1; i < _numa_nodes_indices.size(); ++i) {
      _numa_nodes_indices[i] += _numa_nodes_indices[i - 1];
    }
    ASSERT(_numa_nodes_indices.back() == current_num_nodes);

    // Random shuffle all nodes belonging to the same numa node
    ASSERT(_numa_nodes_indices.size() - 1 == (size_t)TBB::instance().num_used_numa_nodes());
    for (int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node) {
      size_t start = _numa_nodes_indices[node];
      size_t end = _numa_nodes_indices[node + 1];
      utils::Randomize::instance().shuffleVector(_nodes, start, end, sched_getcpu());
    }
  }

  void labelPropagation(const parallel::scalable_vector<HypernodeID>& refinement_nodes,
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
        converged = labelPropagationRound(refinement_nodes, start, end, node, i == 0);
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

  bool labelPropagationRound(const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                             const size_t start,
                             const size_t end,
                             const int node,
                             const bool is_first_round) {
    unused(node);
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
                             _gain.computeDeltaForHyperedge(edge_weight, edge_size,
                                                            pin_count_in_from_part_after, pin_count_in_to_part_after);
                           };

    bool converged = true;
    if ( _context.refinement.label_propagation.execute_sequential ) {
      size_t local_iteration_cnt = 0;
      for ( size_t j = start; j < end; ++j ) {
        const HypernodeID hn = refinement_nodes[j];
        converged &= !moveVertex(hn, local_iteration_cnt, node, is_first_round, objective_delta);
      }
    } else {
      tbb::enumerable_thread_specific<size_t> iteration_cnt(0);
      tbb::parallel_for(start, end, [&](const size_t& j) {
          size_t& local_iteration_cnt = iteration_cnt.local();
          const HypernodeID hn = refinement_nodes[j];
          converged &= !moveVertex(hn, local_iteration_cnt, node, is_first_round, objective_delta);
        });
    }
    return converged;
  }

  template<typename F>
  bool moveVertex(const HypernodeID hn,
                  size_t& local_iteration_cnt,
                  const int node,
                  const bool is_first_round,
                  const F& objective_delta) {
    unused(node);
    bool is_moved = false;
    const HypernodeID original_id = _hg.originalNodeID(hn);
    ASSERT(node == -1 || StreamingHyperGraph::get_numa_node_of_vertex(hn) == node);

    // We only compute the max gain move for a node if we are either in the first round of label
    // propagation or if the vertex is still active. A vertex is active, if it changed its block
    // in the last round or one of its neighbors.
    if (is_first_round || _active[0][original_id]) {
      Move best_move = _gain.computeMaxGainMove(hn);

      // We perform a move if it either improves the solution quality or, in case of a
      // zero gain move, the balance of the solution.
      bool perform_move = best_move.gain < 0 ||
                          (_context.refinement.label_propagation.rebalancing &&
                            best_move.gain == 0 &&
                            _hg.localPartWeight(best_move.from) - 1 >
                            _hg.localPartWeight(best_move.to) + 1);
      if (perform_move) {
        PartitionID from = best_move.from;
        PartitionID to = best_move.to;
        ASSERT(from != to);
        ASSERT(_hg.localPartWeight(to) + _hg.nodeWeight(hn) <= _context.partition.max_part_weights[to]);

        Gain delta_before = _gain.localDelta();
        if (_hg.changeNodePart(hn, from, to, objective_delta)) {
          // In case the move to block 'to' was successful, we verify that the "real" gain
          // of the move is either equal to our computed gain or if not, still improves
          // the solution quality.
          Gain move_delta = _gain.localDelta() - delta_before;
          bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
          if (accept_move) {
            DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";

            // Set all neighbors of the vertex to active
            for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
              for (const HypernodeID& pin : _hg.pins(he)) {
                int pin_numa_node = _context.refinement.label_propagation.numa_aware ?
                  StreamingHyperGraph::get_numa_node_of_vertex(pin) : 0;
                _next_active[pin_numa_node].set(_hg.originalNodeID(pin), true);
              }
            }
            is_moved = true;
          } else {
            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
            // In case, the real gain is not equal with the computed gain and
            // worsen the solution quality we revert the move.
            ASSERT(_hg.partID(hn) == to);
            _hg.changeNodePart(hn, to, from, objective_delta);
          }
        }
        _next_active[0].set(original_id, true);
      }
    }

    ++local_iteration_cnt;
    if (local_iteration_cnt % _context.refinement.label_propagation.part_weight_update_frequency == 0) {
      // We frequently update the local block weights of the current threads
      _hg.updateLocalPartInfos();
    }

    return is_moved;
  }

  // ! Adds a hypernode to set of vertices that are
  // ! considered during label propagation
  void addVertex(const HypernodeID hn) {
    _nodes.emplace_back(hn);
    // Swap new vertex to the consecutive range of hypernodes
    // that belong to the same numa node than the new vertex
    int num_numa_nodes = _numa_nodes_indices.size() - 1;
    int node = StreamingHyperGraph::get_numa_node_of_vertex(hn);
    ++_numa_nodes_indices.back();
    size_t current_position = _nodes.size() - 1;
    int current_numa_node = num_numa_nodes - 1;
    while (node < current_numa_node) {
      size_t& numa_node_start = _numa_nodes_indices[current_numa_node--];
      std::swap(_nodes[current_position], _nodes[numa_node_start]);
      current_position = numa_node_start++;
    }
    // Pseudo-Random Shuffling
    // Since the nodes are initial shuffled (see initialize()) and the contraction order
    // gives us also some kind of randomness, we only swap the new vertex to the middle
    // of the node range of its numa node (avoids expensive call to generate random integers)
    size_t swap_pos = (_numa_nodes_indices[node] + _numa_nodes_indices[node + 1]) / 2;
    std::swap(_nodes[current_position], _nodes[swap_pos]);
  }

  HyperGraph& _hg;
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
  // ! Incremental counter of calls to LP
  size_t _current_level;
  // ! Indicate whether LP should be executed or not
  ExecutionPolicy _execution_policy;
  // ! Computes max gain moves
  GainCalculator _gain;
  // ! Indicate which vertices are active and considered for LP
  NumaFasetFlagArray _active;
  NumaFasetFlagArray _next_active;
  std::atomic<size_t> _numa_lp_round_synchronization;
};

template <typename ExecutionPolicy = Mandatory>
using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;
template <typename ExecutionPolicy = Mandatory>
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, CutPolicy>;
}  // namespace kahypar
