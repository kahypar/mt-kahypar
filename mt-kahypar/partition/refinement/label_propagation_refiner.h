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

#include "tbb/blocked_range.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"

#include "kahypar/meta/mandatory.h"

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
          typename ExecutionPolicy,
          template <typename> class GainPolicy>
class LabelPropagationRefinerT final : public IRefiner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using GainCalculator = GainPolicy<HyperGraph>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefinerT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _numa_nodes_indices(),
    _nodes(),
    _current_level(0),
    _execution_policy(context.refinement.label_propagation.execution_policy_alpha),
    _gain(_hg, _context),
    _active(hypergraph.initialNumNodes(), false) {
    initialize();
    _nodes.reserve(_hg.initialNumNodes());
  }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT & operator= (const LabelPropagationRefinerT &) = delete;
  LabelPropagationRefinerT & operator= (LabelPropagationRefinerT &&) = delete;

 private:
  bool refineImpl(const std::vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();

    // We collect all contraction partners of the uncontractions in order
    // to avoid rebuilding and reshuffling the set of enabled nodes each
    // time from scratch if we execute label propagation.
    // NOTE, the assumption is that the LP Refiner is called after
    // each uncontraction and that the contraction partner is on the
    // second position of vector refinement_nodes.
    ASSERT(refinement_nodes.size() == 0  /* only for unit tests */ ||
           refinement_nodes.size() == 2);
    if (refinement_nodes.size() == 2) {
      addVertex(refinement_nodes[1]);
    }

    // Label propagation is not executed on all levels of the n-level hierarchy.
    // If LP should be executed on the current level is determined by the execution policy.
    ++_current_level;
    if (!_execution_policy.execute(_current_level)) {
      return false;
    }

    HEAVY_REFINEMENT_ASSERT([&] {
        // Assertion verifies, that all enabled nodes are contained in _nodes
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
        return true;
      } (), "LP Refiner does not contain all vertices hypergraph");

    if (_context.refinement.label_propagation.numa_aware) {
      // Execute label propagation on all numa nodes
      TBB::instance().execute_parallel_on_all_numa_nodes([&](const int node) {
          labelPropagation(node);
        });
    } else {
      // Execute label propagation
      labelPropagation();
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
    return delta < 0;
  }

  void initialize() {
    HypernodeID current_num_nodes = 0;
    for (const HypernodeID& hn : _hg.nodes()) {
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
    int last_node = 0;
    while (last_node < TBB::instance().num_used_numa_nodes() &&
           _hg.initialNumNodes(last_node) == 0) {
      _numa_nodes_indices.push_back(0);
      ++last_node;
    }
    ASSERT(StreamingHyperGraph::get_numa_node_of_vertex(_nodes[0]) == last_node);
    _numa_nodes_indices.push_back(0);
    for (size_t i = 1; i < _nodes.size(); ++i) {
      int node = StreamingHyperGraph::get_numa_node_of_vertex(_nodes[i]);
      if (node != last_node) {
        _numa_nodes_indices.push_back(i);
        last_node = node;
      }
    }
    _numa_nodes_indices.push_back(_nodes.size());

    // Random shuffle all nodes belonging to the same numa node
    ASSERT(_numa_nodes_indices.size() - 1 == (size_t)TBB::instance().num_used_numa_nodes());
    for (int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node) {
      size_t start = _numa_nodes_indices[node];
      size_t end = _numa_nodes_indices[node + 1];
      utils::Randomize::instance().shuffleVector(_nodes, start, end, sched_getcpu());
    }
  }

  void labelPropagation(int node = -1) {
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
                             _gain.computeDeltaForHyperedge(edge_weight, edge_size,
                                                            pin_count_in_from_part_after, pin_count_in_to_part_after);
                           };

    bool converged = false;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations && !converged; ++i) {
      DBG << "Starting Label Propagation Round" << i << "on NUMA Node" << node;

      HighResClockTimepoint start_time = std::chrono::high_resolution_clock::now();
      converged = true;
      tbb::enumerable_thread_specific<size_t> iteration_cnt(0);
      size_t start = node != -1 ? _numa_nodes_indices[node] : 0;
      size_t end = node != -1 ? _numa_nodes_indices[node + 1] : _nodes.size();
      tbb::parallel_for(tbb::blocked_range<size_t>(start, end),
                        [&](const tbb::blocked_range<size_t>& range) {
          size_t& local_iteration_cnt = iteration_cnt.local();
          for (size_t j = range.begin(); j < range.end(); ++j) {
            const HypernodeID hn = _nodes[j];
            const HypernodeID original_id = _hg.originalNodeID(hn);
            ASSERT(node == -1 || StreamingHyperGraph::get_numa_node_of_vertex(hn) == node);

            // We only compute the max gain move for a node if we are either in the first round of label
            // propagation or if the vertex is still active. A vertex is active, if it changed its block
            // in the last round or one of its neighbors.
            if (i == 0 || _active[original_id]) {
              _active[original_id] = false;
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
                      for (const HypernodeID& pins : _hg.pins(he)) {
                        _active[_hg.originalNodeID(pins)] = true;
                      }
                    }
                    converged = false;
                  } else {
                    DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                        << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
                    // In case, the real gain is not equal with the computed gain and
                    // worsen the solution quality we revert the move.
                    ASSERT(_hg.partID(hn) == to);
                    _hg.changeNodePart(hn, to, from, objective_delta);
                  }
                }
                _active[original_id] = true;
              }
            }

            ++local_iteration_cnt;
            if (local_iteration_cnt % _context.refinement.label_propagation.part_weight_update_frequency == 0) {
              // We frequently update the local block weights of the current threads
              _hg.updateLocalPartInfos();
            }
          }
        });
      HighResClockTimepoint end_time = std::chrono::high_resolution_clock::now();
      mt_kahypar::utils::Timer::instance().update_timing(
        "lp_round_" + std::to_string(i) + (node != -1 ? "_" + std::to_string(node) : ""),
        "Label Propagation Round " + std::to_string(i) + (node != -1 ? " - " + std::to_string(node) : ""),
        "label_propagation", mt_kahypar::utils::Timer::Type::REFINEMENT,
        std::chrono::duration<double>(end_time - start_time).count());
    }
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
  parallel::scalable_vector<bool> _active;
};

template <typename ExecutionPolicy = Mandatory>
using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;
template <typename ExecutionPolicy = Mandatory>
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, CutPolicy>;
}  // namespace kahypar
