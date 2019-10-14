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

#include "tbb/parallel_sort.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/enumerable_thread_specific.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"

namespace mt_kahypar {

template< typename TypeTraits,
          typename ExecutionPolicy,
          template<typename> typename GainPolicy >
class LabelPropagationRefinerT final : public IRefiner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using GainCalculator = GainPolicy<HyperGraph>;

  static constexpr bool debug = false;

 public:

  explicit LabelPropagationRefinerT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _current_level(0),
    _execution_policy(),
    _gain(_hg, _context),
    _active(hypergraph.initialNumNodes(), false) {
    initialize();
  }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT& operator= (const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT& operator= (LabelPropagationRefinerT&&) = delete;

  ~LabelPropagationRefinerT() override = default;

 private:
  bool refineImpl(std::vector<HypernodeID>&,
                  kahypar::Metrics& best_metrics) override final {
    _gain.reset();

    // Label propagation is not executed on all levels of the n-level hierarchy.
    // If LP should be executed on the current level is determined by the execution policy.
    ++_current_level;
    if ( !_execution_policy.execute(_current_level) ) {
      return false;
    }

    if ( _context.refinement.label_propagation.numa_aware ) {
      // Execute label propagation on all numa nodes
      TBB::instance().execute_on_all_numa_nodes([&](const int node) {
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
    ASSERT(delta <= 0, "LP refiner worsen metric");
    ASSERT( current_metric + delta == metrics::objective(_hg, _context.partition.objective),
            V(current_metric) << V(delta) << V(metrics::objective(_hg, _context.partition.objective)) );
    best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);

    return delta < 0;
  }

  void initialize() {
    HypernodeID current_num_nodes = 0;
    for ( const HypernodeID& hn : _hg.nodes() ) {
      unused(hn);
      ++current_num_nodes;
    }
    _execution_policy.initialize(_hg, current_num_nodes);
  }

  void labelPropagation(int node = -1) {

    parallel::scalable_vector<HypernodeID> nodes;
    if ( node == -1 ) {
      // Label propagation is executed on all vertices
      for ( const HypernodeID& hn : _hg.nodes() ) {
        nodes.emplace_back(hn);
        _active[_hg.originalNodeID(hn)] = false;
      }
    } else {
      // Label propagation is executed on all vertices of a NUMA node
      for ( const HypernodeID& hn : _hg.nodes(node) ) {
        nodes.emplace_back(hn);
        _active[_hg.originalNodeID(hn)] = false;
      }
    }


    if ( _context.refinement.label_propagation.use_node_degree_ordering ) {
      // Sort vertices in increasing order of their node degree
      tbb::parallel_sort(nodes.begin(), nodes.end(),
        [&](const HypernodeID& lhs, const HypernodeID& rhs) {
        return _hg.nodeDegree(lhs) < _hg.nodeDegree(rhs);
      });
    } else {
      // Random shuffle vertices
      utils::Randomize::instance().shuffleVector(nodes, nodes.size(), sched_getcpu());
    }

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
    for ( size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations && !converged; ++i ) {
      DBG << "Starting Label Propagation Round" << i << "on NUMA Node" << node;

      converged = true;
      tbb::enumerable_thread_specific<size_t> iteration_cnt(0);
      tbb::parallel_for(tbb::blocked_range<size_t>(0UL, nodes.size()),
        [&](const tbb::blocked_range<size_t>& range) {

        size_t& local_iteration_cnt = iteration_cnt.local();
        for ( size_t j = range.begin(); j < range.end(); ++j ) {
          const HypernodeID hn = nodes[j];
          const HypernodeID original_id = _hg.originalNodeID(hn);

          // We only compute the max gain move for a node if we are either in the first round of label
          // propagation or if the vertex is still active. A vertex is active, if it changed its block
          // in the last round or one of its neighbors.
          if ( i == 0 || _active[original_id] ) {

            Move best_move = _gain.computeMaxGainMove(hn);

            // We perform a move if it either improves the solution quality or, in case of a
            // zero gain move, the balance of the solution.
            bool perform_move =   best_move.gain < 0 ||
                                ( _context.refinement.label_propagation.rebalancing &&
                                  best_move.gain == 0 &&
                                  _hg.localPartWeight(best_move.from) - 1 >
                                  _hg.localPartWeight(best_move.to) + 1 );
            if ( perform_move ) {
              PartitionID from = best_move.from;
              PartitionID to = best_move.to;
              ASSERT(from != to);
              ASSERT(_hg.localPartWeight(to) + _hg.nodeWeight(hn) <= _context.partition.max_part_weights[to]);

              Gain delta_before = _gain.localDelta();
              if ( _hg.changeNodePart(hn, from, to, objective_delta) ) {

                // In case the move to block 'to' was successful, we verify that the "real" gain
                // of the move is either equal to our computed gain or if not, still improves
                // the solution quality.
                Gain move_delta = _gain.localDelta() - delta_before;
                bool accept_move = ( move_delta == best_move.gain || move_delta <= 0 );
                if ( accept_move ) {
                  DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                      << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";

                  // Set all neighbors of the vertex to active
                  for ( const HyperedgeID& he : _hg.incidentEdges(hn) ) {
                    for ( const HypernodeID& pins : _hg.pins(he) ) {
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
            } else {
              _active[original_id] = false;
            }
          }

          ++local_iteration_cnt;
          if ( local_iteration_cnt % _context.refinement.label_propagation.part_weight_update_frequency == 0 ) {
            // We frequently update the local block weights of the current threads
            _hg.updateLocalPartInfos();
          }
        }

      });

    }

  }

  HyperGraph& _hg;
  const Context& _context;
  // ! Incremental counter of calls to LP
  size_t _current_level;
  // ! Indicate whether LP should be executed or not
  ExecutionPolicy _execution_policy;
  // ! Computes max gain moves
  GainCalculator _gain;
  // ! Indicate which vertices are active and considered for LP
  parallel::scalable_vector<bool> _active;
};

template< typename ExecutionPolicy = Mandatory >
using LabelPropagationKm1Refiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, Km1Policy>;
template< typename ExecutionPolicy = Mandatory >
using LabelPropagationCutRefiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy, CutPolicy>;

}  // namespace kahypar
