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

#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"

#include "tbb/parallel_for.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

  template <template <typename> class GainPolicy>
  bool LabelPropagationRefiner<GainPolicy>::refineImpl(
                  PartitionedHypergraph& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  const double)  {
    hypergraph.resetMoveState();
    _gain.reset();
    _next_active.reset();

    // Initialize set of active vertices
    initializeActiveNodes(hypergraph, refinement_nodes);

    // Perform Label Propagation
    labelPropagation(hypergraph);

    // Update global part weight and sizes
    best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

    // Update metrics statistics
    HyperedgeWeight current_metric = best_metrics.getMetric(
    Mode::direct, _context.partition.objective);
    Gain delta = _gain.delta();
    ASSERT(delta <= 0, "LP refiner worsen solution quality");

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
    HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
                            metrics::objective(hypergraph, _context.partition.objective,
                                                !_context.refinement.label_propagation.execute_sequential),
                            V(current_metric) << V(delta) <<
                                              V(metrics::objective(hypergraph, _context.partition.objective,
                                                                    _context.refinement.label_propagation.execute_sequential)));

    best_metrics.updateMetric(current_metric + delta, Mode::direct, _context.partition.objective);
    utils::Utilities::instance().getStats(_context.utility_id).update_stat("lp_improvement", std::abs(delta));
    return delta < 0;
  }


  template <template <typename> class GainPolicy>
  void LabelPropagationRefiner<GainPolicy>::labelPropagation(PartitionedHypergraph& hypergraph) {
    NextActiveNodes next_active_nodes;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i;

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

      if ( _active_nodes.size() == 0 ) {
        break;
      }
    }
  }

  template <template <typename> class GainPolicy>
  bool LabelPropagationRefiner<GainPolicy>::labelPropagationRound(
                              PartitionedHypergraph& hypergraph,
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
        if ( moveVertex(hypergraph, hn, next_active_nodes, objective_delta) ) {
          _active_node_was_moved[j] = uint8_t(true);
        } else {
          converged = false;
        }
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
              _active_nodes, 0UL, _active_nodes.size());

      tbb::parallel_for(0UL, _active_nodes.size(), [&](const size_t& j) {
        const HypernodeID hn = _active_nodes[j];
        if ( moveVertex(hypergraph, hn, next_active_nodes, objective_delta) ) {
          _active_node_was_moved[j] = uint8_t(true);
        } else {
          converged = false;
        }
      });
    }

    if ( ( _context.partition.paradigm == Paradigm::nlevel ||
           _context.refinement.refine_until_no_improvement ) &&
           hypergraph.isGainCacheInitialized() ) {
      auto recompute = [&](size_t j) {
        if ( _active_node_was_moved[j] ) {
          hypergraph.recomputeMoveFromPenalty(_active_nodes[j]);
          _active_node_was_moved[j] = uint8_t(false);
        }
      };

      if ( _context.refinement.label_propagation.execute_sequential ) {
        for (size_t j = 0; j < _active_nodes.size(); ++j) {
          recompute(j);
        }
      } else {
        tbb::parallel_for(0UL, _active_nodes.size(), recompute);
      }
    }

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
    return converged;
  }

  template <template <typename> class GainPolicy>
  void LabelPropagationRefiner<GainPolicy>::initializeImpl(PartitionedHypergraph& hypergraph) {
    ActiveNodes tmp_active_nodes;
    _active_nodes = std::move(tmp_active_nodes);

    if ( _context.refinement.label_propagation.execute_sequential ) {
      // Setup active nodes sequential
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        if ( _context.refinement.label_propagation.rebalancing || hypergraph.isBorderNode(hn) ) {
          _active_nodes.push_back(hn);
        }
      }
    } else {
      // Setup active nodes in parallel
      // A node is active, if it is a border vertex.
      NextActiveNodes tmp_active_nodes;

      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        if ( _context.refinement.label_propagation.rebalancing || hypergraph.isBorderNode(hn) ) {
          tmp_active_nodes.stream(hn);
        }
      });

      _active_nodes = tmp_active_nodes.copy_parallel();
    }
  }

  template <template <typename> class GainPolicy>
  void LabelPropagationRefiner<GainPolicy>::initializeActiveNodes(
                              PartitionedHypergraph& hypergraph,
                              const parallel::scalable_vector<HypernodeID>& refinement_nodes) {
    ActiveNodes tmp_active_nodes;
    _active_nodes = std::move(tmp_active_nodes);

    if ( refinement_nodes.empty() ) {
      if ( _context.refinement.label_propagation.execute_sequential ) {
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
          if ( _next_active.compare_and_set_to_true(hn) ) {
            tmp_active_nodes.stream(hn);
          }
        };

        hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( _context.refinement.label_propagation.rebalancing ||
               hypergraph.isBorderNode(hn) ) {
            add_vertex(hn);
          }
        });

        _active_nodes = tmp_active_nodes.copy_parallel();
      }
    } else {
      _active_nodes = refinement_nodes;
    }

    _next_active.reset();
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  template class LabelPropagationRefiner<Km1Policy>;
  template class LabelPropagationRefiner<CutPolicy>;
}