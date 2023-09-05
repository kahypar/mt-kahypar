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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

  template <typename TypeTraits, typename GainTypes>
  bool LabelPropagationRefiner<TypeTraits, GainTypes>::refineImpl(
                  mt_kahypar_partitioned_hypergraph_t& phg,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  const double)  {
    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    resizeDataStructuresForCurrentK();
    _gain.reset();
    _next_active.reset();
    Gain old_quality = best_metrics.quality;

    // Initialize set of active vertices
    initializeActiveNodes(hypergraph, refinement_nodes);

    // Perform Label Propagation
    labelPropagation(hypergraph, best_metrics);

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality ==
      metrics::quality(hypergraph, _context,
        !_context.refinement.label_propagation.execute_sequential),
      V(best_metrics.quality) << V(metrics::quality(hypergraph, _context,
          !_context.refinement.label_propagation.execute_sequential)));

    // Update metrics statistics
    Gain delta = old_quality - best_metrics.quality;
    ASSERT(delta >= 0, "LP refiner worsen solution quality");
    utils::Utilities::instance().getStats(_context.utility_id).update_stat("lp_improvement", delta);
    return delta > 0;
  }


  template <typename TypeTraits, typename GainTypes>
  void LabelPropagationRefiner<TypeTraits, GainTypes>::labelPropagation(PartitionedHypergraph& hypergraph,
                                                                        Metrics& best_metrics) {
    NextActiveNodes next_active_nodes;
    vec<vec<Move>> rebalance_moves_by_part;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i;

      bool stop = false;
      if ( _active_nodes.size() > 0 ) {
        if (_context.refinement.label_propagation.unconstrained) {
          stop = labelPropagationRound<true>(hypergraph, next_active_nodes,
                                             best_metrics, rebalance_moves_by_part);
        } else {
          stop = labelPropagationRound<false>(hypergraph, next_active_nodes,
                                              best_metrics, rebalance_moves_by_part);
        }
      }

      if ( _context.refinement.label_propagation.execute_sequential ) {
        _active_nodes = next_active_nodes.copy_sequential();
        next_active_nodes.clear_sequential();
      } else {
        _active_nodes = next_active_nodes.copy_parallel();
        next_active_nodes.clear_parallel();
      }

      if ( stop || _active_nodes.size() == 0 ) {
        break;
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template<bool unconstrained>
  bool LabelPropagationRefiner<TypeTraits, GainTypes>::labelPropagationRound(
                              PartitionedHypergraph& hypergraph,
                              NextActiveNodes& next_active_nodes,
                              Metrics& best_metrics,
                              vec<vec<Move>>& rebalance_moves_by_part) {
    Metrics current_metrics = best_metrics;
    const bool should_update_gain_cache = !PartitionedHypergraph::is_graph && _gain_cache.isInitialized();
    _visited_he.reset();
    _next_active.reset();
    _gain.reset();
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };

    // Shuffle Vector
    bool converged = true;
    if ( _context.refinement.label_propagation.execute_sequential ) {
      utils::Randomize::instance().shuffleVector(
              _active_nodes, UL(0), _active_nodes.size(), THREAD_ID);

      for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
        const HypernodeID hn = _active_nodes[j];
        if ( moveVertex<unconstrained>(hypergraph, hn, next_active_nodes, objective_delta) ) {
          if (should_update_gain_cache) { _active_node_was_moved[j] = uint8_t(true); }
        } else {
          converged = false;
        }
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
              _active_nodes, UL(0), _active_nodes.size());

      tbb::parallel_for(UL(0), _active_nodes.size(), [&](const size_t& j) {
        const HypernodeID hn = _active_nodes[j];
        if ( moveVertex<unconstrained>(hypergraph, hn, next_active_nodes, objective_delta) ) {
          if (should_update_gain_cache) { _active_node_was_moved[j] = uint8_t(true); }
        } else {
          converged = false;
        }
      });
    }

    current_metrics.imbalance = metrics::imbalance(hypergraph, _context);
    current_metrics.quality += _gain.delta();

    if ( should_update_gain_cache ) {
      auto recompute = [&](size_t j) {
        if ( _active_node_was_moved[j] ) {
          _gain_cache.recomputeInvalidTerms(hypergraph, _active_nodes[j]);
          _active_node_was_moved[j] = uint8_t(false);
        }
      };

      if ( _context.refinement.label_propagation.execute_sequential ) {
        for (size_t j = 0; j < _active_nodes.size(); ++j) {
          recompute(j);
        }
      } else {
        tbb::parallel_for(UL(0), _active_nodes.size(), recompute);
      }
    }

    if constexpr ( unconstrained ) {
      if (!metrics::isBalanced(hypergraph, _context)) {
        converged = applyRebalancing(hypergraph, next_active_nodes, best_metrics,
                                     current_metrics, rebalance_moves_by_part);
      } else {
        updateNodeData(hypergraph, next_active_nodes, false);
      }
    }

    ASSERT(current_metrics.quality <= best_metrics.quality);
    const Gain old_quality = best_metrics.quality;
    best_metrics = current_metrics;

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
    return converged || old_quality - current_metrics.quality <
                        _context.refinement.label_propagation.relative_improvement_threshold * old_quality;
  }

  template <typename TypeTraits, typename GainTypes>
  bool LabelPropagationRefiner<TypeTraits, GainTypes>::applyRebalancing(PartitionedHypergraph& hypergraph,
                                                                        NextActiveNodes& next_active_nodes,
                                                                        Metrics& best_metrics,
                                                                        Metrics& current_metrics,
                                                                        vec<vec<Move>>& rebalance_moves_by_part) {
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    timer.start_timer("rebalance_lp", "Rebalance");
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(hypergraph);
    _rebalancer.refineAndOutputMoves(phg, {}, rebalance_moves_by_part, current_metrics, 0.0);
    timer.stop_timer("rebalance_lp");
    DBG << "[LP] Imbalance after rebalancing: " << current_metrics.imbalance << ", quality: " << current_metrics.quality;

    const bool should_update_gain_cache = !PartitionedHypergraph::is_graph && _gain_cache.isInitialized();
    if (current_metrics.quality > best_metrics.quality) { // rollback and stop LP
      auto noop_obj_fn = [](const SynchronizedEdgeUpdate&) { };
      current_metrics = best_metrics;

      // rollback all changes and update gain cache
      forEachMovedNode(hypergraph, [&](const HypernodeID hn) {
        const PartitionID old_part = _old_part[hn];
        if (hypergraph.partID(hn) != old_part && old_part != kInvalidPartition) {
          changeNodePart<true>(hypergraph, hn, hypergraph.partID(hn), old_part, noop_obj_fn);
          if (should_update_gain_cache) {
            // slightly hacky: use kInvalidPartition to mark nodes that need to be updated
            _old_part[hn] = kInvalidPartition;
          }
        }
      }, &rebalance_moves_by_part);
      if (should_update_gain_cache) {
        forEachMovedNode(hypergraph, [&](const HypernodeID hn) {
          _gain_cache.recomputeInvalidTerms(hypergraph, hn);
        }, &rebalance_moves_by_part);
      }
      return true;
    }

    updateNodeData(hypergraph, next_active_nodes, should_update_gain_cache, &rebalance_moves_by_part);
    return false;
  }

  template <typename TypeTraits, typename GainTypes>
  void LabelPropagationRefiner<TypeTraits, GainTypes>::updateNodeData(PartitionedHypergraph& hypergraph,
                                                                      NextActiveNodes& next_active_nodes,
                                                                      bool should_update_gain_cache,
                                                                      vec<vec<Move>>* rebalance_moves_by_part) {
    // collect activated nodes and update gain cache
    forEachMovedNode(hypergraph, [&](const HypernodeID hn) {
      activateNodeAndNeighbors(hypergraph, next_active_nodes, hn, false);
      if (should_update_gain_cache) {
        _gain_cache.recomputeInvalidTerms(hypergraph, hn);
      }
    }, rebalance_moves_by_part);
    // store current part of each node (required for rollback)
    forEachMovedNode(hypergraph, [&](const HypernodeID hn) {
      _old_part[hn] = hypergraph.partID(hn);
    }, rebalance_moves_by_part);
  }

  template <typename TypeTraits, typename GainTypes>
  template<typename F>
  void LabelPropagationRefiner<TypeTraits, GainTypes>::forEachMovedNode(const PartitionedHypergraph& hypergraph,
                                                                        F node_fn,
                                                                        const vec<vec<Move>>* rebalance_moves_by_part) {
    if (rebalance_moves_by_part != nullptr) {
      if ( _context.refinement.label_propagation.execute_sequential ) {
        for (const auto& moves: *rebalance_moves_by_part) {
          for (const Move& m: moves) {
            node_fn(m.node);
          }
        }
      } else {
        tbb::parallel_for(UL(0), rebalance_moves_by_part->size(), [&](const size_t index) {
          const auto& moves = (*rebalance_moves_by_part)[index];
          tbb::parallel_for(UL(0), moves.size(), [&](const size_t j) {
            node_fn(moves[j].node);
          });
        }, tbb::static_partitioner());
      }
    }

    // NOTE: we need to handle rebalancing nodes and active nodes in separate steps,
    // otherwise a node might be concurrently updated by two threads
    if ( _context.refinement.label_propagation.execute_sequential ) {
      for (const HypernodeID hn: _active_nodes) {
        if (hypergraph.partID(hn) != _old_part[hn]) {
          node_fn(hn);
        }
      }
    } else {
      tbb::parallel_for(UL(0), _active_nodes.size(), [&](const size_t j) {
        const HyperedgeID hn = _active_nodes[j];
        if (hypergraph.partID(hn) != _old_part[hn]) {
          node_fn(hn);
        }
      });
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void LabelPropagationRefiner<TypeTraits, GainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    ActiveNodes tmp_active_nodes;
    _active_nodes = std::move(tmp_active_nodes);

    if ( _context.refinement.label_propagation.execute_sequential ) {
      // Setup active nodes sequential
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        if ( _context.refinement.label_propagation.rebalancing || hypergraph.isBorderNode(hn) ) {
          _active_nodes.push_back(hn);
        }
        if ( _context.refinement.label_propagation.unconstrained ) {
          _old_part[hn] = hypergraph.partID(hn);
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
        if ( _context.refinement.label_propagation.unconstrained ) {
          _old_part[hn] = hypergraph.partID(hn);
        }
      });

      _active_nodes = tmp_active_nodes.copy_parallel();
    }
    _rebalancer.initialize(phg);
  }

  template <typename TypeTraits, typename GainTypes>
  void LabelPropagationRefiner<TypeTraits, GainTypes>::initializeActiveNodes(
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

  namespace {
  #define LABEL_PROPAGATION_REFINER(X, Y) LabelPropagationRefiner<X, Y>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(LABEL_PROPAGATION_REFINER)
}
