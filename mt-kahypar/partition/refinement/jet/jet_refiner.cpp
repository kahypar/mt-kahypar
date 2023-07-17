/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "mt-kahypar/partition/refinement/jet/jet_refiner.h"

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/datastructures/streaming_vector.h"

namespace mt_kahypar {

  template <typename TypeTraits, typename GainTypes>
  JetRefiner<TypeTraits, GainTypes>::JetRefiner(const HypernodeID num_hypernodes,
                                                const HyperedgeID num_hyperedges,
                                                const Context& context,
                                                GainCache& gain_cache,
                                                IRebalancer& rebalancer) :
    _context(context),
    _gain_cache(gain_cache),
    _precomputed(context.refinement.jet.algorithm == JetAlgorithm::precomputed_ordered),
    _current_k(context.partition.k),
    _top_level_num_nodes(num_hypernodes),
    _current_partition_is_best(true),
    _best_partition(num_hypernodes, kInvalidPartition),
    _current_partition(num_hypernodes, kInvalidPartition),
    _gain(context),
    _active_nodes(),
    _gains_and_target(_precomputed ? num_hypernodes : 0),
    _next_active(_context.refinement.jet.exactly_as_in_jet_paper ? 0 : num_hypernodes),
    _visited_he(_context.refinement.jet.exactly_as_in_jet_paper ? 0 : num_hyperedges),
    _locks(_context.refinement.jet.exactly_as_in_jet_paper ? num_hypernodes : 0),
    _rebalancer(rebalancer) { }

  template <typename TypeTraits, typename GainTypes>
  bool JetRefiner<TypeTraits, GainTypes>::refineImpl(
                  mt_kahypar_partitioned_hypergraph_t& phg,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  const double time_limit)  {
    const HyperedgeWeight input_quality = best_metrics.quality;
    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    Metrics current_metrics = best_metrics;
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    resizeDataStructuresForCurrentK();
    _gain_cache.reset(); // current rebalancer is not capable of using the gain cache

    // Initialize set of active vertices
    timer.start_timer("compute_active_nodes", "Compute Active Nodes");
    if (refinement_nodes.empty()) {
      if (_precomputed) {
        computeActiveNodesFromGraph<true>(hypergraph, true);
      } else {
        computeActiveNodesFromGraph<false>(hypergraph, true);
      }
    } else {
      computeActiveNodesFromVector(hypergraph, refinement_nodes);
    }
    timer.stop_timer("compute_active_nodes");
    DBG << "[JET] initialization done, num_nodes=" << hypergraph.initialNumNodes();

    _current_partition_is_best = true;
    size_t rounds_without_improvement = 0;
    const size_t max_rounds = _context.refinement.jet.fixed_n_iterations;
    for (size_t i = 0; max_rounds == 0 || i < max_rounds; ++i) {
      // We need to store the best partition for rollback and the current partition to determine
      // which nodes have been moved. However, we don't need to write both if the best partition
      // is equal to the current partition.
      if (_current_partition_is_best) {
        storeCurrentPartition(hypergraph, _best_partition);
      } else {
        storeCurrentPartition(hypergraph, _current_partition);
      }
      _gain.reset();

      // Perform Label Propagation
      timer.start_timer("label_propagation_jet", "Label Propagation");
      labelPropagationRound(hypergraph);
      timer.stop_timer("label_propagation_jet");

      // Update metrics statistics
      const HyperedgeWeight old_quality = current_metrics.quality;
      DBG << "[JET] Old imbalance: " << current_metrics.imbalance << ", objective: " << old_quality;
      current_metrics.imbalance = metrics::imbalance(hypergraph, _context);
      Gain delta = _gain.delta();
      DBG << "[JET] New imbalance: " << current_metrics.imbalance << ", tmp delta: " << delta;
      current_metrics.quality += delta;

      // TODO for hypergraphs a recomputePenalties call is needed before rebalancing if rebalancer_v2 is used
      // recomputePenalties(hypergraph, false);

      bool did_rebalance = false;
      if (!metrics::isBalanced(hypergraph, _context)) {
        timer.start_timer("rebalance_jet", "Rebalance");
        rebalance(hypergraph, current_metrics, time_limit, rounds_without_improvement);
        timer.stop_timer("rebalance_jet");

        current_metrics.imbalance = metrics::imbalance(hypergraph, _context);
        delta = current_metrics.quality - old_quality;
        DBG << "[JET] Imbalance after rebalancing: " << current_metrics.imbalance << ", total delta: " << delta;
        did_rebalance = true;
      }
      recomputePenalties(hypergraph, did_rebalance);

      HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
      HEAVY_REFINEMENT_ASSERT(current_metrics.quality ==
        metrics::quality(hypergraph, _context, !_context.refinement.jet.execute_sequential),
        V(current_metrics.quality) << V(metrics::quality(hypergraph, _context,
            !_context.refinement.jet.execute_sequential)));

      ++rounds_without_improvement;
      if (metrics::isBalanced(hypergraph, _context) && (current_metrics.quality < best_metrics.quality
          || (current_metrics.quality == best_metrics.quality && current_metrics.imbalance <= best_metrics.imbalance))) {
        if (best_metrics.quality - current_metrics.quality >
                _context.refinement.jet.relative_improvement_threshold * best_metrics.quality) {
          rounds_without_improvement = 0;
        }
        best_metrics = current_metrics;
        _current_partition_is_best = true;
      } else {
        _current_partition_is_best = false;
      }
      if (rounds_without_improvement >= _context.refinement.jet.num_iterations) {
        break;
      } else {
        // initialize active vertices for next round
        timer.start_timer("compute_active_nodes", "Compute Active Nodes");
        if (_context.refinement.jet.exactly_as_in_jet_paper && _precomputed) {
          computeActiveNodesFromGraph<true>(hypergraph, false);
        } else if (_context.refinement.jet.exactly_as_in_jet_paper) {
          computeActiveNodesFromGraph<false>(hypergraph, false);
        } else if (_precomputed) {
          computeActiveNodesFromPreviousRound<true>(hypergraph);
        } else {
          computeActiveNodesFromPreviousRound<false>(hypergraph);
        }
        timer.stop_timer("compute_active_nodes");
      }
      if (_context.refinement.jet.rollback_after_each_iteration && !_current_partition_is_best) {
        rollbackToBestPartition(hypergraph);
        recomputePenalties(hypergraph, true);
      }
    }

    if (!_current_partition_is_best) {
      DBG << "[JET] Rollback to best partition with value " << best_metrics.quality;
      rollbackToBestPartition(hypergraph);
      recomputePenalties(hypergraph, true);
    }

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality ==
      metrics::quality(hypergraph, _context, !_context.refinement.jet.execute_sequential),
      V(best_metrics.quality) << V(metrics::quality(hypergraph, _context,
          !_context.refinement.jet.execute_sequential)));

    utils::Utilities::instance().getStats(_context.utility_id).update_stat(
      "jet_improvement", input_quality - best_metrics.quality);
    return best_metrics.quality < input_quality;
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::labelPropagationRound(
                  PartitionedHypergraph& hypergraph) {
    _locks.reset();
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

    auto move_node = [&](const size_t j, const bool precomputed) {
      const HypernodeID hn = _active_nodes[j];
      if (precomputed) {
        const PartitionID from = hypergraph.partID(hn);
        const auto [gain, to] = _gains_and_target[hn];
        Gain total_gain = 0;
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          HypernodeID pin_count_in_from_part_after = 0;
          HypernodeID pin_count_in_to_part_after = 1;
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (pin != hn) {
              // Jet uses an order based on the precomputed gain values:
              // If the precomputed gain of another node is better than for the current node
              // (or the gain is equal and the id is smaller), we assume the node is already
              // moved to its target part.
              auto [gain_p, to_p] = _gains_and_target[pin];
              PartitionID part = (gain_p < gain || (gain_p == gain && pin < hn)) ? to_p : hypergraph.partID(pin);
              if (part == from) {
                pin_count_in_from_part_after++;
              } else if (part == to) {
                pin_count_in_to_part_after++;
              }
            }
          }
          total_gain += AttributedGains::gain(he, hypergraph.edgeWeight(he), hypergraph.edgeSize(he),
                                              pin_count_in_from_part_after, pin_count_in_to_part_after);
        }

        if (gain < 0) {
          changeNodePart(hypergraph, hn, from, to, objective_delta);
          if (_context.refinement.jet.exactly_as_in_jet_paper) {
            _locks.set(hn);
          }
        }
      } else {
        bool success = moveVertexGreedily(hypergraph, hn, objective_delta);
        if (success && _context.refinement.jet.exactly_as_in_jet_paper) {
          _locks.set(hn);
        }
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      utils::Randomize::instance().shuffleVector(
              _active_nodes, UL(0), _active_nodes.size(), SCHED_GETCPU);
      // explicit case distinction to enable constant propagation (should be equivalent to template)
      if (_precomputed) {
        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
          move_node(j, true);
        }
      } else {
        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
          move_node(j, false);
        }
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
              _active_nodes, UL(0), _active_nodes.size());
      // explicit case distinction to enable constant propagation (should be equivalent to template)
      if (_precomputed) {
        tbb::parallel_for(UL(0), _active_nodes.size(), [&](size_t j) { move_node(j, true); });
      } else {
        tbb::parallel_for(UL(0), _active_nodes.size(), [&](size_t j) { move_node(j, false); });
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
    _rebalancer.initialize(phg);
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::recomputePenalties(const PartitionedHypergraph& hypergraph,
                                                                          bool did_rebalance) {
    parallel::scalable_vector<PartitionID>& current_parts = _current_partition_is_best ? _best_partition : _current_partition;
    auto recompute = [&](const HypernodeID hn) {
      const bool node_was_moved = (hypergraph.partID(hn) != current_parts[hn]);
      if (node_was_moved) {
        _gain_cache.recomputePenaltyTermEntry(hypergraph, hn);
      } else {
        ASSERT(_gain_cache.penaltyTerm(hn, hypergraph.partID(hn))
              == _gain_cache.recomputePenaltyTerm(hypergraph, hn));
      }
    };

    // TODO this isn't needed for graphs --> skip
    if ( _gain_cache.isInitialized() ) {
      if (did_rebalance) {
        if ( _context.refinement.jet.execute_sequential ) {
          for ( const HypernodeID hn : hypergraph.nodes() ) {
            recompute(hn);
          }
        } else {
          hypergraph.doParallelForAllNodes(recompute);
        }
      } else {
        if ( _context.refinement.jet.execute_sequential ) {
          for (size_t j = 0; j < _active_nodes.size(); ++j) {
            recompute(_active_nodes[j]);
          }
        } else {
          tbb::parallel_for(UL(0), _active_nodes.size(), [&](const size_t j) {
            recompute(_active_nodes[j]);
          });
        }
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template <bool precomputed>
  void JetRefiner<TypeTraits, GainTypes>::computeActiveNodesFromGraph(const PartitionedHypergraph& hypergraph, bool first_round) {
    const bool top_level = (hypergraph.initialNumNodes() == _top_level_num_nodes);
    mt_kahypar::utils::Randomize& randomize = mt_kahypar::utils::Randomize::instance();
    _active_nodes.clear();

    auto process_node = [&](const HypernodeID hn, auto add_node_fn) {
      bool accept_border = !_context.refinement.jet.restrict_to_border_nodes || hypergraph.isBorderNode(hn);
      bool accept_locked = first_round || _context.refinement.jet.vertex_locking == 0.0 || !_locks[hn];
      if (!accept_locked && _context.refinement.jet.vertex_locking < 1.0) {
        accept_locked = randomize.getRandomFloat(0.0, 1.0, SCHED_GETCPU) > _context.refinement.jet.vertex_locking;
      }
      if (accept_border && accept_locked) {
        processNode<precomputed>(hypergraph, hn, add_node_fn, top_level);
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      // setup active nodes sequentially
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        process_node(hn, [&] {
          _active_nodes.push_back(hn);
        });
      }
    } else {
      // setup active nodes in parallel
      ds::StreamingVector<HypernodeID> tmp_active_nodes;
      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        process_node(hn, [&] {
          tmp_active_nodes.stream(hn);
        });
      });
      _active_nodes = tmp_active_nodes.copy_parallel();
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::computeActiveNodesFromVector(const PartitionedHypergraph& hypergraph,
                                                                                    const parallel::scalable_vector<HypernodeID>& refinement_nodes) {
    _active_nodes.clear();

    if (_precomputed) {
      const bool top_level = (hypergraph.initialNumNodes() == _top_level_num_nodes);
      if ( _context.refinement.jet.execute_sequential ) {
        // setup active nodes sequentially
        for ( const HypernodeID hn : refinement_nodes ) {
          processNode<true>(hypergraph, hn, [&]{
            _active_nodes.push_back(hn);
          }, top_level);
        }
      } else {
        // setup active nodes in parallel
        ds::StreamingVector<HypernodeID> tmp_active_nodes;
        tbb::parallel_for(UL(0), refinement_nodes.size(), [&](const size_t j) {
          const HypernodeID hn = refinement_nodes[j];
          processNode<true>(hypergraph, hn, [&]{
            tmp_active_nodes.stream(hn);
          }, top_level);
        });
        _active_nodes = tmp_active_nodes.copy_parallel();
      }
    } else {
      _active_nodes = refinement_nodes;
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template <bool precomputed>
  void JetRefiner<TypeTraits, GainTypes>::computeActiveNodesFromPreviousRound(const PartitionedHypergraph& hypergraph) {
    const bool top_level = (hypergraph.initialNumNodes() == _top_level_num_nodes);
    mt_kahypar::utils::Randomize& randomize = mt_kahypar::utils::Randomize::instance();
    parallel::scalable_vector<PartitionID>& current_parts = _current_partition_is_best ? _best_partition : _current_partition;
    _active_nodes.clear();
    _visited_he.reset();
    _next_active.reset();

    auto activate_neighbors_if_moved = [&](const HypernodeID hn, auto add_node_fn) {
      if (hypergraph.partID(hn) == current_parts[hn]) {
        return; // vertex was not moved
      }

      // Set all neighbors of the vertex to active
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if ( hypergraph.edgeSize(he) <= ID(_context.refinement.jet.hyperedge_size_activation_threshold) ) {
          if ( !_visited_he[he] ) {
            for (const HypernodeID& pin : hypergraph.pins(he)) {
              if ( pin != hn && _next_active.compare_and_set_to_true(pin) ) {
                // check whether vertex locking forbids activating this vertex
                bool accept_locked = (_context.refinement.jet.vertex_locking == 0.0) || hypergraph.partID(pin) == current_parts[pin];
                if (!accept_locked && _context.refinement.jet.vertex_locking < 1.0) {
                  accept_locked = randomize.getRandomFloat(0.0, 1.0, SCHED_GETCPU) > _context.refinement.jet.vertex_locking;
                }
                if (accept_locked) {
                  add_node_fn(pin);
                }
              }
            }
            _visited_he.set(he, true);
          }
        }
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      // setup active nodes sequentially
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        activate_neighbors_if_moved(hn, [&](const HypernodeID neighbor) {
          processNode<precomputed>(hypergraph, neighbor, [&] {
            _active_nodes.push_back(neighbor);
          }, top_level);
        });
      }
    } else {
      // setup active nodes in parallel (use two phases for better load balancing)
      ds::StreamingVector<HypernodeID> tmp_active_nodes;
      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        activate_neighbors_if_moved(hn, [&](const HypernodeID neighbor) {
          tmp_active_nodes.stream(neighbor);
        });
      });
      _active_nodes = tmp_active_nodes.copy_parallel();
      tmp_active_nodes.clear_sequential();
      tbb::parallel_for(UL(0), _active_nodes.size(), [&](const size_t j) {
        processNode<precomputed>(hypergraph, _active_nodes[j], [&] {
          tmp_active_nodes.stream(_active_nodes[j]);
        }, top_level);
      });
      _active_nodes = tmp_active_nodes.copy_parallel();
      tmp_active_nodes.clear_parallel();
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template <bool precomputed, typename F>
  void JetRefiner<TypeTraits, GainTypes>::processNode(const PartitionedHypergraph& hypergraph,
                                                                   const HypernodeID hn, F add_node_fn,
                                                                   const bool top_level) {
    const double gain_factor = top_level ? _context.refinement.jet.negative_gain_factor_fine :
                                           _context.refinement.jet.negative_gain_factor_coarse;
    if constexpr (precomputed) {
      RatingMap& tmp_scores = _gain.localScores();
      Gain isolated_block_gain = 0;
      _gain.precomputeGains(hypergraph, hn, tmp_scores, isolated_block_gain);
      Move best_move = _gain.computeMaxGainMoveForScores(hypergraph, tmp_scores, isolated_block_gain,
                                                          hn, false, false, true);
      tmp_scores.clear();
      bool accept_node = best_move.gain < std::floor(gain_factor * isolated_block_gain);
      if (accept_node) {
        add_node_fn();
        _gains_and_target[hn] = {best_move.gain, best_move.to};
      }
    } else {
      unused(hn);
      add_node_fn();
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::storeCurrentPartition(const PartitionedHypergraph& hypergraph,
                                                                             parallel::scalable_vector<PartitionID>& parts) {
    if ( _context.refinement.jet.execute_sequential ) {
      for (const HypernodeID hn : hypergraph.nodes()) {
        parts[hn] = hypergraph.partID(hn);
      }
    } else {
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        parts[hn] = hypergraph.partID(hn);
      });
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::rollbackToBestPartition(PartitionedHypergraph& hypergraph) {
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

    auto reset_node = [&](const HypernodeID hn) {
      const PartitionID part_id = hypergraph.partID(hn);
      if (part_id != _best_partition[hn]) {
        ASSERT(_best_partition[hn] != kInvalidPartition);
        changeNodePart(hypergraph, hn, part_id, _best_partition[hn], objective_delta);
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      for (const HypernodeID hn : hypergraph.nodes()) {
        reset_node(hn);
      }
    } else {
      hypergraph.doParallelForAllNodes(reset_node);
    }
    _current_partition_is_best = true;
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRefiner<TypeTraits, GainTypes>::rebalance(PartitionedHypergraph& hypergraph,
                                                    Metrics& current_metrics,
                                                    double time_limit,
                                                    size_t& rounds_without_improvement) {
    ASSERT(!_context.partition.deterministic);
    unused(time_limit);
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(hypergraph);
    _rebalancer.jetRebalance(phg, current_metrics, rounds_without_improvement);
  }

  namespace {
  #define JET_REFINER(X, Y) JetRefiner<X, Y>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(JET_REFINER)
}
