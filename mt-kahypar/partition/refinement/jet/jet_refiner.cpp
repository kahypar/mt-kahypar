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

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  bool JetRefiner<TypeTraits, GainTypes, precomputed>::refineImpl(
                  mt_kahypar_partitioned_hypergraph_t& phg,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  const double)  {
    PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
    resizeDataStructuresForCurrentK();
    _gain.reset();

    // Initialize set of active vertices
    initializeActiveNodes(hypergraph, refinement_nodes);

    // Perform Label Propagation
    labelPropagationRound(hypergraph);


    // Update global part weight and sizes
    DBG << "[JET] Old imbalance: " << best_metrics.imbalance;
    best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

    // Update metrics statistics
    Gain delta = _gain.delta();
    DBG << "[JET] New imbalance: " << best_metrics.imbalance << ", delta: " << delta;

    recomputePenalties(hypergraph);

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta ==
      metrics::quality(hypergraph, _context,
        !_context.refinement.label_propagation.execute_sequential),
      V(best_metrics.quality) << V(delta) << V((best_metrics.quality + delta))
        << V(metrics::quality(hypergraph, _context,
          !_context.refinement.label_propagation.execute_sequential)));

    best_metrics.quality += delta;
    utils::Utilities::instance().getStats(_context.utility_id).update_stat("jet_improvement", std::abs(delta));
    return delta < 0;
  }

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  void JetRefiner<TypeTraits, GainTypes, precomputed>::labelPropagationRound(
                  PartitionedHypergraph& hypergraph) {
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

    auto move_node = [&](const size_t j) {
      const HypernodeID hn = _active_nodes[j];
      if constexpr (precomputed) {
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
          _active_node_was_moved[hn] = uint8_t(true);
        }
      } else {
        if ( moveVertexGreedily(hypergraph, hn, objective_delta) ) {
          _active_node_was_moved[hn] = uint8_t(true);
        }
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
        move_node(j);
      }
    } else {
      tbb::parallel_for(UL(0), _active_nodes.size(), move_node);
    }

    HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation(_gain_cache));
  }

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  void JetRefiner<TypeTraits, GainTypes, precomputed>::rollback(
                  PartitionedHypergraph& hypergraph) {
    auto reset_node = [&](const HypernodeID hn) {
      const PartitionID part_id = hypergraph.partID(hn);
      if (part_id != _old_parts[hn]) {
        bool success = false;
        if ( _context.forceGainCacheUpdates() && _gain_cache.isInitialized() ) {
          success = hypergraph.changeNodePart(_gain_cache, hn, part_id, _old_parts[hn]);
        } else {
          success = hypergraph.changeNodePart(hn, part_id, _old_parts[hn]);
        }
        ASSERT(success);
      }
    };

    if ( _context.refinement.jet.execute_sequential ) {
      for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
        reset_node(j);
      }
    } else {
      tbb::parallel_for(UL(0), _active_nodes.size(), reset_node);
    }
    // TODO: update gain cache?!
  }

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  void JetRefiner<TypeTraits, GainTypes, precomputed>::initializeImpl(mt_kahypar_partitioned_hypergraph_t&) {
    // PartitionedHypergraph& hypergraph = utils::cast<PartitionedHypergraph>(phg);
  }

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  void JetRefiner<TypeTraits, GainTypes, precomputed>::recomputePenalties(const PartitionedHypergraph& hypergraph) {
    if ( _context.forceGainCacheUpdates() && _gain_cache.isInitialized() ) {
      auto recompute = [&](size_t j) {
        const HypernodeID hn = _active_nodes[j];
        if ( _active_node_was_moved[hn] ) {
          _gain_cache.recomputePenaltyTermEntry(hypergraph, hn);
          if (!_context.refinement.jet.vertex_locking) {
            _active_node_was_moved[hn] = uint8_t(false);
          }
        }
      };

      if ( _context.refinement.jet.execute_sequential ) {
        for (size_t j = 0; j < _active_nodes.size(); ++j) {
          recompute(j);
        }
      } else {
        tbb::parallel_for(UL(0), _active_nodes.size(), recompute);
      }
    }
  }

  template <typename TypeTraits, typename GainTypes, bool precomputed>
  void JetRefiner<TypeTraits, GainTypes, precomputed>::initializeActiveNodes(
                              PartitionedHypergraph& hypergraph,
                              const parallel::scalable_vector<HypernodeID>& refinement_nodes) {
    ASSERT(_active_node_was_moved.size() == hypergraph.initialNumNodes());
    _active_nodes.clear();
    _gains_and_target.clear();

    auto process_node = [&](const HypernodeID hn, auto add_node_fn) {
      bool accept_border = !_context.refinement.jet.restrict_to_border_nodes || hypergraph.isBorderNode(hn);
      bool accept_locked = !_context.refinement.jet.vertex_locking || !_active_node_was_moved[hn];
      if ( accept_border && accept_locked ) {
        const PartitionID from = hypergraph.partID(hn);
        if constexpr (precomputed) {
          RatingMap& tmp_scores = _gain.localScores();
          Gain isolated_block_gain = 0;
          _gain.precomputeGains(hypergraph, hn, tmp_scores, isolated_block_gain);
          Move best_move = _gain.computeMaxGainMoveForScores(hypergraph, tmp_scores, isolated_block_gain,
                                                             hn, false, false, true);
          tmp_scores.clear();
          // TODO: fine factor?
          bool accept_node = best_move.gain <  std::floor(
                _context.refinement.jet.negative_gain_factor_coarse * isolated_block_gain);
          if (accept_node) {
            add_node_fn(hn);
            _gains_and_target[hn] = {best_move.gain, best_move.to};
          }
        } else {
          add_node_fn(hn);
        }
        _old_parts[hn] = from;
      } else if (!accept_locked) {
        ASSERT(_context.refinement.jet.vertex_locking);
        _active_node_was_moved[hn] = false;
      }
    };

    if ( refinement_nodes.empty() ) {
        // setup active nodes sequentially
      if ( _context.refinement.jet.execute_sequential ) {
        for ( const HypernodeID hn : hypergraph.nodes() ) {
          process_node(hn, [&](const HypernodeID hn) {
            _active_nodes.push_back(hn);
          });
        }
      } else {
        // setup active nodes in parallel
        ds::StreamingVector<HypernodeID> tmp_active_nodes;
        hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
          process_node(hn, [&](const HypernodeID hn) {
            tmp_active_nodes.stream(hn);
          });
        });

        _active_nodes = tmp_active_nodes.copy_parallel();
      }
    } else {
      ALWAYS_ASSERT(false); // TODO: rollback
      if constexpr (precomputed) {
        ALWAYS_ASSERT(false); // TODO
      } else {
        _active_nodes = refinement_nodes;
      }
    }

    if ( _context.refinement.jet.execute_sequential ) {
      utils::Randomize::instance().shuffleVector(
              _active_nodes, UL(0), _active_nodes.size(), SCHED_GETCPU);
    } else {
      utils::Randomize::instance().parallelShuffleVector(
              _active_nodes, UL(0), _active_nodes.size());
    }
  }

  namespace {
  #define JET_REFINER_GREEDY(X, Y) JetRefiner<X, Y, false>
  #define JET_REFINER_PRE(X, Y) JetRefiner<X, Y, true>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(JET_REFINER_GREEDY)
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(JET_REFINER_PRE)
}
