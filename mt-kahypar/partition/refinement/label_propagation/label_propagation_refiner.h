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

#pragma once

#include "kahypar-resources/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"


namespace mt_kahypar {
template <typename TypeTraits, typename GainTypes>
class LabelPropagationRefiner final : public IRefiner {
 private:
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using GainCalculator = typename GainTypes::GainComputation;
  using ActiveNodes = parallel::scalable_vector<HypernodeID>;
  using NextActiveNodes = ds::StreamingVector<HypernodeID>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefiner(const HypernodeID num_hypernodes,
                                   const HyperedgeID num_hyperedges,
                                   const Context& context,
                                   GainCache& gain_cache,
                                   IRebalancer& rb) :
    _context(context),
    _gain_cache(gain_cache),
    _current_k(context.partition.k),
    _current_num_nodes(kInvalidHypernode),
    _current_num_edges(kInvalidHyperedge),
    _gain(context),
    _active_nodes(),
    _active_node_was_moved(num_hypernodes, uint8_t(false)),
    _old_part(_context.refinement.label_propagation.unconstrained ? num_hypernodes : 0, kInvalidPartition),
    _next_active(num_hypernodes),
    _visited_he(Hypergraph::is_graph ? 0 : num_hyperedges),
    _rebalancer(rb) { }

  explicit LabelPropagationRefiner(const HypernodeID num_hypernodes,
                                   const HyperedgeID num_hyperedges,
                                   const Context& context,
                                   gain_cache_t gain_cache,
                                   IRebalancer& rb) :
    LabelPropagationRefiner(num_hypernodes, num_hyperedges, context,
      GainCachePtr::cast<GainCache>(gain_cache), rb) { }

  LabelPropagationRefiner(const LabelPropagationRefiner&) = delete;
  LabelPropagationRefiner(LabelPropagationRefiner&&) = delete;

  LabelPropagationRefiner & operator= (const LabelPropagationRefiner &) = delete;
  LabelPropagationRefiner & operator= (LabelPropagationRefiner &&) = delete;

 private:
  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  double) final ;

  void labelPropagation(PartitionedHypergraph& phg, Metrics& best_metrics);

  template<bool unconstrained>
  bool labelPropagationRound(PartitionedHypergraph& hypergraph,
                             NextActiveNodes& next_active_nodes,
                             Metrics& best_metrics,
                             vec<vec<Move>>& rebalance_moves_by_part);

  bool applyRebalancing(PartitionedHypergraph& hypergraph,
                        NextActiveNodes& next_active_nodes,
                        Metrics& best_metrics,
                        Metrics& current_metrics,
                        vec<vec<Move>>& rebalance_moves_by_part);

  void updateNodeData(PartitionedHypergraph& hypergraph,
                      NextActiveNodes& next_active_nodes,
                      vec<vec<Move>>* rebalance_moves_by_part = nullptr);

  template<bool unconstrained, typename F>
  bool moveVertex(PartitionedHypergraph& hypergraph,
                  const HypernodeID hn,
                  NextActiveNodes& next_active_nodes,
                  const F& objective_delta) {
    bool is_moved = false;
    ASSERT(hn != kInvalidHypernode);
    if ( hypergraph.isBorderNode(hn) && !hypergraph.isFixed(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));

      Move best_move = _gain.computeMaxGainMove(hypergraph, hn, false, false, unconstrained);
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
        bool changed_part = changeNodePart<unconstrained>(hypergraph, hn, from, to, objective_delta);
        ASSERT(!unconstrained || changed_part);
        is_moved = true;
        if (unconstrained || changed_part) {
          // In case the move to block 'to' was successful, we verify that the "real" gain
          // of the move is either equal to our computed gain or if not, still improves
          // the solution quality.
          Gain move_delta = _gain.localDelta() - delta_before;
          bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
          if (accept_move) {
            DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";
            if constexpr (!unconstrained) {
              // in unconstrained case, we don't want to activate neighbors if the move is undone
              // by the rebalancing
              activateNodeAndNeighbors(hypergraph, next_active_nodes, hn, true);
            }
          } else {
            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
            // In case, the real gain is not equal with the computed gain and
            // worsen the solution quality we revert the move.
            ASSERT(hypergraph.partID(hn) == to);
            changeNodePart<unconstrained>(hypergraph, hn, to, from, objective_delta);
          }
        }
      }
    }

    return is_moved;
  }

  void initializeActiveNodes(PartitionedHypergraph& hypergraph,
                             const parallel::scalable_vector<HypernodeID>& refinement_nodes);

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final;

  template<bool unconstrained, typename F>
  bool changeNodePart(PartitionedHypergraph& phg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to,
                      const F& objective_delta) {
    HypernodeWeight max_weight = unconstrained ? std::numeric_limits<HypernodeWeight>::max()
                                                 : _context.partition.max_part_weights[to];
    if ( _gain_cache.isInitialized() ) {
      return phg.changeNodePart(_gain_cache, hn, from, to, max_weight, []{}, objective_delta);
    } else {
      return phg.changeNodePart(hn, from, to, max_weight, []{}, objective_delta);
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void activateNodeAndNeighbors(PartitionedHypergraph& hypergraph,
                                NextActiveNodes& next_active_nodes,
                                const HypernodeID hn,
                                bool activate_moved) {
    auto activate = [&](const HypernodeID hn) {
      if (activate_moved || hypergraph.partID(hn) == _old_part[hn]) {
        if ( _next_active.compare_and_set_to_true(hn) ) {
          next_active_nodes.stream(hn);
        }
      }
    };

    // Set all neighbors of the vertex to active
    if constexpr (Hypergraph::is_graph) {
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        activate(hypergraph.edgeTarget(he));
      }
    } else {
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if ( hypergraph.edgeSize(he) <=
              ID(_context.refinement.label_propagation.hyperedge_size_activation_threshold) ) {
          if ( !_visited_he[he] ) {
            for (const HypernodeID& pin : hypergraph.pins(he)) {
              activate(pin);
            }
            _visited_he.set(he, true);
          }
        }
      }
    }

    if ( activate_moved && _next_active.compare_and_set_to_true(hn) ) {
      next_active_nodes.stream(hn);
    }
  }

  void resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if ( _current_k != _context.partition.k ) {
      _current_k = _context.partition.k;
      _gain.changeNumberOfBlocks(_current_k);
    }
  }

  const Context& _context;
  GainCache& _gain_cache;
  PartitionID _current_k;
  HypernodeID _current_num_nodes;
  HyperedgeID _current_num_edges;
  GainCalculator _gain;
  ActiveNodes _active_nodes;
  parallel::scalable_vector<uint8_t> _active_node_was_moved;
  parallel::scalable_vector<PartitionID> _old_part;
  ds::ThreadSafeFastResetFlagArray<> _next_active;
  kahypar::ds::FastResetFlagArray<> _visited_he;
  IRebalancer& _rebalancer;
};

}  // namespace kahypar
