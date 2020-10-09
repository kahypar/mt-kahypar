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

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"


namespace mt_kahypar {
template <template <typename> class GainPolicy>
class LabelPropagationRefiner final : public IRefiner {
 private:
  using GainCalculator = GainPolicy<PartitionedHypergraph>;
  using ActiveNodes = parallel::scalable_vector<HypernodeID>;
  using NextActiveNodes = ds::StreamingVector<HypernodeID>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit LabelPropagationRefiner(Hypergraph& hypergraph,
                                   const Context& context,
                                   const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _current_num_nodes(kInvalidHypernode),
    _current_num_edges(kInvalidHyperedge),
    _gain(context),
    _active_nodes(),
    _active_node_was_moved(hypergraph.initialNumNodes(), uint8_t(false)),
    _next_active(hypergraph.initialNumNodes()),
    _visited_he(hypergraph.initialNumEdges()) { }

  LabelPropagationRefiner(const LabelPropagationRefiner&) = delete;
  LabelPropagationRefiner(LabelPropagationRefiner&&) = delete;

  LabelPropagationRefiner & operator= (const LabelPropagationRefiner &) = delete;
  LabelPropagationRefiner & operator= (LabelPropagationRefiner &&) = delete;

 private:

  bool refineImpl(PartitionedHypergraph& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics,
                  double) final ;

  void labelPropagation(PartitionedHypergraph& hypergraph);

  bool labelPropagationRound(PartitionedHypergraph& hypergraph, NextActiveNodes& next_active_nodes);

  template<typename F>
  bool moveVertex(PartitionedHypergraph& hypergraph,
                  const HypernodeID hn,
                  NextActiveNodes& next_active_nodes,
                  const F& objective_delta) {
    bool is_moved = false;
    ASSERT(hn != kInvalidHypernode);
    if ( hypergraph.isBorderNode(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));

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
        bool changed_part = changeNodePart(hypergraph, hn, from, to, objective_delta);
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
                if ( !_visited_he[he] ) {
                  for (const HypernodeID& pin : hypergraph.pins(he)) {
                    if ( _next_active.compare_and_set_to_true(pin) ) {
                      next_active_nodes.stream(pin);
                    }
                  }
                  _visited_he.set(he, true);
                }
              }
            }
            if ( _next_active.compare_and_set_to_true(hn) ) {
              next_active_nodes.stream(hn);
            }
            is_moved = true;
          } else {
            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
            // In case, the real gain is not equal with the computed gain and
            // worsen the solution quality we revert the move.
            ASSERT(hypergraph.partID(hn) == to);
            changeNodePart(hypergraph, hn, to, from, objective_delta);
          }
        }
      }
    }

    return is_moved;
  }

  void initializeActiveNodes(PartitionedHypergraph& hypergraph,
                             const parallel::scalable_vector<HypernodeID>& refinement_nodes);

  void initializeImpl(PartitionedHypergraph&) final;

  template<typename F>
  bool changeNodePart(PartitionedHypergraph& phg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to,
                      const F& objective_delta) {
    bool success = false;
    if ( _context.partition.paradigm == Paradigm::nlevel && phg.isGainCacheInitialized()) {
      success = phg.changeNodePartWithGainCacheUpdate(hn, from, to,
                                                      _context.partition.max_part_weights[to], [] { }, objective_delta);
    } else {
      success = phg.changeNodePart(hn, from, to,
        _context.partition.max_part_weights[to], []{}, objective_delta);
    }
    return success;
  }

  const Context& _context;
  const TaskGroupID _task_group_id;
  HypernodeID _current_num_nodes;
  HyperedgeID _current_num_edges;
  GainCalculator _gain;
  ActiveNodes _active_nodes;
  parallel::scalable_vector<uint8_t> _active_node_was_moved;
  ds::ThreadSafeFastResetFlagArray<> _next_active;
  kahypar::ds::FastResetFlagArray<> _visited_he;
};

using LabelPropagationKm1Refiner = LabelPropagationRefiner<Km1Policy>;
using LabelPropagationCutRefiner = LabelPropagationRefiner<CutPolicy>;
}  // namespace kahypar
