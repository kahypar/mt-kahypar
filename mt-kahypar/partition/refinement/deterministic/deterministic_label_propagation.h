/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

#include "mt-kahypar/partition/refinement/fm/strategies/km1_gains.h"
#include "mt-kahypar/utils/reproducible_random.h"

namespace mt_kahypar {

class DeterministicLabelPropagationRefiner final : public IRefiner {
public:
  explicit DeterministicLabelPropagationRefiner(Hypergraph& hypergraph, const Context& context, TaskGroupID )
  :
    context(context),
    moves(hypergraph.initialNumNodes()),   // make smaller --> max round size
    sorted_moves(hypergraph.initialNumNodes())
  {

  }
  
private:

  bool refineImpl(PartitionedHypergraph& hypergraph, const vec<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics, double) final ;

  void initializeImpl(PartitionedHypergraph&) final { /* nothing to do */ }


  // can generate a sub-round division
  size_t coloring(PartitionedHypergraph& phg);

  // functions to apply moves from a sub-round
  void applyAllMoves(PartitionedHypergraph& phg);
  void applyMovesByMaximalPrefixesInBlockPairs(PartitionedHypergraph& phg);

  vec<size_t> aggregateDirectionBucketsInplace();

  /*
   * for configs where we don't know exact gains --> have to trace the overall improvement with attributed gains
   * called from applyAllMoves() for example
   */
  Gain performMoveWithAttributedGain(PartitionedHypergraph& phg, const Move& m) {
    Gain attributed_gain = 0;
    auto objective_delta = [&](HyperedgeID he, HyperedgeWeight edge_weight, HypernodeID edge_size,
                               HypernodeID pin_count_in_from_part_after, HypernodeID pin_count_in_to_part_after) {
      attributed_gain -= km1Delta(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
    };
    phg.changeNodePart(m.node, m.from, m.to, objective_delta);
    return attributed_gain;
  }

  /*
   * for configs where we know the exact gains already
   * e.g., if coloring is used or exact gain recalculation (or hybrid)
   */
  void performMove(PartitionedHypergraph& phg, const Move& m) {
    phg.changeNodePart(m.node, m.from, m.to);
  }

  void calculateAndSaveBestMove(PartitionedHypergraph& phg, HypernodeID u) {
    auto [to, gain] = compute_gains.local().computeBestTargetBlockIgnoringBalance(phg, u);
    if (gain > 0 && to != kInvalidPartition) {    // depending on apply moves function we might do gain >= 0
      moves[moves_back.fetch_add(1, std::memory_order_relaxed)] = {phg.partID(u), to};
    }
  }

  const Context& context;
  tbb::enumerable_thread_specific<Km1GainComputer> compute_gains;
  vec<Move> moves, sorted_moves;
  std::atomic<size_t> moves_back = {0};

  std::mt19937 prng;
  utils::ParallelPermutation<HypernodeID> permutation;  // gets memory only once used
  utils::FeistelPermutation feistel_permutation;        // uses barely any memory
};

}

