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
    compute_gains(context.partition.k),
    moves(hypergraph.initialNumNodes()),   // make smaller --> max round size
    sorted_moves(hypergraph.initialNumNodes())
  {

  }
  
private:
  static constexpr bool debug = false;

  bool refineImpl(PartitionedHypergraph& hypergraph, const vec<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics, double) final ;

  void initializeImpl(PartitionedHypergraph&) final { /* nothing to do */ }


  // functions to apply moves from a sub-round
  Gain applyMovesSortedByGainAndRevertUnbalanced(PartitionedHypergraph& phg);
  Gain applyMovesByMaximalPrefixesInBlockPairs(PartitionedHypergraph& phg);

  vec<size_t> aggregateDirectionBucketsInplace();

  void calculateAndSaveBestMove(PartitionedHypergraph& phg, HypernodeID u) {
    assert(u < phg.initialNumNodes());
    auto [to, gain] = compute_gains.local().computeBestTargetBlock(phg, u, context.partition.max_part_weights);
    // auto [to, gain] = compute_gains.local().computeBestTargetBlockIgnoringBalance(phg, u);
    if (gain > 0 && to != kInvalidPartition) {    // depending on apply moves function we might do gain >= 0
      assert(to >= 0 && to < phg.k());
      size_t pos = moves_back.fetch_add(1, std::memory_order_relaxed);
      assert(pos < moves.size());
      moves[pos] = { phg.partID(u), to, u, gain };
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

