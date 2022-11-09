/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

#include "mt-kahypar/partition/refinement/fm/strategies/km1_gains.h"
#include "mt-kahypar/utils/reproducible_random.h"

namespace mt_kahypar {

class DeterministicLabelPropagationRefiner final : public IRefiner {
public:
  explicit DeterministicLabelPropagationRefiner(Hypergraph& hypergraph,
                                                const Context& context) :
      context(context),
      compute_gains(context.partition.k),
      moves(hypergraph.initialNumNodes()),
      sorted_moves(hypergraph.initialNumNodes()),
      prng(context.partition.seed),
      active_nodes(0),
      ets_recalc_data( vec<RecalculationData>(context.partition.k) ),
      max_num_nodes(hypergraph.initialNumNodes()),
      max_num_edges(hypergraph.initialNumEdges())
  {
    if (context.refinement.deterministic_refinement.use_active_node_set) {
      active_nodes.adapt_capacity(hypergraph.initialNumNodes());
      last_moved_in_round.resize(hypergraph.initialNumNodes() + hypergraph.initialNumEdges(), CAtomic<uint32_t>(0));
    }

  }

private:
  static constexpr bool debug = false;

  bool refineImpl(PartitionedHypergraph& hypergraph, const vec<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics, double) final ;

  void initializeImpl(PartitionedHypergraph&) final { /* nothing to do */ }

  // functions to apply moves from a sub-round
  Gain applyMovesSortedByGainAndRevertUnbalanced(PartitionedHypergraph& phg);
  Gain applyMovesByMaximalPrefixesInBlockPairs(PartitionedHypergraph& phg);
  Gain applyMovesSortedByGainWithRecalculation(PartitionedHypergraph& phg);
  Gain performMoveWithAttributedGain(PartitionedHypergraph& phg, const Move& m, bool activate_neighbors);
  template<typename Predicate>
  Gain applyMovesIf(PartitionedHypergraph& phg, const vec<Move>& moves, size_t end, Predicate&& predicate);

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void calculateAndSaveBestMove(PartitionedHypergraph& phg, HypernodeID u) {
    assert(u < phg.initialNumNodes());
    if (!phg.isBorderNode(u)) return;
    //auto [to, gain] = compute_gains.local().computeBestTargetBlock(phg, u, context.partition.max_part_weights);
    auto [to, gain] = compute_gains.local().computeBestTargetBlockIgnoringBalance(phg, u);
    if (gain > 0 && to != kInvalidPartition) {    // depending on apply moves function we might do gain >= 0
      moves.push_back_buffered( { phg.partID(u), to, u, gain } );
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void calculateAndSaveBestMoveTwoWay(PartitionedHypergraph& phg, HypernodeID u) {
    if (!phg.isBorderNode(u)) return;
    const Gain gain = TwoWayGainComputer::gainToOtherBlock(phg, u);
    if (gain > 0) {
      moves.push_back_buffered({ phg.partID(u), 1 - phg.partID(u), u, gain });
    }
  }

  struct RecalculationData {
    MoveID first_in, last_out;
    HypernodeID remaining_pins;
    RecalculationData() :
            first_in(std::numeric_limits<MoveID>::max()),
            last_out(std::numeric_limits<MoveID>::min()),
            remaining_pins(0)
    { }
  };

  const Context& context;
  tbb::enumerable_thread_specific<Km1GainComputer> compute_gains;
  ds::BufferedVector<Move> moves;
  vec<Move> sorted_moves;

  std::mt19937 prng;
  utils::ParallelPermutation<HypernodeID> permutation;  // gets memory only once used
  ds::BufferedVector<HypernodeID> active_nodes;
  vec<CAtomic<uint32_t>> last_moved_in_round;
  uint32_t round = 0;


  tbb::enumerable_thread_specific< vec<RecalculationData> > ets_recalc_data;
  vec<CAtomic<uint32_t>> last_recalc_round;
  vec<MoveID> move_pos_of_node;
  uint32_t recalc_round = 1;
  size_t max_num_nodes = 0, max_num_edges = 0;
};

}

