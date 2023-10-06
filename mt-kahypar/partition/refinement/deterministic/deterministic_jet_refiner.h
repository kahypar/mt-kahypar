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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/reproducible_random.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
class DeterministicJetRefiner final : public IRefiner {

  using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
  using GainComputation = typename GraphAndGainTypes::GainComputation;
  using AttributedGains = typename GraphAndGainTypes::AttributedGains;

public:
  explicit DeterministicJetRefiner(const HypernodeID num_hypernodes,
                                                const HyperedgeID num_hyperedges,
                                                const Context& context,
                                                gain_cache_t /* only relevant for other refiners */,
                                                IRebalancer& /* only relevant for other refiners */) :
    DeterministicJetRefiner(num_hypernodes, num_hyperedges, context) { }

  explicit DeterministicJetRefiner(const HypernodeID num_hypernodes,
                                                const HyperedgeID num_hyperedges,
                                                const Context& context) :
      context(context),
      gain_computation(context, true /* disable_randomization */),
      cumulative_node_weights(num_hypernodes),
      moves(num_hypernodes),
      sorted_moves(num_hypernodes),
      current_k(context.partition.k),
      prng(context.partition.seed),
      active_nodes(0) {
    if (context.refinement.deterministic_refinement.use_active_node_set) {
      active_nodes.adapt_capacity(num_hypernodes);
      last_moved_in_round.resize(num_hypernodes + num_hyperedges, CAtomic<uint32_t>(0));
    }
  }

private:
  static constexpr bool debug = false;
  static constexpr size_t invalid_pos = std::numeric_limits<size_t>::max() / 2;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics, double) final ;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final { /* nothing to do */ }

  const Context& context;
  GainComputation gain_computation;
  vec<HypernodeWeight> cumulative_node_weights;
  ds::BufferedVector<Move> moves;
  vec<Move> sorted_moves;

  PartitionID current_k;
  std::mt19937 prng;
  utils::ParallelPermutation<HypernodeID> permutation;
  ds::BufferedVector<HypernodeID> active_nodes;
  vec<CAtomic<uint32_t>> last_moved_in_round;
  uint32_t round = 0;
};

}
