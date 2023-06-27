/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"

namespace mt_kahypar {

template <typename TypeTraits, typename GainTypes>
class RebalancerV2 final : public IRebalancer {
private:
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using GainCalculator = typename GainTypes::GainComputation;
  using AttributedGains = typename GainTypes::AttributedGains;
  using RatingMap = typename GainCalculator::RatingMap;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

public:

  explicit RebalancerV2(const Context& context,
                         GainCache& gain_cache);

  explicit RebalancerV2(const Context& context,
                         gain_cache_t gain_cache);

private:
  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  double);

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) final;

  void setMaxPartWeightsForRoundImpl(const std::vector<HypernodeWeight>& max_part_weights) final;

  bool refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                const vec<HypernodeID>& refinement_nodes,
                                vec<vec<Move>>& moves_by_part,
                                Metrics& best_metrics,
                                const double);

  bool refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                      vec<vec<Move>>* moves_by_part,
                      Metrics& best_metric);


  const Context& _context;
  const HypernodeWeight* _max_part_weights;
  GainCache& _gain_cache;
  PartitionID _current_k;
  GainCalculator _gain;
};

}  // namespace mt_kahypar
