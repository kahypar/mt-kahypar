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
 * The above copyright notice and this permission notice shall be included in
 *all copies or substantial portions of the Software.
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

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {
template <typename TypeTraits, typename GainTypes>
class KaminparRebalancer final : public IRebalancer {
 private:
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using GainCalculator = typename GainTypes::GainComputation;
  using AttributedGains = typename GainTypes::AttributedGains;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit KaminparRebalancer(const Context& context, GainCache& gain_cache)
      : _context(context), _gain_cache(gain_cache), _max_weight(0) { }

  explicit KaminparRebalancer(const Context& context, gain_cache_t gain_cache)
      : KaminparRebalancer(context,
                    GainCachePtr::cast<GainCache>(gain_cache)) {}

  KaminparRebalancer(const KaminparRebalancer&) = delete;
  KaminparRebalancer(KaminparRebalancer&&) = delete;

  KaminparRebalancer& operator=(const KaminparRebalancer&) = delete;
  KaminparRebalancer& operator=(KaminparRebalancer&&) = delete;

 private:
  bool refineImpl(
      mt_kahypar_partitioned_hypergraph_t& hypergraph,
      const parallel::scalable_vector<HypernodeID>& refinement_nodes,
      Metrics& best_metrics, double) final;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final { };

  void setMaxPartWeightsForRoundImpl(const std::vector<HypernodeWeight>& max_part_weights) final {
    ASSERT(max_part_weights[0] == max_part_weights[1]);
    _max_weight = max_part_weights[0] * _context.partition.k;
  }

  bool refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t&,
                                const vec<HypernodeID>&,
                                vec<vec<Move>>&,
                                Metrics&,
                                const double) {
    ALWAYS_ASSERT(false, "not implemented");
  }

  bool isBalanced(const PartitionedHypergraph& phg) {
    ASSERT(_max_weight != 0);
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (phg.partWeight(i) * _context.partition.k > _max_weight) {
        return false;
      }
    }
    return true;
  }

  const Context& _context;
  GainCache& _gain_cache;
  HypernodeWeight _max_weight;
};

}  // namespace mt_kahypar
