/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"


namespace mt_kahypar {

class GlobalRollback {
  static constexpr bool enable_heavy_assert = false;
public:
  explicit GlobalRollback(const Hypergraph& hg, const Context& context) :
          context(context),
          max_part_weight_scaling(context.refinement.fm.rollback_balance_violation_factor),
          num_parts(context.partition.k),
          ets_recalc_data(vec<RecalculationData>(num_parts)),
          last_recalc_round(),
          round(1)
  {
    if (context.refinement.fm.iter_moves_on_recalc && context.refinement.fm.rollback_parallel) {
      last_recalc_round.resize(hg.initialNumEdges(), CAtomic<uint32_t>(0));
    }
  }


  template<bool update_gain_cache>
  HyperedgeWeight revertToBestPrefix(PartitionedHypergraph& phg,
                                     FMSharedData& sharedData,
                                     const vec<HypernodeWeight>& partWeights);

  template<bool update_gain_cache>
  HyperedgeWeight revertToBestPrefixParallel(PartitionedHypergraph& phg,
                                             FMSharedData& sharedData,
                                             const vec<HypernodeWeight>& partWeights,
                                             const std::vector<HypernodeWeight>& maxPartWeights);

  void recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData);

  template<bool update_gain_cache>
  HyperedgeWeight revertToBestPrefixSequential(PartitionedHypergraph& phg,
                                               FMSharedData& sharedData,
                                               const vec<HypernodeWeight>&,
                                               const std::vector<HypernodeWeight>& maxPartWeights);

  template<bool update_gain_cache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void moveVertex(PartitionedHypergraph& phg, HypernodeID u, PartitionID from, PartitionID to) {
    if constexpr (update_gain_cache) {
      phg.changeNodePartWithGainCacheUpdate(u, from, to);
    } else {
      phg.changeNodePart(u, from, to);
    }
  }
  template<bool update_gain_cache>
  bool verifyGains(PartitionedHypergraph& phg, FMSharedData& sharedData);

private:
  const Context& context;

  // ! Factor to multiply max part weight with, in order to relax or disable the balance criterion. Set to zero for disabling
  double max_part_weight_scaling;

  PartitionID num_parts;

  struct RecalculationData {
    MoveID first_in, last_out;
    HypernodeID remaining_pins;
    RecalculationData() :
      first_in(std::numeric_limits<MoveID>::max()),
      last_out(std::numeric_limits<MoveID>::min()),
      remaining_pins(0)
      { }
  };

  tbb::enumerable_thread_specific< vec<RecalculationData> > ets_recalc_data;
  vec<CAtomic<uint32_t>> last_recalc_round;
  uint32_t round;
};

}