/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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