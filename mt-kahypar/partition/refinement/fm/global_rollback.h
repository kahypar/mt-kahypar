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
  explicit GlobalRollback(const Hypergraph& hg, const Context& context, PartitionID numParts) :
    context(context),
    maxPartWeightScaling(context.refinement.fm.rollback_balance_violation_factor),
    numParts(numParts),
    remaining_original_pins(),
    first_move_in(),
    last_move_out() {

    // In case we perform parallel rollback we need
    // some additional data structures
    if ( context.refinement.fm.revert_parallel ) {
      tbb::parallel_invoke([&] {
        remaining_original_pins.resize("Refinement", "remaining_original_pins",
          static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
      }, [&] {
        first_move_in.resize("Refinement", "first_move_in",
          static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
      }, [&] {
        last_move_out.resize("Refinement", "last_move_out",
          static_cast<size_t>(hg.numNonGraphEdges()) * numParts);
      });
      resetStoredMoveIDs();
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

  MoveID lastMoveOut(HyperedgeID he, PartitionID block) const {
    return last_move_out[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  MoveID firstMoveIn(HyperedgeID he, PartitionID block) const {
    return first_move_in[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  HypernodeID remainingPinsFromBeginningOfMovePhase(HyperedgeID he, PartitionID block) const {
    return remaining_original_pins[size_t(he) * numParts + block].load(std::memory_order_relaxed);
  }

  bool verifyPinCountMatches(const PartitionedHypergraph& phg) const ;

  template<bool update_gain_cache>
  bool verifyGains(PartitionedHypergraph& phg, FMSharedData& sharedData);

  void resetStoredMoveIDs();

  void setRemainingOriginalPins(PartitionedHypergraph& phg);

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

  const Context& context;

  // ! Factor to multiply max part weight with, in order to relax or disable the balance criterion. Set to zero for disabling
  double maxPartWeightScaling;

  PartitionID numParts;

  // ! For each hyperedge and block, the number of pins_in_part
  // ! at the beginning of a move phase minus the number of moved out pins
  // ! A value of zero means at some point a gain improvement might have occured, which then has to be checked
  // ! with the move ID flagging
  ds::Array<CAtomic<HypernodeID>> remaining_original_pins;
  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block
  ds::Array<CAtomic<MoveID>> first_move_in;
  // ! For each hyperedge and each block, the ID of the last move to remove a pin from that block
  ds::Array<CAtomic<MoveID>> last_move_out;
};

}