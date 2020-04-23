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

#include <atomic>
#include <mt-kahypar/partition/metrics.h>

#include "fm_commons.h"

namespace mt_kahypar {
namespace refinement {

struct BestIndexReduceBody {
  const vec<Gain>& gains;
  MoveID best_index = std::numeric_limits<MoveID>::max();
  HyperedgeWeight sum = 0, best_sum = 0;

  BestIndexReduceBody(const vec<Gain>& gains) : gains(gains) { }

  BestIndexReduceBody(BestIndexReduceBody& b, tbb::split) : gains(b.gains) { }

  void operator()(const tbb::blocked_range<MoveID>& r) {
    if (best_index == std::numeric_limits<MoveID>::max()) {
      best_index = r.begin(); // initialize if the range was split
    }
    for (MoveID i = r.begin(); i < r.end(); ++i) {
      sum += gains[i];
      if (sum > best_sum) {
        best_sum = sum;
        best_index = i + 1;
      }
    }
  }

  void join(BestIndexReduceBody& rhs) {
    const HyperedgeWeight rhs_sum = sum + rhs.sum;
    const HyperedgeWeight rhs_best_sum = sum + rhs.best_sum;
    if (rhs_best_sum > best_sum) {
      best_index = rhs.best_index;
      best_sum = rhs_best_sum;
    }
    sum = rhs_sum;
  }

};


class GlobalRollback {
  static constexpr bool enable_heavy_assert = false;
public:

  explicit GlobalRollback(size_t numNodes, size_t numHyperedges, PartitionID numParts) :
          numParts(numParts),
          remaining_original_pins(numHyperedges * numParts),
          first_move_in(numHyperedges * numParts),
          last_move_out(numHyperedges * numParts),
          gains(numNodes, 0)
  {
    resetStoredMoveIDs();
  }

  HyperedgeWeight globalRollbackToBestPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData) {

    recalculateGains(phg, sharedData);

    const MoveID numMoves = sharedData.moveTracker.numPerformedMoves();
    BestIndexReduceBody b(gains);
    tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, numMoves), b, tbb::static_partitioner()); // find best index

    const auto& move_order = sharedData.moveTracker.moveOrder;

    // revert rejected moves
    tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
      const Move& m = move_order[moveID];
      phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{/* do nothing */});
      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        remaining_original_pins[e * numParts + m.from].fetch_add(1, std::memory_order_relaxed);
      }
    });

    // apply updates to remaining original pins
    tbb::parallel_for(MoveID(0), b.best_index, [&](const MoveID moveID) {
      const PartitionID to = move_order[moveID].to;
      for (HyperedgeID e : phg.incidentEdges(move_order[moveID].node)) {
        remaining_original_pins[e * numParts + to].fetch_add(1, std::memory_order_relaxed);
      }
    });

    // recompute moveFromBenefit values since they are potentially invalid
    tbb::parallel_for(0U, sharedData.moveTracker.numPerformedMoves(), [&](MoveID localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    });

    if (sharedData.moveTracker.reset()) {
      resetStoredMoveIDs();
    }

    HEAVY_REFINEMENT_ASSERT(std::equal(phg.getPinCountInPartVector().begin(), phg.getPinCountInPartVector().end(), remaining_original_pins.begin()));
    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
    return b.best_sum;
  }


  void recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    const vec<Move>& move_order = sharedData.moveTracker.moveOrder;
    MoveID firstMoveID = sharedData.moveTracker.firstMoveID;

#ifndef NDEBUG
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }
    phg.checkTrackedPartitionInformation();
#endif


    tbb::parallel_for(0U, sharedData.moveTracker.numPerformedMoves(), [&](MoveID localMoveID) {
      const Move& m = move_order[localMoveID];
      const MoveID globalMoveID = localMoveID + firstMoveID;

      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        CAtomic<MoveID>& fmi = first_move_in[e * numParts + m.to];
        MoveID expected = fmi.load(std::memory_order_acq_rel);
        while ((sharedData.moveTracker.isIDStale(expected) || expected > globalMoveID)
               && !fmi.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

        CAtomic<MoveID>& lmo = last_move_out[e * numParts + m.from];
        expected = lmo.load(std::memory_order_acq_rel);
        while (expected < globalMoveID && !lmo.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

        remaining_original_pins[e * numParts + m.from].fetch_sub(1, std::memory_order_relaxed);
      }
    });

    tbb::parallel_for(0U, sharedData.moveTracker.numPerformedMoves(), [&](MoveID localMoveID) {
      MoveID globalMoveID = firstMoveID + localMoveID;
      Gain gain = 0;
      const Move& m = move_order[localMoveID];
      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        const MoveID firstInFrom = firstMoveIn(e, m.from);
        if (remainingPinsFromBeginningOfMovePhase(e, m.from) == 0  && lastMoveOut(e, m.from) == globalMoveID
            && (firstInFrom > globalMoveID || firstInFrom < firstMoveID)) {
          gain += phg.edgeWeight(e);
        }

        const MoveID lastOutTo = lastMoveOut(e, m.to);
        if (remainingPinsFromBeginningOfMovePhase(e, m.to) == 0 && firstMoveIn(e, m.to) == globalMoveID && lastOutTo < globalMoveID) {
          gain -= phg.edgeWeight(e);
        }
      }
      gains[localMoveID] = gain;
    });


#ifndef NDEBUG
    // recheck all gains
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{});
    }

    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }

    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      const Gain estimated_gain = phg.km1Gain(m.node, m.from, m.to);
      assert(phg.moveFromBenefit(m.node) == phg.moveFromBenefitRecomputed(m.node));
      assert(phg.moveToPenalty(m.node, m.to) == phg.moveToPenaltyRecomputed(m.node, m.to));
      if (tbb::this_task_arena::max_concurrency() == 1) {
        assert(estimated_gain == m.gain);
      }
      const HyperedgeWeight km1_before_move = metrics::km1(phg, false);
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(), []{});
      const HyperedgeWeight km1_after_move = metrics::km1(phg, false);
      assert(km1_after_move + estimated_gain == km1_before_move);
      assert(km1_after_move + gains[localMoveID] == km1_before_move);
      assert(estimated_gain == gains[localMoveID]);
    }

    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }
#endif
  }




  MoveID lastMoveOut(HyperedgeID he, PartitionID block) const {
    return last_move_out[he * numParts + block].load(std::memory_order_relaxed);
  }

  MoveID firstMoveIn(HyperedgeID he, PartitionID block) const {
    return first_move_in[he * numParts + block].load(std::memory_order_relaxed);
  }


  HypernodeID remainingPinsFromBeginningOfMovePhase(HyperedgeID he, PartitionID block) const {
    return remaining_original_pins[he * numParts + block].load(std::memory_order_relaxed);
  }

  void resetStoredMoveIDs() {
    for (auto &x : last_move_out)
      x.store(0, std::memory_order_relaxed);
    for (auto &x : first_move_in)
      x.store(0, std::memory_order_relaxed);
  }

  void setRemainingOriginalPins(PartitionedHypergraph& phg) {
    assert(remaining_original_pins.size() >= phg.getPinCountInPartVector().size());
    size_t n = phg.getPinCountInPartVector().size();
    std::copy_n(phg.getPinCountInPartVector().begin(), n, remaining_original_pins.begin());
  }

  PartitionID numParts;

  // ! For each hyperedge and block, the number of pins_in_part at the beginning of a move phase minus the number of moved out pins
  vec<CAtomic<HypernodeID>> remaining_original_pins;
  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block / the last move to remove a pin from that block
  vec<CAtomic<MoveID>> first_move_in, last_move_out;

  vec<Gain> gains;

};
}
}