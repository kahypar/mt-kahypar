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
        best_index = i;
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


class GlobalRollBack {
public:

  explicit GlobalRollBack(size_t numNodes) : gains(numNodes, 0) { }

  HyperedgeWeight globalRollbackToBestPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    recalculateGains(phg, sharedData);

    const MoveID numMoves = sharedData.moveTracker.numPerformedMoves();
    BestIndexReduceBody b(gains);
    // find best index
    tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, numMoves), b, tbb::static_partitioner());

    const auto& move_order = sharedData.moveTracker.globalMoveOrder;
    const PartitionID numParts = sharedData.numParts;

    // revert rejected moves
    tbb::parallel_for(b.best_index + 1, numMoves, [&](const MoveID moveID) {
      const Move& m = move_order[moveID];
      phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{/* do nothing */});
      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        sharedData.remaining_original_pins[e * numParts + m.from].fetch_add(1, std::memory_order_relaxed);
      }
    });

    // apply updates to remaining original pins
    tbb::parallel_for(MoveID(0), b.best_index + 1, [&](const MoveID moveID) {
      const PartitionID to = move_order[moveID].to;
      for (HyperedgeID e : phg.incidentEdges(move_order[moveID].node)) {
        sharedData.remaining_original_pins[e * numParts + to].fetch_add(1, std::memory_order_relaxed);
      }
    });
    assert(sharedData.remaining_original_pins == phg.getPinCountInPartVector());

    if (sharedData.moveTracker.reset()) {
      sharedData.resetStoredMoveIDs();
    }

    return b.best_sum;
  }


  void recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    const auto& move_order = sharedData.moveTracker.globalMoveOrder;
    const MoveID firstMoveID = sharedData.moveTracker.firstMoveID;

    tbb::parallel_for(0U, sharedData.moveTracker.numPerformedMoves(), [&](MoveID localMoveID) {
      MoveID moveID = firstMoveID + localMoveID;
      Gain gain = 0;
      const Move& m = move_order[localMoveID];
      const HypernodeID u = m.node;
      for (HyperedgeID e : phg.incidentEdges(u)) {
        const MoveID firstInFrom = sharedData.firstMoveIn(e, m.from);
        if (sharedData.remainingPinsFromBeginningOfMovePhase(e, m.from) == 0  && sharedData.lastMoveOut(e, m.from) == moveID
            && (firstInFrom > moveID || firstInFrom < firstMoveID)) {
          gain += phg.edgeWeight(e);
        }

        const MoveID lastOutTo = sharedData.lastMoveOut(e, m.to);
        if (sharedData.remainingPinsFromBeginningOfMovePhase(e, m.to) == 0 && sharedData.firstMoveIn(e, m.to) == moveID
            && lastOutTo >= firstMoveID && lastOutTo < moveID) {
          gain -= phg.edgeWeight(e);
        }
      }
      gains[localMoveID] = gain;
    });
  }

  vec<Gain> gains;

};
}
}