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

#include <mt-kahypar/definitions.h>

namespace mt_kahypar {
namespace refinement {


struct GlobalMoveTracker {
  vec<Move> globalMoveOrder;
  parallel::IntegralAtomicWrapper<MoveID> runningMoveID;
  MoveID firstMoveID = 1;

  explicit GlobalMoveTracker(size_t numNodes) : globalMoveOrder(numNodes), runningMoveID(1) {}

  // Returns true if stored move IDs should be reset
  bool reset() {
    if (runningMoveID.load() >= std::numeric_limits<MoveID>::max() - globalMoveOrder.size() - 20) {
      firstMoveID = 1;
      runningMoveID.store(1);
      return true;
    } else {
      firstMoveID = ++runningMoveID;
      return false;
    }
  }

  MoveID insertMove(Move &m) {
    const MoveID move_id = runningMoveID.fetch_add(1, std::memory_order_relaxed);
    globalMoveOrder[move_id - firstMoveID] = m;
    return move_id;
  }

  MoveID numPerformedMoves() const {
    return runningMoveID.load(std::memory_order_relaxed) - firstMoveID;
  }

  bool isIDStale(const MoveID move_id) const {
    return move_id < firstMoveID;
  }
};

struct NodeTracker {
  vec<std::atomic<SearchID>> searchOfNode;

  SearchID lowestActiveSearchID, highestActiveSearchID;

  explicit NodeTracker(size_t numNodes) : searchOfNode(numNodes), lowestActiveSearchID(1), highestActiveSearchID(0) {}

  bool isSearchInactive(SearchID search_id) const {
    return search_id < lowestActiveSearchID;
  }

  bool isNodeInMySearch(HypernodeID u, SearchID search) const {
    return searchOfNode[u].load(std::memory_order_acq_rel) == search;
  }

  bool isNodeInAnActiveSearch(HypernodeID u) const {
    return searchOfNode[u].load(std::memory_order_acq_rel) >= lowestActiveSearchID;
  }

  bool claimNode(HypernodeID u, SearchID search) {
    SearchID old_search = searchOfNode[u].load(std::memory_order_acq_rel);
    if (isSearchInactive(old_search)) {
      return searchOfNode[u].compare_exchange_strong(old_search, search, std::memory_order_acq_rel);  // one shot
    }
    return false;
  }

  void requestNewSearches(SearchID num_searches) {
    if (highestActiveSearchID >= std::numeric_limits<SearchID>::max() - num_searches - 20) {
      for (auto &x : searchOfNode) x.store(0, std::memory_order_relaxed);
      highestActiveSearchID = 0;
    }
    lowestActiveSearchID = highestActiveSearchID + 1;
    highestActiveSearchID += num_searches;
  }
};


struct FMSharedData {
  // ! For each hyperedge and block, the number of pins_in_part at the beginning of a move phase minus the number of moved out pins
  vec<std::atomic<HypernodeID>> remaining_original_pins;

  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block / the last move to remove a pin from that block
  vec<std::atomic<MoveID>> first_move_in, last_move_out;

  PartitionID numParts;

  GlobalMoveTracker moveTracker;

  NodeTracker nodeTracker;

  FMSharedData(size_t numNodes, size_t numHyperedges, PartitionID numParts) :
          remaining_original_pins(numHyperedges * numParts),
          first_move_in(numHyperedges * numParts),
          last_move_out(numHyperedges * numParts),
          numParts(numParts),
          moveTracker(numNodes),
          nodeTracker(numNodes) {

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
    for (auto &x : last_move_out) x.store(0, std::memory_order_relaxed);    // should be called very rarely
    for (auto &x : first_move_in) x.store(0, std::memory_order_relaxed);
  }

  void setRemainingOriginalPins(PartitionedHypergraph &phg) {
    assert(remaining_original_pins.size() == phg.getPinCountInPartVector().size());
    // let's try if this works
    std::memcpy(remaining_original_pins.data(), phg.getPinCountInPartVector().data(), phg.getPinCountInPartVector().size());
  }
};

}
}