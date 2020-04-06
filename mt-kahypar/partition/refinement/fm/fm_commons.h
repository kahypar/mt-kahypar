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
#include <mt-kahypar/datastructures/priority_queue.h>
#include "partition_weight_budgets.h"

namespace mt_kahypar {
namespace refinement {

using BlockPriorityQueue = ds::ExclusiveHandleHeap< ds::MaxHeap<Gain, PartitionID> >;

struct MoveTo { // TODO should we store the target block in the PQ as well, or should we perform another call to bestDestinationBlock upon extraction
  HypernodeID node = invalidNode;
  PartitionID toBlock = kInvalidPartition;
};

using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;    // these need external handles


struct GlobalMoveTracker {
  vec<Move> globalMoveOrder;
  CAtomic<MoveID> runningMoveID;
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

  Move& getMove(MoveID move_id) {
    return globalMoveOrder[move_id - firstMoveID];
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

  static constexpr SearchID deactivatedNodeMarker = std::numeric_limits<SearchID>::max();
  SearchID lowestActiveSearchID = 1;
  std::atomic<SearchID> highestActiveSearchID { 0 };

  explicit NodeTracker(size_t numNodes) : searchOfNode(numNodes) {
    for (auto& x : searchOfNode) {
      x.store(0, std::memory_order_relaxed);
    }
  }

  // only the search that owns u is allowed to call this
  void deactivateNode(HypernodeID u, SearchID search_id) {
    assert(searchOfNode[u].load() == search_id);
    searchOfNode[u].store(deactivatedNodeMarker, std::memory_order_acq_rel);
  }

  // should not be called when searches try to claim nodes
  void activateNode(HypernodeID u) {
    searchOfNode[u].store(0, std::memory_order_relaxed);
  }

  bool isSearchInactive(SearchID search_id) const {
    return search_id < lowestActiveSearchID;
  }

  bool canNodeStartNewSearch(HypernodeID u) const {
    return isSearchInactive( searchOfNode[u].load(std::memory_order_acq_rel) );
  }

  void requestNewSearches(SearchID max_num_searches) {
    if (highestActiveSearchID.load(std::memory_order_relaxed) >= std::numeric_limits<SearchID>::max() - max_num_searches - 20) {
      for (auto& x : searchOfNode) {
        x.store(0, std::memory_order_relaxed);
      }
      highestActiveSearchID = 0;
    }
    lowestActiveSearchID = highestActiveSearchID + 1;
  }
};


struct FMSharedData {

  //PartitionWeightBudgets partition_weight_budgets;

  // ! For each hyperedge and block, the number of pins_in_part at the beginning of a move phase minus the number of moved out pins
  vec<CAtomic<HypernodeID>> remaining_original_pins;

  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block / the last move to remove a pin from that block
  vec<CAtomic<MoveID>> first_move_in, last_move_out;

  vec<PosT> vertexPQHandles;

  PartitionID numParts;

  GlobalMoveTracker moveTracker;

  NodeTracker nodeTracker;

  FMSharedData(size_t numNodes, size_t numHyperedges, PartitionID numParts, size_t maxNumThreads) :
          //partition_weight_budgets(static_cast<size_t>(numParts), maxNumThreads),
          remaining_original_pins(numHyperedges * numParts),
          first_move_in(numHyperedges * numParts),
          last_move_out(numHyperedges * numParts),
          vertexPQHandles(numNodes, invalid_position),
          numParts(numParts),
          moveTracker(numNodes),
          nodeTracker(numNodes)
  {
    unused(maxNumThreads);
    resetStoredMoveIDs();
  }

  ~FMSharedData() {
    tbb::parallel_invoke(
            [&]() { parallel::parallel_free(remaining_original_pins, first_move_in, last_move_out); },
            [&]() { parallel::parallel_free(vertexPQHandles, moveTracker.globalMoveOrder); }
    );

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
    assert(remaining_original_pins.size() == phg.getPinCountInPartVector().size());
    size_t n = phg.getPinCountInPartVector().size();
    std::copy_n(phg.getPinCountInPartVector().begin(), n, remaining_original_pins.begin());
  }

  void performHyperedgeSpecificMoveUpdates(MoveID move_id, HyperedgeID e) {
    // TODO more pruning!
    Move& m = moveTracker.getMove(move_id);
    // update first move in
    CAtomic<MoveID>& fmi = first_move_in[e * numParts + m.to];
    MoveID expected = fmi.load(std::memory_order_acq_rel);
    while ((moveTracker.isIDStale(expected) || expected > move_id)
           && !fmi.compare_exchange_weak(expected, move_id, std::memory_order_acq_rel)) { }

    // update last move out
    CAtomic<MoveID>& lmo = last_move_out[e * numParts + m.from];
    expected = lmo.load(std::memory_order_acq_rel);
    while (expected < move_id && !lmo.compare_exchange_weak(expected, move_id, std::memory_order_acq_rel)) { }

    remaining_original_pins[e * numParts + m.from].fetch_sub(1, std::memory_order_relaxed);
  }
};

}
}