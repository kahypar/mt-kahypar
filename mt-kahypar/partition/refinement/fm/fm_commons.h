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

#include <tbb/parallel_for.h>

namespace mt_kahypar {
namespace refinement {

using BlockPriorityQueue = ds::ExclusiveHandleHeap< ds::MaxHeap<Gain, PartitionID> >;
using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;    // these need external handles

struct GlobalMoveTracker {
  vec<Move> moveOrder;
  CAtomic<MoveID> runningMoveID;
  MoveID firstMoveID = 1;

  explicit GlobalMoveTracker(size_t numNodes) : moveOrder(numNodes), runningMoveID(1) {}

  // Returns true if stored move IDs should be reset
  bool reset() {
    if (runningMoveID.load() >= std::numeric_limits<MoveID>::max() - moveOrder.size() - 20) {
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
    assert(move_id - firstMoveID < moveOrder.size());
    moveOrder[move_id - firstMoveID] = m;
    return move_id;
  }

  Move& getMove(MoveID move_id) {
    assert(move_id - firstMoveID < moveOrder.size());
    return moveOrder[move_id - firstMoveID];
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

  SearchID deactivatedNodeMarker = 1;
  CAtomic<SearchID> highestActiveSearchID { 1 };

  explicit NodeTracker(size_t numNodes) : searchOfNode(numNodes) {
    for (auto& x : searchOfNode) {
      x.store(0, std::memory_order_relaxed);
    }
  }

  // only the search that owns u is allowed to call this
  void deactivateNode(HypernodeID u, SearchID search_id) {
    assert(searchOfNode[u].load() == search_id);
    unused(search_id);
    searchOfNode[u].store(deactivatedNodeMarker, std::memory_order_acq_rel);
  }

  // should not be called when searches try to claim nodes
  void activateNode(HypernodeID u) {
    searchOfNode[u].store(0, std::memory_order_relaxed);
  }

  bool isSearchInactive(SearchID search_id) const {
    return search_id < deactivatedNodeMarker;
  }

  bool canNodeStartNewSearch(HypernodeID u) const {
    return isSearchInactive( searchOfNode[u].load(std::memory_order_acq_rel) );
  }

  void requestNewSearches(SearchID max_num_searches) {
    if (highestActiveSearchID.load(std::memory_order_relaxed) >= std::numeric_limits<SearchID>::max() - max_num_searches - 20) {
      for (auto& x : searchOfNode) {
        x.store(0, std::memory_order_relaxed);
      }
      highestActiveSearchID.store(0, std::memory_order_relaxed);
    }
    deactivatedNodeMarker = ++highestActiveSearchID;
  }
};


struct FMSharedData {

  //PartitionWeightBudgets partition_weight_budgets;


  vec<PosT> vertexPQHandles;

  PartitionID numParts;

  GlobalMoveTracker moveTracker;

  NodeTracker nodeTracker;

  FMSharedData(size_t numNodes = 0, PartitionID numParts = 0, size_t maxNumThreads = 0) :
          //partition_weight_budgets(static_cast<size_t>(numParts), maxNumThreads),
          vertexPQHandles(numNodes, invalid_position),
          numParts(numParts),
          moveTracker(numNodes),
          nodeTracker(numNodes)
  {
    unused(maxNumThreads);
  }

  /*
  ~FMSharedData() {
    tbb::parallel_invoke(
            [&]() { parallel::parallel_free(remaining_original_pins, first_move_in, last_move_out); },
            [&]() { parallel::parallel_free(vertexPQHandles, moveTracker.moveOrder); }
    );

  }
  */

};

}
}