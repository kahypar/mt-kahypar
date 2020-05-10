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
#include <mt-kahypar/parallel/work_stack.h>

#include <tbb/parallel_for.h>

namespace mt_kahypar {

using BlockPriorityQueue = ds::ExclusiveHandleHeap< ds::MaxHeap<Gain, PartitionID> >;
using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;    // these need external handles

struct GlobalMoveTracker {
  vec<Move> moveOrder;
  vec<MoveID> moveOfNode;
  CAtomic<MoveID> runningMoveID;
  MoveID firstMoveID = 1;

  explicit GlobalMoveTracker(size_t numNodes) :
          moveOrder(numNodes),
          moveOfNode(numNodes, 0),
          runningMoveID(1) { }

  // Returns true if stored move IDs should be reset
  bool reset() {
    if (runningMoveID.load() >= std::numeric_limits<MoveID>::max() - moveOrder.size() - 20) {
      tbb::parallel_for(0UL, moveOfNode.size(), [&](size_t i) { moveOfNode[i] = 0; }, tbb::static_partitioner());
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

  bool isMoveStillValid(MoveID move_id) {
    return getMove(move_id).gain != invalidGain;
  }

  bool isMoveStillValid(const Move& m) const {
    return m.gain != invalidGain;
  }

  void invalidateMove(MoveID move_id) {
    getMove(move_id).gain = invalidGain;
  }

  void invalidateMove(Move& m) {
    m.gain = invalidGain;
  }

  MoveID numPerformedMoves() const {
    return runningMoveID.load(std::memory_order_relaxed) - firstMoveID;
  }

  bool isMoveStale(const MoveID move_id) const {
    return move_id < firstMoveID;
  }
};

struct NodeTracker {
  vec<std::atomic<SearchID>> searchOfNode;

  SearchID deactivatedNodeMarker = 1;
  CAtomic<SearchID> highestActiveSearchID { 1 };

  explicit NodeTracker(size_t numNodes) :
          searchOfNode(numNodes)
  {
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

  bool isLocked(HypernodeID u) {
    return searchOfNode[u].load(std::memory_order_relaxed) == deactivatedNodeMarker;
  }

  // should not be called when searches try to claim nodes
  void releaseNode(HypernodeID u) {
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
      tbb::parallel_for(0UL, searchOfNode.size(), [&](const size_t i) {
        searchOfNode[i].store(0, std::memory_order_relaxed);
      });
      highestActiveSearchID.store(0, std::memory_order_relaxed);
    }
    deactivatedNodeMarker = ++highestActiveSearchID;
  }
};


struct FMSharedData {

  // ! Nodes to initialize the localized FM searches with
  WorkContainer<HypernodeID> refinementNodes;

  // ! PQ handles shared by all threads (each vertex is only held by one thread)
  vec<PosT> vertexPQHandles;

  // ! num parts
  PartitionID numParts;

  // ! Stores the sequence of performed moves and assigns IDs to moves that can be used in the global rollback code
  GlobalMoveTracker moveTracker;

  // ! Tracks the current search of a node, and if a node can still be added to an active search
  NodeTracker nodeTracker;

  FMSharedData(size_t numNodes = 0, PartitionID numParts = 0) :
          refinementNodes(numNodes),
          vertexPQHandles(numNodes, invalid_position),
          numParts(numParts),
          moveTracker(numNodes),
          nodeTracker(numNodes) { }

  FMSharedData(size_t numNodes, const Context& context) : FMSharedData(numNodes, context.partition.k) { }

};

struct FMStats {
  size_t retries = 0;
  size_t extractions = 0;
  size_t pushes = 0;
  size_t moves = 0;
  size_t local_reverts = 0;
  Gain estimated_improvement = 0;


  void clear() {
    retries = 0;
    extractions = 0;
    pushes = 0;
    moves = 0;
    local_reverts = 0;
    estimated_improvement = 0;
  }

  void merge(FMStats& other) {
    other.retries += retries;
    other.extractions += extractions;
    other.pushes += pushes;
    other.moves += moves;
    other.local_reverts += local_reverts;
    other.estimated_improvement += estimated_improvement;
    clear();
  }

  std::string serialize() {
    std::stringstream os;
    os << V(retries) << " " << V(extractions) << " " << V(pushes) << " " << V(moves) << " " << V(local_reverts);
    return os.str();
  }
};

}