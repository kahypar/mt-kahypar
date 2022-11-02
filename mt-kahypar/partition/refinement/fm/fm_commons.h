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

#include <mt-kahypar/definitions.h>
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/parallel/work_stack.h>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"

#include <tbb/parallel_for.h>

namespace mt_kahypar {


struct GlobalMoveTracker {
  vec<Move> moveOrder;
  vec<MoveID> moveOfNode;
  CAtomic<MoveID> runningMoveID;
  MoveID firstMoveID = 1;

  explicit GlobalMoveTracker(size_t numNodes = 0) :
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
    moveOrder[move_id - firstMoveID].gain = 0;      // set to zero so the recalculation can safely distribute
    moveOfNode[m.node] = move_id;
    return move_id;
  }

  Move& getMove(MoveID move_id) {
    assert(move_id - firstMoveID < moveOrder.size());
    return moveOrder[move_id - firstMoveID];
  }

  bool wasNodeMovedInThisRound(HypernodeID u) const {
    const MoveID m_id = moveOfNode[u];
    return m_id >= firstMoveID
           && m_id < runningMoveID.load(std::memory_order_relaxed)  // active move ID
           && moveOrder[m_id - firstMoveID].isValid();      // not reverted already
  }

  MoveID numPerformedMoves() const {
    return runningMoveID.load(std::memory_order_relaxed) - firstMoveID;
  }

  bool isMoveStale(const MoveID move_id) const {
    return move_id < firstMoveID;
  }
};

struct NodeTracker {
  vec<CAtomic<SearchID>> searchOfNode;

  SearchID releasedMarker = 1;
  SearchID deactivatedNodeMarker = 2;
  CAtomic<SearchID> highestActiveSearchID { 2 };

  explicit NodeTracker(size_t numNodes = 0) : searchOfNode(numNodes, CAtomic<SearchID>(0)) { }

  // only the search that owns u is allowed to call this
  void deactivateNode(HypernodeID u, SearchID search_id) {
    assert(searchOfNode[u].load() == search_id);
    unused(search_id);
    searchOfNode[u].store(deactivatedNodeMarker, std::memory_order_release);
  }

  bool isLocked(HypernodeID u) {
    return searchOfNode[u].load(std::memory_order_relaxed) == deactivatedNodeMarker;
  }

  void releaseNode(HypernodeID u) {
    searchOfNode[u].store(releasedMarker, std::memory_order_relaxed);
  }

  bool isSearchInactive(SearchID search_id) const {
    return search_id < deactivatedNodeMarker;
  }

  bool canNodeStartNewSearch(HypernodeID u) const {
    return isSearchInactive( searchOfNode[u].load(std::memory_order_relaxed) );
  }

  bool tryAcquireNode(HypernodeID u, SearchID new_search) {
    SearchID current_search = searchOfNode[u].load(std::memory_order_relaxed);
    return isSearchInactive(current_search)
            && searchOfNode[u].compare_exchange_strong(current_search, new_search, std::memory_order_acq_rel);
  }

  void requestNewSearches(SearchID max_num_searches) {
    if (highestActiveSearchID.load(std::memory_order_relaxed) >= std::numeric_limits<SearchID>::max() - max_num_searches - 20) {
      tbb::parallel_for(0UL, searchOfNode.size(), [&](const size_t i) {
        searchOfNode[i].store(0, std::memory_order_relaxed);
      });
      highestActiveSearchID.store(1, std::memory_order_relaxed);
    }
    deactivatedNodeMarker = ++highestActiveSearchID;
    releasedMarker = deactivatedNodeMarker - 1;
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

  // ! Stores the designated target part of a vertex, i.e. the part with the highest gain to which moving is feasible
  vec<PartitionID> targetPart;

  // ! Stop parallel refinement if finishedTasks > finishedTasksLimit to avoid long-running single searches
  CAtomic<size_t> finishedTasks;
  size_t finishedTasksLimit = std::numeric_limits<size_t>::max();

  // ! Switch to applying moves directly if the use of local delta partitions exceeded a memory limit
  bool deltaExceededMemoryConstraints = false;
  size_t deltaMemoryLimitPerThread = 0;

  bool release_nodes = true;
  bool perform_moves_global = true;

  FMSharedData(size_t numNodes = 0, PartitionID numParts = 0, size_t numThreads = 0, size_t numPQHandles = 0) :
          refinementNodes(), //numNodes, numThreads),
          vertexPQHandles(), //numPQHandles, invalid_position),
          numParts(numParts),
          moveTracker(), //numNodes),
          nodeTracker(), //numNodes),
          targetPart()
  {
    finishedTasks.store(0, std::memory_order_relaxed);

    // 128 * 3/2 GB --> roughly 1.5 GB per thread on our biggest machine
    deltaMemoryLimitPerThread = 128UL * (1UL << 30) * 3 / ( 2 * std::max(1UL, numThreads) );

    tbb::parallel_invoke([&] {
      moveTracker.moveOrder.resize(numNodes);
    }, [&] {
      moveTracker.moveOfNode.resize(numNodes);
    }, [&] {
      nodeTracker.searchOfNode.resize(numNodes, CAtomic<SearchID>(0));
    }, [&] {
      vertexPQHandles.resize(numPQHandles, invalid_position);
    }, [&] {
      refinementNodes.tls_queues.resize(numThreads);
    }, [&] {
      targetPart.resize(numNodes, kInvalidPartition);
    });
  }

  FMSharedData(size_t numNodes, const Context& context) :
        FMSharedData(
                numNodes,
                context.partition.k,
                TBBInitializer::instance().total_number_of_threads(),
                getNumberOfPQHandles(context, numNodes)
                )  { }


  size_t getNumberOfPQHandles(const Context& context, size_t numNodes) {
    if (context.refinement.fm.algorithm == FMAlgorithm::fm_gain_delta) {
      return numNodes * context.partition.k;
    } else {
      return numNodes;
    }
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* shared_fm_data_node = parent->addChild("Shared FM Data");

    utils::MemoryTreeNode* pq_handles_node = shared_fm_data_node->addChild("PQ Handles");
    pq_handles_node->updateSize(vertexPQHandles.capacity() * sizeof(PosT));
    utils::MemoryTreeNode* move_tracker_node = shared_fm_data_node->addChild("Move Tracker");
    move_tracker_node->updateSize(moveTracker.moveOrder.capacity() * sizeof(Move) +
                                  moveTracker.moveOfNode.capacity() * sizeof(MoveID));
    utils::MemoryTreeNode* node_tracker_node = shared_fm_data_node->addChild("Node Tracker");
    node_tracker_node->updateSize(nodeTracker.searchOfNode.capacity() * sizeof(SearchID));
    refinementNodes.memoryConsumption(shared_fm_data_node);
  }
};

struct FMStats {
  size_t retries = 0;
  size_t extractions = 0;
  size_t pushes = 0;
  size_t moves = 0;
  size_t local_reverts = 0;
  size_t task_queue_reinsertions = 0;
  size_t best_prefix_mismatch = 0;
  Gain estimated_improvement = 0;


  void clear() {
    retries = 0;
    extractions = 0;
    pushes = 0;
    moves = 0;
    local_reverts = 0;
    task_queue_reinsertions = 0;
    best_prefix_mismatch = 0;
    estimated_improvement = 0;
  }

  void merge(FMStats& other) {
    other.retries += retries;
    other.extractions += extractions;
    other.pushes += pushes;
    other.moves += moves;
    other.local_reverts += local_reverts;
    other.task_queue_reinsertions += task_queue_reinsertions;
    other.best_prefix_mismatch += best_prefix_mismatch;
    other.estimated_improvement += estimated_improvement;
    clear();
  }

  std::string serialize() const {
    std::stringstream os;
    os  << V(retries) << " " << V(extractions) << " " << V(pushes) << " "
        << V(moves) << " " << V(local_reverts) << " " << V(best_prefix_mismatch);
    return os.str();
  }
};


}