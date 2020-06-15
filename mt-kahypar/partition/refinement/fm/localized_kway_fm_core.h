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
#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"

namespace mt_kahypar {

class LocalizedKWayFM {

 struct FMLocalData {

   void clear() {
     seedVertices.clear();
     localMoves.clear();
     localMoveIDs.clear();
     runStats.clear();
   }

   void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    utils::MemoryTreeNode* local_data_node = parent->addChild("Local FM Data");
    utils::MemoryTreeNode* seed_vertices_node = local_data_node->addChild("Seed Vertices");
    seed_vertices_node->updateSize(seedVertices.capacity() * sizeof(HypernodeID));
    utils::MemoryTreeNode* local_moves_node = local_data_node->addChild("Local Moves");
    local_moves_node->updateSize(localMoves.capacity() * sizeof(Move));
    utils::MemoryTreeNode* local_move_ids_node = local_data_node->addChild("Local Move IDs");
    local_move_ids_node->updateSize(localMoveIDs.capacity() * sizeof(MoveID));
   }

   // ! Contains all seed vertices of the current local search
   vec<HypernodeID> seedVertices;
   // ! Contains all moves performed during the current local search
   vec<Move> localMoves;
   // ! Contains all move IDs of all commited moves of the current local search
   vec<MoveID> localMoveIDs;
   // ! Stats of the current local search
   FMStats runStats;
 };

 public:
  explicit LocalizedKWayFM(const Context& context, HypernodeID numNodes, PosT* pq_handles) :
          context(context),
          thisSearch(0),
          k(context.partition.k),
          localData(),
          deltaPhg(context.partition.k),
          blockPQ(static_cast<size_t>(k)),
          vertexPQs(static_cast<size_t>(k), VertexPriorityQueue(pq_handles, numNodes)),
          updateDeduplicator(),
          validHyperedges() { }


  bool findMovesUsingFullBoundary(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    localData.clear();
    validHyperedges.clear();
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    for (HypernodeID u : sharedData.refinementNodes.safely_inserted_range()) {
      insertPQ(phg, u, sharedData);
    }
    for (PartitionID i = 0; i < k; ++i) {
      updateBlock(i);
    }

    // this is boundary FM, so it's sequential. no need for delta hypergraph
    internalFindMovesOnGlobalHypergraph(phg, sharedData);
    return true;
  }

  bool findMovesLocalized(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t taskID) {
    localData.clear();
    validHyperedges.clear();

    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    const size_t nSeeds = context.refinement.fm.num_seed_nodes;
    HypernodeID seedNode;
    while (localData.runStats.pushes < nSeeds && sharedData.refinementNodes.try_pop(seedNode, taskID)) {
      if (insertPQ(phg, seedNode, sharedData)) {
        localData.seedVertices.push_back(seedNode);
      }
    }
    for (PartitionID i = 0; i < k; ++i) {
      updateBlock(i);
    }

    if (localData.runStats.pushes > 0) {
      if (!context.refinement.fm.perform_moves_global
          && deltaPhg.combinedMemoryConsumption() > sharedData.deltaMemoryLimitPerThread) {
        sharedData.deltaExceededMemoryConstraints = true;
      }

      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
      }

      if (context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints) {
        internalFindMovesOnGlobalHypergraph(phg, sharedData);
      } else {
        deltaPhg.clear();
        deltaPhg.setPartitionedHypergraph(&phg);
        internalFindMovesOnDeltaHypergraph(phg, sharedData);
      }
      return true;
    } else {
      return false;
    }

  }

private:

  // ! Performs localized FM local search on the delta partitioned hypergraph.
  // ! Moves made by this search are not immediately visible to other concurrent local searches.
  // ! The best prefix of moves is applied to the global partitioned hypergraph after the search finishes.
  void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg,
                                          FMSharedData& sharedData) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    auto hes_to_update_func = [&](const HyperedgeID he,
                                  const HyperedgeWeight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situation.
      // In such cases, we mark the hyperedge as invalid and update the gain of all
      // pins afterwards.
      if ( pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
           pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2 ) {
        validHyperedges[he] = false;
      }
    };

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
           && findNextMove(deltaPhg, move)) {
      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);

      bool moved = false;
      if (move.to != kInvalidPartition) {
        heaviestPartWeight = metrics::heaviestPartAndWeight(deltaPhg).second;
        fromWeight = deltaPhg.partWeight(move.from);
        toWeight = deltaPhg.partWeight(move.to);
        moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
          context.partition.max_part_weights[move.to], hes_to_update_func);
      }

      if (moved) {
        localData.runStats.moves++;
        estimatedImprovement += move.gain;
        localData.localMoves.push_back(move);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localData.localMoves.size();
        }

        insertOrUpdateNeighbors(deltaPhg, sharedData, move);
      }

      for (PartitionID i = 0; i < k; ++i) {
        updateBlock(i);
      }
    }

    std::tie(bestImprovement, bestImprovementIndex) =
      applyMovesOnGlobalHypergraph(phg, sharedData, bestImprovementIndex, bestImprovement);
    localData.runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    localData.runStats.merge(stats);
  }

  // ! Performs FM local search on the global partitioned hypergraph (localized or full boundary).
  // ! Moves made by this search are immediately visible to other concurrent local searches.
  // ! After the search finishes, the moves are rolled back to the best prefix of moves found.
  void internalFindMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                           FMSharedData& sharedData) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;


    auto hes_to_update_func = [&](const HyperedgeID he,
                                  const HyperedgeWeight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      // In such cases, we mark the hyperedge as invalid and update the gain of all
      // pins afterwards.
      if ( pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
           pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2 ) {
        validHyperedges[he] = false;
      }
    };

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
           && findNextMove(phg, move)) {
      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      auto report_success = [&] { move_id = sharedData.moveTracker.insertMove(move); };

      bool moved = false;
      if (move.to != kInvalidPartition) {
        heaviestPartWeight = metrics::heaviestPartAndWeight(phg).second;
        fromWeight = phg.partWeight(move.from);
        toWeight = phg.partWeight(move.to);
        moved = phg.changeNodePartFullUpdate(move.node, move.from, move.to,
          context.partition.max_part_weights[move.to],
          report_success, hes_to_update_func);
      }

      if (moved) {
        localData.runStats.moves++;
        estimatedImprovement += move.gain;
        localData.localMoveIDs.push_back(move_id);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localData.localMoveIDs.size();
        }

        insertOrUpdateNeighbors(phg, sharedData, move);
      }

      for (PartitionID i = 0; i < k; ++i) {
        updateBlock(i);
      }
    }

    revertToBestLocalPrefix(phg, sharedData, bestImprovementIndex);
    localData.runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    localData.runStats.merge(stats);
  }

  void clearPQs(FMSharedData& sharedData,
                const size_t bestImprovementIndex) {
    // release all nodes that were not moved
    // reinsert into task queue only if we're doing multitry and at least one node was moved
    // unless a node was moved, only seed nodes are in the pqs
    const bool release = context.refinement.fm.release_nodes
                         && context.refinement.fm.algorithm == FMAlgorithm::fm_multitry
                         && localData.runStats.moves > 0;
    const bool reinsert_seeds = bestImprovementIndex > 0;

    if (release) {
      if (!reinsert_seeds) {
        for (HypernodeID u : localData.seedVertices) {
          sharedData.fruitlessSeed.set(u, true);
        }
      }

      for (PartitionID i = 0; i < k; ++i) {
        for (PosT j = 0; j < vertexPQs[i].size(); ++j) {
          const HypernodeID node = vertexPQs[i].at(j);
          sharedData.nodeTracker.releaseNode(node);
          if (!sharedData.fruitlessSeed[node] && sharedData.refinementNodes.was_pushed_and_removed(node)) {
            sharedData.refinementNodes.concurrent_push(node);
            localData.runStats.task_queue_reinsertions++;
          }
        }
      }
    }

    for (PartitionID i = 0; i < k; ++i) {
      vertexPQs[i].clear();
    }
    blockPQ.clear();
  }

  void updateBlock(PartitionID i) {
    if (!vertexPQs[i].empty()) {
      blockPQ.insertOrAdjustKey(i, vertexPQs[i].topKey());
    } else if (blockPQ.contains(i)) {
      blockPQ.remove(i);
    }
  }

  template<typename PHG>
  void insertOrUpdateNeighbors(const PHG& phg,
                               FMSharedData& sharedData,
                               const Move& move) {
    for (HyperedgeID e : phg.incidentEdges(move.node)) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold && !validHyperedges[e]) {
        for (HypernodeID v : phg.pins(e)) {
          if (!updateDeduplicator.contains(v)) {
            updateDeduplicator[v] = { };  // insert
            insertOrUpdatePQ(phg, v, sharedData, move);
          }
        }
        validHyperedges[e] = true;
      }
    }
    updateDeduplicator.clear();
  }

  bool insertPQ(const PartitionedHypergraph& phg, const HypernodeID v, FMSharedData& sharedData) {
    NodeTracker& nt = sharedData.nodeTracker;
    SearchID searchOfV = nt.searchOfNode[v].load(std::memory_order_acq_rel);
    if (nt.isSearchInactive(searchOfV)) {
      if (nt.searchOfNode[v].compare_exchange_strong(searchOfV, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID pv = phg.partID(v);
        auto [target, gain] = bestDestinationBlock(phg, v);
        sharedData.targetPart[v] = target;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
        localData.runStats.pushes++;
        return true;
      }
    }
    return false;
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool insertOrUpdatePQ(const PHG& phg,
                        const HypernodeID v,
                        FMSharedData& sharedData,
                        const Move& move) {
    assert(move.from != kInvalidPartition && move.to != kInvalidPartition);
    NodeTracker& nt = sharedData.nodeTracker;
    SearchID searchOfV = nt.searchOfNode[v].load(std::memory_order_acq_rel);
    // Note. Deactivated nodes have a special active search ID so that neither branch is executed
    if (nt.isSearchInactive(searchOfV)) {
      if (nt.searchOfNode[v].compare_exchange_strong(searchOfV, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID pv = phg.partID(v);
        auto [target, gain] = bestDestinationBlock(phg, v);
        sharedData.targetPart[v] = target;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
        localData.runStats.pushes++;
        return true;
      }
    } else if (searchOfV == thisSearch) {
      const PartitionID pv = phg.partID(v);
      assert(vertexPQs[pv].contains(v));
      const PartitionID designatedTargetV = sharedData.targetPart[v];
      Gain gain = 0;
      PartitionID newTarget = kInvalidPartition;

      if (k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
        // moveToPenalty of designatedTargetV is affected.
        // and may now be greater than that of other blocks --> recompute full
        std::tie(newTarget, gain) = bestDestinationBlock(phg, v);
      } else {
        // moveToPenalty of designatedTargetV is not affected.
        // only move.from and move.to may be better
        std::tie(newTarget, gain) = bestOfThree(phg, v, pv, { designatedTargetV, move.from, move.to });
      }
      sharedData.targetPart[v] = newTarget;
      vertexPQs[pv].adjustKey(v, gain);

      return true;
    }
    return false;
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> bestDestinationBlock(const PHG& phg,
                                                               const HypernodeID u) {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const PartitionID from = phg.partID(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i = 0; i < k; ++i) {
      if (i != from) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if ( ( penalty < to_penalty || ( penalty == to_penalty && to_weight < best_to_weight ) ) &&
              to_weight + wu <= context.partition.max_part_weights[i] ) {
          to_penalty = penalty;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    const Gain gain = to != kInvalidPartition ? phg.moveFromBenefit(u) - to_penalty
                      : std::numeric_limits<HyperedgeWeight>::min();
    return std::make_pair(to, gain);
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> bestOfThree(const PHG& phg, HypernodeID u, PartitionID from,
                                                      std::array<PartitionID, 3> parts) {

    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i : parts) {
      if (i != from && i != kInvalidPartition) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if ( ( penalty < to_penalty || ( penalty == to_penalty && to_weight < best_to_weight ) ) &&
             to_weight + wu <= context.partition.max_part_weights[i] ) {
          to_penalty = penalty;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    const Gain gain = to != kInvalidPartition ? phg.moveFromBenefit(u) - to_penalty
                                              : std::numeric_limits<HyperedgeWeight>::min();
    return std::make_pair(to, gain);
  }

  template<typename PHG>
  bool findNextMove(const PHG& phg, Move& m) {
    if (blockPQ.empty()) {
      return false;
    }
    while (true) {
      const PartitionID from = blockPQ.top();
      const HypernodeID u = vertexPQs[from].top();
      const Gain estimated_gain = vertexPQs[from].topKey();
      assert(estimated_gain == blockPQ.topKey());
      auto [to, gain] = bestDestinationBlock(phg, u);
      if (gain >= estimated_gain) { // accept any gain that is at least as good
        m.node = u; m.to = to; m.from = from;
        m.gain = gain;
        localData.runStats.extractions++;
        vertexPQs[from].deleteTop();  // blockPQ updates are done later, collectively.
        return true;
      } else {
        localData.runStats.retries++;
        vertexPQs[from].adjustKey(u, gain);
        if (vertexPQs[from].topKey() != blockPQ.keyOf(from)) {
          blockPQ.adjustKey(from, vertexPQs[from].topKey());
        }
      }
    }
  }

  // ! Makes moves applied on delta hypergraph visible on the global partitioned hypergraph.
  std::pair<Gain, size_t> applyMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                                       FMSharedData& sharedData,
                                                       const size_t bestGainIndex,
                                                       const Gain bestEstimatedImprovement) {
    ASSERT(localData.localMoveIDs.empty());
    Gain estimatedImprovement = 0;
    Gain lastGain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      lastGain += km1Delta(he, edge_weight, edge_size,
        pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain bestImprovement = 0;
    size_t bestIndex = 0;
    for ( size_t i = 0; i < bestGainIndex; ++i ) {
      Move& m = localData.localMoves[i];
      MoveID m_id = std::numeric_limits<MoveID>::max();
      lastGain = 0;
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(),
        [&] { m_id = sharedData.moveTracker.insertMove(m); }, delta_gain_func);
      lastGain = -lastGain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      estimatedImprovement += lastGain;
      ASSERT(m_id != std::numeric_limits<MoveID>::max());
      Move& move = sharedData.moveTracker.getMove(m_id);
      move.gain = lastGain; // Update gain value based on hypergraph delta
      localData.localMoveIDs.push_back(m_id);
      if ( estimatedImprovement >= bestImprovement ) {  // TODO also incorporate balance into this
        bestImprovement = estimatedImprovement;
        bestIndex = i;
      }
    }

    // Kind of double rollback, if gain values are not correct
    ASSERT(localData.localMoveIDs.size() == bestGainIndex);
    if ( estimatedImprovement < 0 ) {
      localData.runStats.local_reverts += localData.localMoves.size() - bestIndex;
      for ( size_t i = bestIndex + 1; i < bestGainIndex; ++i ) {
        Move& m = sharedData.moveTracker.getMove(localData.localMoveIDs[i]);
        phg.changeNodePartFullUpdate(m.node, m.to, m.from);
        sharedData.moveTracker.invalidateMove(m);
      }
      localData.runStats.local_reverts += localData.localMoves.size() - bestGainIndex;
      return std::make_pair(bestImprovement, bestIndex);
    } else {
      return std::make_pair(bestEstimatedImprovement, bestGainIndex);
    }
  }

  // ! Rollback to the best improvement found during local search in case we applied moves
  // ! directly on the global partitioned hypergraph.
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex) {
    localData.runStats.local_reverts += localData.localMoveIDs.size() - bestGainIndex;
    while (localData.localMoveIDs.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localData.localMoveIDs.back());
      phg.changeNodePartFullUpdate(m.node, m.to, m.from);
      sharedData.moveTracker.invalidateMove(m);
      localData.localMoveIDs.pop_back();
    }
  }

 public:
  FMStats stats;

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode* deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(updateDeduplicator.size_in_bytes());
    utils::MemoryTreeNode* valid_hyperedges_node = localized_fm_node->addChild("Valid Hyperedges");
    valid_hyperedges_node->updateSize(validHyperedges.size_in_bytes());
    utils::MemoryTreeNode* block_pq_node = localized_fm_node->addChild("Block PQ");
    block_pq_node->updateSize(blockPQ.size_in_bytes());
    utils::MemoryTreeNode* vertex_pq_node = localized_fm_node->addChild("Vertex PQ");
    for ( const VertexPriorityQueue& pq : vertexPQs ) {
      vertex_pq_node->updateSize(pq.size_in_bytes());
    }
    localData.memoryConsumption(localized_fm_node);
    deltaPhg.memoryConsumption(localized_fm_node);
  }

 private:

  const Context& context;

  // ! Unique search id associated with the current local search
  SearchID thisSearch;

  // ! Number of blocks
  PartitionID k;

  // ! Local data members required for one localized search run
  FMLocalData localData;

  // ! Wrapper around the global partitioned hypergraph, that allows
  // ! to perform moves non-visible for other local searches
  DeltaPartitionedHypergraph deltaPhg;

  // ! Priority Queue that contains for each block of the partition
  // ! the vertex with the best gain value
  BlockPriorityQueue blockPQ;

  // ! From PQs -> For each block it contains the vertices (contained
  // ! in that block) touched by the current local search associated
  // ! with their gain values
  vec<VertexPriorityQueue> vertexPQs;

  // ! After a move it collects all neighbors of the moved vertex
  ds::DynamicSparseSet<HypernodeID> updateDeduplicator;

  // ! Marks all hyperedges that are visited during the local search
  // ! and where the gain of its pin is expected to be equal to gain value
  // ! inside the PQs. A hyperedge can become invalid if a move changes the
  // ! gain values of its pins.
  ds::DynamicSparseMap<HyperedgeID, bool> validHyperedges;
};

}