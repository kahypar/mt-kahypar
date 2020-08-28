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

#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"
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
   // ! Contains all move IDs of all committed moves of the current local search
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


  bool findMovesUsingFullBoundary(PartitionedHypergraph& phg, FMSharedData& sharedData);

  bool findMovesLocalized(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t taskID);

private:

  // ! Performs localized FM local search on the delta partitioned hypergraph.
  // ! Moves made by this search are not immediately visible to other concurrent local searches.
  // ! The best prefix of moves is applied to the global partitioned hypergraph after the search finishes.
  void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg,
                                          FMSharedData& sharedData);


  void internalFindMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                           FMSharedData& sharedData);

  void clearPQs(FMSharedData& sharedData, size_t bestImprovementIndex);

  void updateBlock(PartitionID i) {
    if (!vertexPQs[i].empty()) {
      blockPQ.insertOrAdjustKey(i, vertexPQs[i].topKey());
    } else if (blockPQ.contains(i)) {
      blockPQ.remove(i);
    }
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void insertOrUpdateNeighbors(const PHG& phg,
                               FMSharedData& sharedData,
                               const Move& move) {
    for (HyperedgeID e : phg.incidentEdges(move.node)) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold && !validHyperedges[e]) { // TODO why not just store them in a vector to later iterate over?
        for (HypernodeID v : phg.pins(e)) {
          if (!updateDeduplicator.contains(v)) {
            insertOrUpdatePQ(phg, v, sharedData, move);
            updateDeduplicator[v] = { };  // insert
          }
        }
        validHyperedges[e] = true;
      }
    }
    updateDeduplicator.clear();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
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
      if ( vertexPQs[pv].contains(v) ) {
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
      } else {
        // these are reinserts of vertices that changed their status from boundary to non-boundary vertex and back again
        auto [target, gain] = bestDestinationBlock(phg, v);
        sharedData.targetPart[v] = target;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
        localData.runStats.boundary_toggles++;
      }

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
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
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
                                                       size_t bestGainIndex,
                                                       Gain bestEstimatedImprovement);

  // ! Rollback to the best improvement found during local search in case we applied moves
  // ! directly on the global partitioned hypergraph.
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex);

 public:
  FMStats stats;

  void memoryConsumption(utils::MemoryTreeNode* parent) const ;

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
  ds::DeltaPartitionedHypergraph<PartitionedHypergraph> deltaPhg;

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