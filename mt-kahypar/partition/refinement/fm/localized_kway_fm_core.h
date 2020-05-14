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

#include "fm_commons.h"
#include "clearlist.hpp"
#include "stop_rule.h"

namespace mt_kahypar {

class LocalizedKWayFM {
public:
  explicit LocalizedKWayFM(const Context& context, HypernodeID numNodes, PosT* pq_handles) :
          numParts(context.partition.k),
          delta_phg(context.partition.k),
          blockPQ(static_cast<size_t>(numParts)),
          vertexPQs(static_cast<size_t>(numParts), VertexPriorityQueue(pq_handles, numNodes)),
          updateDeduplicator(numNodes),
          context(context)
  {
    maxPartWeight = context.partition.max_part_weights[0];
  }

  bool findMoves(PartitionedHypergraph& phg, FMSharedData& sharedData, vec<HypernodeID>& seedNodes) {
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    for (HypernodeID u : seedNodes) {
      insertOrUpdatePQ(phg, u, sharedData.nodeTracker);
    }
    for (PartitionID i = 0; i < numParts; ++i) {
      updateBlock(i);
    }

    // this is boundary FM, so it's sequential. no need for delta hypergraph
    internalFindMovesOnGlobalHypergraph(phg, sharedData);
    return true;
  }



  bool findMoves(PartitionedHypergraph& phg, FMSharedData& sharedData) {

    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
    const size_t nSeeds = context.refinement.fm.num_seed_nodes;
    HypernodeID seedNode;
    while (runStats.pushes < nSeeds && sharedData.refinementNodes.try_pop(seedNode)) {
      if (!updateDeduplicator.contains(seedNode) && insertOrUpdatePQ(phg, seedNode, sharedData.nodeTracker)) {
        seeds.push_back(seedNode);
      }
    }
    updateBlocks(phg, kInvalidPartition);

    if (runStats.pushes > 0) {
      if ( context.refinement.fm.perform_moves_global ) {
        internalFindMovesOnGlobalHypergraph(phg, sharedData);
      } else {
        delta_phg.clear();
        delta_phg.setPartitionedHypergraph(&phg);
        internalFindMovesOnDeltaHypergraph(phg, sharedData);
      }
      return true;
    } else {
      return false;
    }

  }

private:

  // ! Starts a localized FM search on the delta partitioned hypergraph. Moves
  // ! that are made by this local search are not immediatly visible to other
  // ! concurrent running local searches. Moves are applied to global hypergraph,
  // ! if search yield to an improvement.
  void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg,
                                          FMSharedData& sharedData) {
    localMoves.clear();
    StopRule stopRule(phg.initialNumNodes());
    Move m;

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight toWeight = 0;

    while (!stopRule.searchShouldStop() && findNextMove(delta_phg, m)) {
      sharedData.nodeTracker.deactivateNode(m.node, thisSearch);

      bool moved = false;
      if (m.to != kInvalidPartition) {
        heaviestPartWeight = metrics::heaviestPartAndWeight(delta_phg).second;
        const HypernodeWeight fromWeight = delta_phg.partWeight(m.from);
        toWeight = delta_phg.partWeight(m.to);
        moved = delta_phg.changeNodePart(m.node, m.from, m.to, std::max(maxPartWeight, fromWeight));
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += m.gain;
        localMoves.push_back(m);
        stopRule.update(m.gain);

        // Check if move improves current best solution
        bool move_improved_quality = false;
        if ( context.refinement.fm.allow_zero_gain_moves ) {
          move_improved_quality = estimatedImprovement >= bestImprovement;
        } else {
          const bool improved_km1 = estimatedImprovement > bestImprovement;
          const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement &&
                                                       toWeight + phg.nodeWeight(m.node) < heaviestPartWeight;
          move_improved_quality = improved_km1 || improved_balance_less_equal_km1;
        }

        if (move_improved_quality) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();
        }

        insertOrUpdateNeighbors(delta_phg, sharedData, m.node);
      }
      updateBlocks(delta_phg, m.from);
    }

    std::tie(bestImprovement, bestImprovementIndex) =
      applyMovesOnGlobalHypergraph(phg, sharedData, bestImprovementIndex, bestImprovement);
    runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    runStats.merge(stats);
  }

  // ! Starts a localized FM search on the global partitioned hypergraph. Moves
  // ! that are made by this local search are immediatly visible to other concurrent
  // ! running local searches. Moves are rollbacked to best improvement during
  // ! that search.
  void internalFindMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                           FMSharedData& sharedData) {
    localAppliedMoves.clear();
    StopRule stopRule(phg.initialNumNodes());
    Move m;

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight toWeight = 0;

    while (!stopRule.searchShouldStop() && findNextMove(phg, m)) {
      sharedData.nodeTracker.deactivateNode(m.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      auto report_success = [&] { move_id = sharedData.moveTracker.insertMove(m); };

      bool moved = false;
      if (m.to != kInvalidPartition) {
        heaviestPartWeight = metrics::heaviestPartAndWeight(phg).second;
        const HypernodeWeight fromWeight = phg.partWeight(m.from);
        toWeight = phg.partWeight(m.to);
        moved = phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::max(maxPartWeight, fromWeight), report_success);
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += m.gain;
        localAppliedMoves.push_back(move_id);
        stopRule.update(m.gain);

        // Check if move improves current best solution
        bool move_improved_quality = false;
        if ( context.refinement.fm.allow_zero_gain_moves ) {
          move_improved_quality = estimatedImprovement >= bestImprovement;
        } else {
          const bool improved_km1 = estimatedImprovement > bestImprovement;
          const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement &&
                                                      toWeight + phg.nodeWeight(m.node) < heaviestPartWeight;
          move_improved_quality = improved_km1 || improved_balance_less_equal_km1;
        }

        if (move_improved_quality) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localAppliedMoves.size();
        }

        insertOrUpdateNeighbors(phg, sharedData, m.node);
      }
      updateBlocks(phg, m.from);
    }

    revertToBestLocalPrefix(phg, sharedData, bestImprovementIndex);
    runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    runStats.merge(stats);
  }

  void clearPQs(FMSharedData& sharedData, size_t bestImprovementIndex) {
    // release all nodes that were not moved
    // reinsert into task queue only if we're doing multitry and at least one node was moved
    // unless a node was moved, only seed nodes are in the pqs
    const bool release = context.refinement.fm.algorithm == FMAlgorithm::fm_multitry && runStats.moves > 0;
    const bool reinsert_seeds = bestImprovementIndex > 0;

    if (release) {
      if (!reinsert_seeds) {
        for (HypernodeID u : seeds) {
          sharedData.fruitlessSeed.set(u, true);
        }
      }

      for (PartitionID i = 0; i < numParts; ++i) {
        for (PosT j = 0; j < vertexPQs[i].size(); ++j) {
          const HypernodeID node = vertexPQs[i].at(j);
          sharedData.nodeTracker.releaseNode(node);
          if (!sharedData.fruitlessSeed[node] && sharedData.refinementNodes.was_pushed_and_removed(node)) {
            sharedData.refinementNodes.push(node);
          }
        }
      }
    }

    seeds.clear();
    for (PartitionID i = 0; i < numParts; ++i) {
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
  void updateBlocks(const PHG& phg, const PartitionID moved_from) {
    if (moved_from == kInvalidPartition || updateDeduplicator.size() >= size_t(numParts)) {
      for (PartitionID i = 0; i < numParts; ++i) {
        updateBlock(i);
      }
    } else {
      updateBlock(moved_from);
      for (const HypernodeID v : updateDeduplicator.keys()) {
        updateBlock(phg.partID(v));
      }
    }
    updateDeduplicator.clear();
  }

  template<typename PHG>
  void insertOrUpdateNeighbors(const PHG& phg,
                               FMSharedData& sharedData,
                               const HypernodeID u) {
    for (HyperedgeID e : phg.incidentEdges(u)) {
      if (phg.edgeSize(e) < context.partition.hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (!updateDeduplicator.contains(v)) {
            updateDeduplicator.insert(v);
            insertOrUpdatePQ(phg, v, sharedData.nodeTracker);
          }
        }
      }
    }
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool insertOrUpdatePQ(const PHG& phg,
                        const HypernodeID v,
                        NodeTracker& nt) {
    SearchID searchOfV = nt.searchOfNode[v].load(std::memory_order_acq_rel);
    // Note. Deactivated nodes have a special active search ID so that neither branch is executed
    if (nt.isSearchInactive(searchOfV)) {
      if (nt.searchOfNode[v].compare_exchange_strong(searchOfV, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID pv = phg.partID(v);
        const Gain gain = bestDestinationBlock(phg, v).second;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
        runStats.pushes++;
        return true;
      }
    } else if (searchOfV == thisSearch) {
      const PartitionID pv = phg.partID(v);
      assert(vertexPQs[pv].contains(v));
      const Gain gain = bestDestinationBlock(phg, v).second;
      vertexPQs[pv].adjustKey(v, gain);
      // if pv == move.from or pv == move.to only the gains of move.from and move.to could change
      // if these gains are better than vertexPQ[pv].keyOf(v), we could increase the key
      // however, this is incorrect if this entry is for the target move.from or move.to.
      // since we don't store this information (on purpose), we can't easily figure that out
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
    for (PartitionID i = 0; i < numParts; ++i) {
      if (i != from) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if ( ( penalty < to_penalty || ( penalty == to_penalty && to_weight < best_to_weight ) ) &&
              ( to_weight + wu <= maxPartWeight || to_weight < std::max(best_to_weight, maxPartWeight + 1) ) ) {
          to_penalty = penalty;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    return std::make_pair(to, phg.moveFromBenefit(u) - to_penalty);
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
        runStats.extractions++;
        vertexPQs[from].deleteTop();  // blockPQ updates are done later, collectively.
        return true;
      } else {
        runStats.retries++;
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
    localAppliedMoves.clear();
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
      Move& m = localMoves[i];
      MoveID m_id = std::numeric_limits<MoveID>::max();
      lastGain = 0;
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(),
        [&] { m_id = sharedData.moveTracker.insertMove(m); }, delta_gain_func);
      lastGain = -lastGain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      estimatedImprovement += lastGain;
      ASSERT(m_id != std::numeric_limits<MoveID>::max());
      Move& move = sharedData.moveTracker.getMove(m_id);
      move.gain = lastGain; // Update gain value based on hypergraph delta
      localAppliedMoves.push_back(m_id);
      if ( estimatedImprovement >= bestImprovement ) {
        bestImprovement = estimatedImprovement;
        bestIndex = i;
      }
    }

    runStats.local_reverts += localMoves.size() - bestIndex;

    // Kind of double rollback, if gain values are not correct
    ASSERT(localAppliedMoves.size() == bestGainIndex);
    if ( estimatedImprovement < 0 ) {
      for ( size_t i = bestIndex + 1; i < bestGainIndex; ++i ) {
        Move& m = sharedData.moveTracker.getMove(localAppliedMoves[i]);
        phg.changeNodePartFullUpdate(m.node, m.to, m.from);
        sharedData.moveTracker.invalidateMove(m);
      }
      return std::make_pair(bestImprovement, bestIndex);
    } else {
      return std::make_pair(bestEstimatedImprovement, bestGainIndex);
    }
  }

  // ! Rollback to the best improvement found during local search in case we applied moves
  // ! directly on the global partitioned hypergraph.
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex) {
    runStats.local_reverts += localAppliedMoves.size() - bestGainIndex;
    while (localAppliedMoves.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localAppliedMoves.back());
      phg.changeNodePartFullUpdate(m.node, m.to, m.from);
      sharedData.moveTracker.invalidateMove(m);
      localAppliedMoves.pop_back();
    }
  }

  SearchID thisSearch;
  PartitionID numParts;

  DeltaPartitionedHypergraph delta_phg;

  BlockPriorityQueue blockPQ;
  vec<VertexPriorityQueue> vertexPQs;

  // Note: prefer ClearListSet over SparseSet because
  // ClearListSet takes numNodes + numInsertedNodes*32 bit
  // SparseSet takes 2 * numNodes * 32 bit
  // where numInsertedNodes is presumably much smaller than numNodes
  ldc::ClearListSet<HypernodeID> updateDeduplicator;

  const Context& context;
  HypernodeWeight maxPartWeight = 0;
  FMStats runStats;
  std::vector<HypernodeID> seeds;
  vec<Move> localMoves;
  vec<MoveID> localAppliedMoves;
public:
  FMStats stats;
};

}