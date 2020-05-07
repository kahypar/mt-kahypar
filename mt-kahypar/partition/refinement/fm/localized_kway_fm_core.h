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
  // unfortunately the compiler thinks we're trying to pass a const-ref for the pq_handles, which we don't. therefore it had to be a pointer :(
  explicit LocalizedKWayFM(const Context& context, HypernodeID numNodes, PosT* pq_handles) :
          numParts(context.partition.k),
          delta_phg(context.partition.k),
          blockPQ(static_cast<size_t>(numParts)),
          vertexPQs(static_cast<size_t>(numParts), VertexPriorityQueue(pq_handles, numNodes)),
          updateDeduplicator(numNodes),
          context(context)
  {
    maxPartWeight = context.partition.max_part_weights[0];
    perfectBalancePartWeight = context.partition.perfect_balance_part_weights[0];
    minPartWeight = static_cast<HypernodeWeight>(std::floor(perfectBalancePartWeight * (1 - context.partition.epsilon)));
  }

  bool findMoves(PartitionedHypergraph& phg, FMSharedData& sharedData, vec<HypernodeID>& initialNodes) {
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
    delta_phg.clear();
    delta_phg.setPartitionedHypergraph(&phg);

    for (HypernodeID u : initialNodes) {
      insertOrUpdatePQ(phg, u, sharedData.nodeTracker);
    }

    for (PartitionID i = 0; i < numParts; ++i) {
      updateBlock(i);
    }

    internalFindMoves(phg, sharedData);

    return true;
  }

  bool findMoves(PartitionedHypergraph& phg, FMSharedData& sharedData, HypernodeID initialBorderNode = invalidNode) {
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
    delta_phg.clear();
    delta_phg.setPartitionedHypergraph(&phg);

    if (initialBorderNode == invalidNode) {
      const size_t nSeeds = numberOfSeedNodes(phg.initialNumNodes());
      while (runStats.pushes <= nSeeds && sharedData.refinementNodes.try_pop(initialBorderNode)) {
        if (!updateDeduplicator.contains(initialBorderNode)
            && insertOrUpdatePQ(phg, initialBorderNode, sharedData.nodeTracker)
            && context.refinement.fm.init_localized_search_with_neighbors) {
          updateDeduplicator.insert(initialBorderNode);
          insertOrUpdateNeighbors(phg, sharedData, initialBorderNode);
        }
      }
      updateBlocks(kInvalidPartition);
    } else {
      if (insertOrUpdatePQ(phg, initialBorderNode, sharedData.nodeTracker)) {
        if (context.refinement.fm.init_localized_search_with_neighbors) {
          updateDeduplicator.insert(initialBorderNode);
          insertOrUpdateNeighbors(phg, sharedData, initialBorderNode);
          updateBlocks(delta_phg.partID(initialBorderNode));
        } else {
          updateBlock(delta_phg.partID(initialBorderNode));
        }
      }
    }

    if (runStats.pushes > 0) {
      internalFindMoves(phg, sharedData);
      return true;
    } else {
      return false;
    }
  }

private:

  // assumes the PQs have been initialized
  void internalFindMoves(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    localMoves.clear();
    StopRule stopRule(phg.initialNumNodes());
    Move m;
    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0, bestImprovement = 0;
    while (!stopRule.searchShouldStop() && findNextMove(phg, m)) {
      sharedData.nodeTracker.deactivateNode(m.node, thisSearch);
      const bool moved = m.to != kInvalidPartition &&
        delta_phg.changeNodePart(m.node, m.from, m.to, maxPartWeight);
      if (moved) {
        runStats.moves++;
        insertOrUpdateNeighbors(phg, sharedData, m.node);
        estimatedImprovement += m.gain;
        localMoves.push_back(m);
        stopRule.update(m.gain);
        if (estimatedImprovement >= bestImprovement) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();
        }
      }
      updateBlocks(m.from);
    }

    applyMovesToGlobalHypergraph(phg, sharedData, bestImprovementIndex);
    runStats.estimated_improvement = bestImprovement;

    blockPQ.clear();
    for (PartitionID i = 0; i < numParts; ++i) {
      for (PosT j = 0; j < vertexPQs[i].size(); ++j) {
        sharedData.nodeTracker.releaseNode(vertexPQs[i].at(j));
      }
      vertexPQs[i].clear();
    }
    runStats.merge(stats);
  }

  void updateBlock(PartitionID i) {
    const bool underloaded = vertexPQs[i].empty() || delta_phg.partWeight(i) <= minPartWeight;
    if (!underloaded) {
      blockPQ.insertOrAdjustKey(i, vertexPQs[i].topKey());
    } else if (blockPQ.contains(i)) {
      blockPQ.remove(i);
    }
  }

  void updateBlocks(PartitionID moved_from) {
    if (moved_from == kInvalidPartition || updateDeduplicator.size() >= size_t(numParts)) {
      for (PartitionID i = 0; i < numParts; ++i) {
        updateBlock(i);
      }
    } else {
      updateBlock(moved_from);
      for (const HypernodeID v : updateDeduplicator.keys()) {
        updateBlock(delta_phg.partID(v));
      }
    }
    updateDeduplicator.clear();
  }

  void insertOrUpdateNeighbors(PartitionedHypergraph& phg, FMSharedData& sharedData, HypernodeID u) {
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

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> bestDestinationBlock(PartitionedHypergraph& phg, HypernodeID u) {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const PartitionID pu = delta_phg.partID(u);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
    for (PartitionID i = 0; i < numParts; ++i) {
      if (i != pu) {
        const HyperedgeWeight penalty = delta_phg.moveToPenalty(u, i);
        if (penalty < to_penalty && delta_phg.partWeight(i) + wu <= maxPartWeight) {
          to_penalty = penalty;
          to = i;
        }
      }
    }
    return std::make_pair(to, delta_phg.moveFromBenefit(u) - to_penalty);
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool insertOrUpdatePQ(PartitionedHypergraph& phg, HypernodeID v, NodeTracker& nt) {
    SearchID searchOfV = nt.searchOfNode[v].load(std::memory_order_acq_rel);
    // Note. Deactivated nodes have a special active search ID so that neither branch is executed
    if (nt.isSearchInactive(searchOfV)) {
      if (nt.searchOfNode[v].compare_exchange_strong(searchOfV, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID pv = delta_phg.partID(v);
        const Gain gain = bestDestinationBlock(phg, v).second;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
        runStats.pushes++;
        return true;
      }
    } else if (searchOfV == thisSearch) {
      const PartitionID pv = delta_phg.partID(v);
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

  bool findNextMove(PartitionedHypergraph& phg, Move& m) {
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

  void applyMovesToGlobalHypergraph(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex) {
    localAppliedMoves.clear();
    runStats.local_reverts += localMoves.size() - bestGainIndex;
    Gain estimatedImprovement = 0;
    Gain lastGain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      lastGain = km1Delta(he, edge_weight, edge_size,
        pin_count_in_from_part_after, pin_count_in_to_part_after);
      estimatedImprovement -= lastGain;
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain bestImprovement = 0;
    size_t bestIndex = 0;
    for ( size_t i = 0; i < bestGainIndex; ++i ) {
      Move& m = localMoves[i];
      MoveID m_id = std::numeric_limits<MoveID>::max();
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(),
        [&] { m_id = sharedData.moveTracker.insertMove(m); }, delta_gain_func);
      ASSERT(m_id != std::numeric_limits<MoveID>::max());
      Move& move = sharedData.moveTracker.getMove(m_id);
      move.gain = lastGain; // Update gain value based on hypergraph delta
      localAppliedMoves.push_back(m_id);
      if ( estimatedImprovement >= bestImprovement ) {
        bestImprovement = estimatedImprovement;
        bestIndex = i;
      }
    }

    // Kind of double rollback, if gain values are not correct
    ASSERT(localAppliedMoves.size() == bestGainIndex);
    for ( size_t i = bestIndex + 1; i < bestGainIndex; ++i ) {
      Move& m = sharedData.moveTracker.getMove(localAppliedMoves[i]);
      phg.changeNodePartFullUpdate(m.node, m.to, m.from,
        std::numeric_limits<HypernodeWeight>::max(), []{/* do nothing */});
      sharedData.moveTracker.invalidateMove(m);
    }
  }

  size_t numberOfSeedNodes(HypernodeID numNodes) {
    if (context.refinement.fm.use_seed_node_fraction) {
      const double x = context.refinement.fm.seed_node_fraction * numNodes / context.shared_memory.num_threads;
      return std::max(size_t(50), size_t(std::ceil(x)));
    } else {
      return context.refinement.fm.num_seed_nodes;
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
  // TODO if we're sharing PQ positions which causes HEAVY cache invalidation problems, why aren't we also sharing the bitset?
  // consider using cache-friendly hashmaps for heap positions?

  const Context& context;
  HypernodeWeight maxPartWeight = 0, perfectBalancePartWeight = 0, minPartWeight = 0;
  FMStats runStats;
  vec<Move> localMoves;
  vec<MoveID> localAppliedMoves;
public:
  FMStats stats;
};

}