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


namespace mt_kahypar {
namespace refinement {

class LocalizedKWayFM {
public:
  // unfortunately the compiler thinks we're trying to pass a const-ref for the pq_handles, which we don't. therefore it had to be a pointer :(
  explicit LocalizedKWayFM(const Context& context, HypernodeID numNodes, PosT* pq_handles) :
          numParts(context.partition.k),
          blockPQ(static_cast<size_t>(numParts)),
          vertexPQs(static_cast<size_t>(numParts), VertexPriorityQueue(pq_handles, numNodes)),
          updateDeduplicator(numNodes),
          context(context),
          max_part_weight(context.partition.max_part_weights[0]),
          perfect_balance_part_weight(context.partition.perfect_balance_part_weights[0]),
          min_part_weight(static_cast<HypernodeWeight>(std::floor(perfect_balance_part_weight * (1 - context.partition.epsilon))))
  {

  }


  void findMoves(PartitionedHypergraph& phg, const HypernodeID initialBorderNode, FMSharedData& sharedData, SearchID search_id) {
    thisSearch = search_id;
    localMoves.clear();
    insertOrUpdatePQ(phg, initialBorderNode, sharedData.nodeTracker);
    updateBlock(phg, phg.partID(initialBorderNode));

    Move m;
    size_t consecutiveNonPositiveGainMoves = 0, consecutiveMovesWithNegativeOverallGain = 0, bestImprovementIndex = 0, queue_extractions = 0;
    Gain estimatedImprovement = 0, bestImprovement = 0;
    while (consecutiveNonPositiveGainMoves < context.refinement.fm.max_number_of_fruitless_moves && findNextMove(phg, m)) {
      ++queue_extractions;
      sharedData.nodeTracker.deactivateNode(m.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      const bool moved = m.to != kInvalidPartition
                         && phg.changeNodePartFullUpdate(m.node, m.from, m.to, max_part_weight, [&] { move_id = sharedData.moveTracker.insertMove(m); });
      if (moved) {
        updateAfterSuccessfulMove(phg, sharedData, m);
        estimatedImprovement += m.gain;
        localMoves.push_back(move_id);
        if (estimatedImprovement >= bestImprovement) {
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();
        }
        consecutiveNonPositiveGainMoves = m.gain > 0 ? 0 : consecutiveNonPositiveGainMoves + 1;
        consecutiveMovesWithNegativeOverallGain = estimatedImprovement >= 0 ? 0 : consecutiveMovesWithNegativeOverallGain + 1;

        LOG << V(m.gain) << V(estimatedImprovement) << V(bestImprovement);
      }
      updateAfterMoveExtraction(phg, m);
    }

    revertToBestLocalPrefix(phg, sharedData, bestImprovementIndex);

    blockPQ.clear();
    for (PartitionID i = 0; i < numParts; ++i) {
      vertexPQs[i].clear();
    }
    LOG << V(bestImprovement) << V(bestImprovementIndex) << V(consecutiveNonPositiveGainMoves) << V(queue_extractions);
  }

  void updateBlock(PartitionedHypergraph& phg, PartitionID i) {
    const bool underloaded = vertexPQs[i].empty() || phg.partWeight(i) <= min_part_weight;
    if (!underloaded) {
      blockPQ.insertOrAdjustKey(i, vertexPQs[i].topKey());
    } else if (blockPQ.contains(i)) {
      blockPQ.remove(i);
    }
  }

  void updateAfterMoveExtraction(PartitionedHypergraph& phg, Move& m) {
    if (updateDeduplicator.size() >= size_t(numParts)) {
      for (PartitionID i = 0; i < numParts; ++i) {
        updateBlock(phg, i);
      }
    } else {
      updateBlock(phg, m.from);
      for (const HypernodeID v : updateDeduplicator.keys()) {
        updateBlock(phg, phg.partID(v));
      }
    }
    updateDeduplicator.clear();
  }

  void updateAfterSuccessfulMove(PartitionedHypergraph& phg, FMSharedData& sharedData, Move& m) {
    const HypernodeID u = m.node;
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
    const PartitionID pu = phg.partID(u);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
    for (PartitionID i = 0; i < numParts; ++i) {
      if (i != pu) {
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if (penalty < to_penalty && phg.partWeight(i) + wu <= max_part_weight) {
          to_penalty = penalty;
          to = i;
        }
      }
    }
    return std::make_pair(to, phg.moveFromBenefit(u) - to_penalty);
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void insertOrUpdatePQ(PartitionedHypergraph& phg, HypernodeID v, NodeTracker& nt) {
    SearchID searchOfV = nt.searchOfNode[v].load(std::memory_order_acq_rel);
    // Note. Deactivated nodes have a special active search ID so that neither branch is executed
    if (nt.isSearchInactive(searchOfV)) {
      if (nt.searchOfNode[v].compare_exchange_strong(searchOfV, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID pv = phg.partID(v);
        const Gain gain = bestDestinationBlock(phg, v).second;
        vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
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
    }
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
        vertexPQs[from].deleteTop();  // blockPQ updates are done later, collectively.
        return true;
      } else {
        stats.retries++;
        vertexPQs[from].adjustKey(u, gain);
        if (vertexPQs[from].topKey() != blockPQ.keyOf(from)) {
          blockPQ.adjustKey(from, vertexPQs[from].topKey());
        }
      }
    }
  }

  void revertToBestLocalPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex) {
    while (localMoves.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localMoves.back());
      phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{/* do nothing */});
      sharedData.moveTracker.invalidateMove(m);
      localMoves.pop_back();
    }
  }

private:

  SearchID thisSearch;
  PartitionID numParts;

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
  HypernodeWeight max_part_weight, perfect_balance_part_weight, min_part_weight;

public:
  vec<MoveID> localMoves;

  FMStats stats;
};

}
}