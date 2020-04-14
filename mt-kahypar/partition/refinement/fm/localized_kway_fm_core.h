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

#include "fm_commons.h"
#include "clearlist.hpp"


namespace mt_kahypar {
namespace refinement {


class LocalizedKWayFM {
public:
  // unfortunately the compiler thinks we're trying to pass a const-ref for the pq_handles, which we don't. therefore it had to be a pointer :(
  explicit LocalizedKWayFM(const Context& context, HypernodeID numNodes, vec<PosT>* pq_handles) :
          numParts(context.partition.k),
          blockPQ(static_cast<size_t>(numParts)),
          vertexPQs(static_cast<size_t>(numParts), VertexPriorityQueue(*pq_handles)),
          updateDeduplicator(numNodes),
          context(context),
          max_part_weight(context.partition.max_part_weights[0]),
          perfect_balance_part_weight(context.partition.perfect_balance_part_weights[0]),
          min_part_weight(static_cast<HypernodeWeight>(std::floor(perfect_balance_part_weight * (1 - context.partition.epsilon))))
  {

  }


  void findMoves(PartitionedHypergraph& phg, const HypernodeID initialBorderNode, FMSharedData& sharedData, SearchID search_id) {
    LOG << V(min_part_weight) << V(perfect_balance_part_weight) << V(max_part_weight) << V(phg.totalWeight()) << V(context.partition.epsilon);
    this->thisSearch = search_id;
    reinitialize();
    uint32_t movesWithNonPositiveGain = 0;
    insertOrUpdatePQ(phg, initialBorderNode, sharedData.nodeTracker);

    Move m;
    while (movesWithNonPositiveGain < context.refinement.fm.max_number_of_fruitless_moves && findNextMove(phg, m)) {
      MoveID move_id = 0;
      const bool moved = phg.changeNodePartFullUpdate(m.node, m.from, m.to, max_part_weight, [&] { move_id = sharedData.moveTracker.insertMove(m); });
      if (moved) {
        updateAfterSuccessfulMove(phg, sharedData, m, move_id);
        movesWithNonPositiveGain = m.gain > 0 ? 0 : movesWithNonPositiveGain + 1;
      }
      updateAfterMoveExtraction(phg, sharedData, m, moved);
    }
    LOG << V(movesWithNonPositiveGain) << V(context.refinement.fm.max_number_of_fruitless_moves);
  }

  void updateAfterMoveExtraction(PartitionedHypergraph& phg, FMSharedData& sharedData, Move& m, bool moved) {
    sharedData.nodeTracker.deactivateNode(m.node, thisSearch);
    deactivatedNodes.push_back(m.node);

    if (updateDeduplicator.size() >= numParts) {
      for (PartitionID i = 0; i < numParts; ++i) {
        if (blockPQ.contains(i)) {
          if (vertexPQs[i].empty() || phg.partWeight(i) <= min_part_weight) {
            blockPQ.remove(i);
          } else {
            blockPQ.adjustKey(i, vertexPQs[i].topKey());
          }
        } else if (!vertexPQs[i].empty() && phg.partWeight(i) > min_part_weight) {
          blockPQ.insert(i, vertexPQs[i].topKey());
        }
      }
    } else {
      // only update from PQ no matter whether the move happened or not.
      PartitionID from = m.from;
      if (vertexPQs[from].empty()) {
        if (blockPQ.contains(from)) { // the additional check is necessary, in case an earlier update already deleted the block
          blockPQ.remove(from);
        }
      } else {
        blockPQ.adjustKey(from, vertexPQs[from].topKey());
      }
    }

    updateDeduplicator.clear();
  }

  void updateAfterSuccessfulMove(PartitionedHypergraph& phg, FMSharedData& sharedData, Move& m, MoveID move_id) {
    const HypernodeID u = m.node;
    for (HyperedgeID e : phg.incidentEdges(u)) {
      sharedData.performHyperedgeSpecificMoveUpdates(m, move_id, e);
      // activate neighbors of u and update their gains
      if (true || phg.edgeSize(e) < context.partition.hyperedge_size_threshold) {
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
    return std::make_pair(to, phg.moveFromBenefit(u, phg.partID(u)) - to_penalty);
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void insertOrUpdatePQ(PartitionedHypergraph& phg, HypernodeID u, NodeTracker& nt) {
    SearchID searchOfU = nt.searchOfNode[u].load(std::memory_order_acq_rel);
    // Note. Deactivated nodes have a special inactive search ID

    if (nt.isSearchInactive(searchOfU)) {
      if (nt.searchOfNode[u].compare_exchange_strong(searchOfU, thisSearch, std::memory_order_acq_rel)) {
        const PartitionID from = phg.partID(u);
        auto [to, gain] = bestDestinationBlock(phg, u);
        vertexPQs[from].insert(u, gain);

        // optimization: skip blockPQ updates if more than numParts update calls were done.
        // in this case we do them once for each block after all calls to this function are done
        if (updateDeduplicator.size() < numParts) {
          if (!blockPQ.contains(from)) {
            if (phg.partWeight(from) > min_part_weight) {
              blockPQ.insert(from, gain);
              LOG << "insert" << V(from) << "into block PQ";
            }
          } else if (gain > blockPQ.keyOf(from)) {
            blockPQ.increaseKey(from, gain);
          }
        }
      }

    } else if (searchOfU == thisSearch) {
      const PartitionID from = phg.partID(u);
      // TODO is this call to bestDestinationBlock(..) really necessary?
      // it provides decent updates from other cores but seems somewhat expensive
      // maybe do it infrequently
      // can we prune updates?
      // only gains to from and to part can change --> it makes sense to store the best destination block per vertex and update it
      // or if phg.partID(u) == from/to part
      auto [to, gain] = bestDestinationBlock(phg, u);
      vertexPQs[from].adjustKey(u, gain);

      // optimization as above
      if (updateDeduplicator.size() < numParts) {
        if (blockPQ.contains(from) && gain > blockPQ.keyOf(from)) {
          blockPQ.increaseKey(from, gain);
        }
      }
    }

  }

  bool findNextMove(PartitionedHypergraph& phg, Move& m) {
    if (blockPQ.empty()) {
      LOG << "Block queue empty";
      for (PartitionID i = 0; i < numParts; ++i) {
        LOG << V(vertexPQs[i].size()) << V(phg.partWeight(i)) << V(min_part_weight);
      }
      return false;
    }
    while (true) {
      const PartitionID from = blockPQ.top();
      const HypernodeID u = vertexPQs[from].top();
      const Gain estimated_gain = vertexPQs[from].topKey();
      auto [to, gain] = bestDestinationBlock(phg, u);
      if (gain != estimated_gain) { LOG << "false estimate" << V(gain) << V(estimated_gain); }
      if (gain >= estimated_gain) { // accept any gain that is at least as good
        m.node = u; m.to = to; m.from = from;
        m.gain = gain;
        vertexPQs[from].deleteTop();  // update block PQ later due to potential updates
        return true;
      } else {
        vertexPQs[from].adjustKey(u, gain);
        if (vertexPQs[from].topKey() != blockPQ.keyOf(from)) {
          blockPQ.adjustKey(from, vertexPQs[from].topKey());
        }
      }
    }
  }

  void reinitialize() {
    localMoves.clear();
  }

  void revertToBestLocalPrefix(PartitionedHypergraph &phg, size_t bestGainIndex) {
    while (localMoves.size() > bestGainIndex) {
      Move& m = localMoves.back();
      phg.changeNodePart(m.node, m.to, m.from);
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

  const Context& context;
  HypernodeWeight max_part_weight, perfect_balance_part_weight, min_part_weight;


public:
  vec<Move> localMoves;
  vec<HypernodeID> deactivatedNodes;
};

}
}