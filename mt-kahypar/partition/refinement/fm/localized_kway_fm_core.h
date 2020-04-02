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
          context(context)
  {

  }


  void findMoves(PartitionedHypergraph& phg, const HypernodeID u, FMSharedData& sharedData, SearchID search_id) {
    /*  NOTE (Lars): only for the version with local rollbacks
    HyperedgeWeight bestGain = 0;
    size_t bestGainIndex = 0;
    HyperedgeWeight overallGain = 0;
    */
    this->thisSearch = search_id;
    reinitialize();
    uint32_t movesWithNonPositiveGain = 0;
    insertOrUpdatePQ(phg, u, sharedData.nodeTracker);

    Move m;
    while (movesWithNonPositiveGain < context.refinement.fm.max_number_of_fruitless_moves && findNextMove(phg, m)) {

      sharedData.nodeTracker.deactivateNode(m.node, thisSearch);
      deactivatedNodes.push_back(m.node);

      if (phg.changeNodePartWithBalanceCheckAndGainUpdatesAndPartWeightUpdates(
              m.node, m.from, m.to, context.partition.max_part_weights[m.to])) {

        performSharedDataUpdates(m, phg, sharedData);
        movesWithNonPositiveGain = m.gain > 0 ? 0 : movesWithNonPositiveGain + 1;

        // activate neighbors of m.node and update their gains
        for (HyperedgeID e : phg.incidentEdges(m.node)) {
          if (phg.edgeSize(e) < context.partition.hyperedge_size_threshold) {
            for (HypernodeID v : phg.pins(e)) {
              if (!updateDeduplicator.contains(u)) {
                updateDeduplicator.insert(u);
                insertOrUpdatePQ(phg, v, sharedData.nodeTracker);
              }
            }
          }
        }
        updateDeduplicator.clear();

        /*  NOTE (Lars): only for the version with local rollbacks
        localMoves.push_back(m);
        overallGain += m.gain;
        if (overallGain < bestGain) {
          bestGain = overallGain;
          bestGainIndex = localMoves.size();
        }
        */
      }
    }

    // revertToBestLocalPrefix(phg, bestGainIndex);   NOTE (Lars): only for the version with local rollbacks

  }


private:

  void performSharedDataUpdates(Move& m, PartitionedHypergraph& phg, FMSharedData& sd) {
    const HypernodeID u = m.node;

    // TODO either sync this with the balance decision (requires double CAS operation)
    // or run another prefix sum on the final move order that checks balance
    // And at least get the move ID increment closer to the balance decision
    const MoveID move_id = sd.moveTracker.insertMove(m);

    for (HyperedgeID he : phg.incidentEdges(u)) {
      // update first move in
      std::atomic<MoveID>& fmi = sd.first_move_in[he * numParts + m.to];
      MoveID expected = fmi.load(std::memory_order_acq_rel);
      while ((sd.moveTracker.isIDStale(expected) || expected > move_id) && !fmi.compare_exchange_weak(expected, move_id, std::memory_order_acq_rel)) {  }

      // update last move out
      std::atomic<MoveID>& lmo = sd.last_move_out[he * numParts + m.from];
      expected = lmo.load(std::memory_order_acq_rel);
      while (expected < move_id && !lmo.compare_exchange_weak(expected, move_id, std::memory_order_acq_rel)) { }
    }
  }

  void insertOrUpdatePQ(PartitionedHypergraph& phg, HypernodeID u, NodeTracker& nt) {
    SearchID searchOfU = nt.searchOfNode[u].load(std::memory_order_acq_rel);

    if (nt.isSearchInactive(searchOfU)) {   // both branches already exclude deactivated nodes

      // try to claim node u
      if (nt.searchOfNode[u].compare_exchange_strong(searchOfU, thisSearch, std::memory_order_acq_rel)) {

        // if that was successful, we insert it into the vertices ready to move in its current block
        const PartitionID pu = phg.partID(u);
        auto [best_to_block, move_to_penalty] = phg.bestDestinationBlock(u);
        const Gain gain = phg.moveFromBenefit(u, pu) - move_to_penalty;
        if (vertexPQs[u].empty()) {
          blockPQ.insert(pu, gain);
        }
        vertexPQs[pu].insert(u, gain);

        // and we update the gain of moving the best vertex from that block
        if (blockPQ.contains(pu) && gain > blockPQ.keyOf(pu)) {
          blockPQ.increaseKey(pu, gain);
        }
      }

    } else if (searchOfU == thisSearch) {

      // update PQ entries
      const PartitionID pu = phg.partID(u);
      auto [best_to_block, move_to_penalty] = phg.bestDestinationBlock(u);
      const Gain gain = phg.moveFromBenefit(u, pu) - move_to_penalty;
      vertexPQs[pu].adjustKey(u, gain);

      if (blockPQ.contains(pu) && gain > blockPQ.keyOf(pu)) {
        blockPQ.increaseKey(pu, gain);
      }
    }

  }

  bool findNextMove(PartitionedHypergraph& phg, Move& m) {
    // TODO activate overloaded blocks and deactivate underloaded. this becomes substantially more expensive than in sequential
    // because sequentially this could be determined from just the last block that was moved from

    while (!blockPQ.empty()) {
      const PartitionID from = blockPQ.top();
      const HypernodeID u = vertexPQs[from].top();
      const HypernodeWeight wu = phg.nodeWeight(u);
      const Gain estimated_gain = vertexPQs[from].topKey();

      PartitionID to = kInvalidPartition;
      HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
      for (PartitionID i = 0; i < phg.k(); ++i) {
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if (penalty < to_penalty && phg.partWeight(to) + wu <= context.partition.max_part_weights[i]) {
          to_penalty = penalty;
          to = i;
        }
      }
      const Gain gain = phg.moveFromBenefit(u, from) - to_penalty;

      if (gain >= estimated_gain) {
        vertexPQs[from].deleteTop();
        m.node = u; m.to = to; m.from = from;
        m.gain = phg.km1Gain(u, from, to);
        return true;
      } else {
        vertexPQs[from].adjustKey(u, gain);
        if (vertexPQs[from].topKey() != blockPQ.keyOf(from)) {
          blockPQ.adjustKey(from, vertexPQs[from].topKey());
        }
      }
    }

    return false;
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


  SearchID thisSearch;
  vec<Move> localMoves;
  PartitionID numParts;

  BlockPriorityQueue blockPQ;
  vec<VertexPriorityQueue> vertexPQs;

  // Note: prefer ClearListSet over SparseSet because
  // ClearListSet takes numNodes + numInsertedNodes*32 bit
  // SparseSet takes 2 * numNodes * 32 bit
  // where numInsertedNodes is presumably much smaller than numNodes
  ldc::ClearListSet<HypernodeID> updateDeduplicator;

  const Context& context;

public:
  vec<HypernodeID> deactivatedNodes;
};

}
}