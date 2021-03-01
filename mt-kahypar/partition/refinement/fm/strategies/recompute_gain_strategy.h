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

#include "km1_gains.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"


namespace mt_kahypar {

  class RecomputeGainStrategy {
  public:

    using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;

    static constexpr bool uses_gain_cache = false;
    static constexpr bool maintain_gain_cache_between_rounds = false;

    RecomputeGainStrategy(const Context& context,
                      HypernodeID numNodes,
                      FMSharedData& sharedData,
                      FMStats& runStats) :
            context(context),
            runStats(runStats),
            sharedData(sharedData),
            pq(VertexPriorityQueue(sharedData.vertexPQHandles.data(), numNodes)),
            gc(sharedData.numParts)
    { }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void insertIntoPQ(const PHG& phg, const HypernodeID v, const SearchID ) {
      auto [target, gain] = gc.computeBestTargetBlock(phg, v, context.partition.max_part_weights);
      sharedData.targetPart[v] = target;
      pq.insert(v, gain);
      runStats.pushes++;
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateGain(const PHG& phg, const HypernodeID v, const Move& /*move*/) {
      ASSERT(pq.contains(v));
      auto [target, gain] = gc.computeBestTargetBlock(phg, v, context.partition.max_part_weights);
      sharedData.targetPart[v] = target;
      pq.adjustKey(v, gain);
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool findNextMove(const PHG& phg, Move& m) {
      if (pq.empty()) {
        return false;
      }
      const HypernodeID u = pq.top();
      auto [to, gain] = gc.computeBestTargetBlock(phg, u, context.partition.max_part_weights);
      m.node = u; m.to = to; m.from = phg.partID(u);
      m.gain = gain;
      runStats.extractions++;
      pq.deleteTop();
      return true;
    }

    void clearPQs(const size_t /* bestImprovementIndex */ ) {
      // release all nodes that were not moved
      const bool release = sharedData.release_nodes
                           && runStats.moves > 0;

      if (release) {
        // Release all nodes contained in PQ
        for (PosT j = 0; j < pq.size(); ++j) {
          const HypernodeID v = pq.at(j);
          sharedData.nodeTracker.releaseNode(v);
        }
      }

      pq.clear();
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void deltaGainUpdates(PHG& , const HyperedgeID , const HyperedgeWeight ,
                          const PartitionID , const HypernodeID ,
                          const PartitionID , const HypernodeID ) {
      // do nothing!
    }

    void memoryConsumption(utils::MemoryTreeNode *parent) const {
      parent->addChild("PQs", pq.size_in_bytes());
      parent->addChild("Initial Gain Comp", gc.gains.size() * sizeof(Gain));
    }

  private:
    const Context& context;

    FMStats& runStats;

    FMSharedData& sharedData;

    VertexPriorityQueue pq;

    Km1GainComputer gc;
  };


}