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