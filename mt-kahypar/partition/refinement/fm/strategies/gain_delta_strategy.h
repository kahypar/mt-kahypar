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
  class GainDeltaStrategy {
  public:

    using BlockPriorityQueue = ds::ExclusiveHandleHeap<ds::MaxHeap<Gain, PartitionID>>;
    using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;

    static constexpr bool uses_gain_cache = false;

    GainDeltaStrategy(const Context& context,
                          HypernodeID numNodes,
                          FMSharedData& sharedData,
                          FMStats& runStats) :
            context(context),
            k(context.partition.k),
            runStats(runStats),
            sharedData(sharedData),
            blockPQ(sharedData.numParts),
            vertexPQs(),
            gc(sharedData.numParts)
    {
      vertexPQs.reserve(k);
      for (PartitionID i = 0; i < k; ++i) {
        vertexPQs.emplace_back(sharedData.vertexPQHandles.data() + i * k, numNodes);
      }
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void insertIntoPQ(const PHG& phg, const HypernodeID v) {
      gc.computeGainsFromScratch(phg, v);
      for (PartitionID i = 0; i < k; ++i) {
        if (i != phg.partID(v)) {
          vertexPQs[i].insert(v, gc.gains[i]);
        }
      }
      runStats.pushes++;
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateGain(const PHG& /*phg*/, const HypernodeID /*v*/, const Move& /*move*/) {
      // do nothing
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool findNextMove(const PHG& phg, Move& m) {
      if (blockPQ.empty()) {
        return false;
      }
      const PartitionID target = blockPQ.top();
      const HypernodeID u = vertexPQs[target].top();
      const Gain estimated_gain = vertexPQs[target].topKey();
      vertexPQs[target].deleteTop();
      ASSERT(estimated_gain == blockPQ.topKey());
      m.node = u; m.to = target; m.from = phg.partID(u);
      m.gain = estimated_gain;
      runStats.extractions++;
      for (PartitionID i = 0; i < k; ++i) {
        if (i != m.from && i != target) {
          vertexPQs[i].remove(u);
        }
      }
      return true;
    }

    void clearPQs(const size_t /* bestImprovementIndex */ ) {
      // release all nodes that were not moved
      const bool release = sharedData.release_nodes
                           && runStats.moves > 0;

      if (release) {
        // Release all nodes contained in the search

        for (PosT j = 0; j < vertexPQs[1].size(); ++j) {
          const HypernodeID v = vertexPQs[1].at(j);
          // we're not storing nodes in pqs for the block they're currently in --> have to check two pqs and deduplicate
          if (!vertexPQs[0].contains(v)) {
            sharedData.nodeTracker.releaseNode(v);
          }
        }
        for (PosT j = 0; j < vertexPQs[0].size(); ++j) {
          sharedData.nodeTracker.releaseNode(vertexPQs[0].at(j));
        }

      }

      for (PartitionID i = 0; i < k; ++i) {
        vertexPQs[i].clear();
      }
      blockPQ.clear();
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updatePQs(const PHG& phg) {
      for (PartitionID i = 0; i < k; ++i) {
        updateBlock(phg, i);
      }
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void deltaGainUpdates(PHG& phg, const HyperedgeID he, const HyperedgeWeight edge_weight,
                          const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                          const PartitionID to, const HypernodeID pin_count_in_to_part_after) {


      // perform delta gain updates for vertices that are in our search

      auto some_other_block = [&](const PartitionID i) {
        return k - 1 - i;
      };

      auto in_search = [&](const HypernodeID u) {
        return vertexPQs[some_other_block(phg.partID(u))].contains(u);
      };

      auto increase = [&](const HypernodeID u, const PartitionID i) {
        vertexPQs[i].increaseKey(u, vertexPQs[i].getKey(u) + edge_weight);
      };

      auto decrease = [&](const HypernodeID u, const PartitionID i) {
        vertexPQs[i].decreaseKey(u, vertexPQs[i].getKey(u) - edge_weight);
      };

      //  TODO this might be the place where a bucket pq works well

      // gain = moveFromBenefit - moveToPenalty

      if (pin_count_in_from_part_after == 1) {
        for (HypernodeID u : phg.pins(he)) {
          if (phg.partID(u) == from && in_search(u)) {
            // move from benefit increased --> gain increased
            for (PartitionID i = 0; i < k; ++i) {
              if (i != from) {
                increase(u, i);
              }
            }
          }
        }
      } else if (pin_count_in_from_part_after == 0) {
        for (HypernodeID u : phg.pins(he)) {
          // moveToPenalty increased --> gain decreased
          if (in_search(u)) {
            decrease(u, from);
          }
        }
      }

      if (pin_count_in_to_part_after == 1) {
        for (HypernodeID u : phg.pins(he)) {
          // moveToPenalty decreased --> gain increased
          if (in_search(u))  {
            increase(u, to);
          }
        }
      } else if (pin_count_in_to_part_after == 2) {
        for (HypernodeID u : phg.pins(he)) {
          if (phg.partID(u) == to && in_search(u)) {
            // move from benefit decreased --> gain decreased
            for (PartitionID i = 0; i < k; ++i) {
              if (i != to) {
                decrease(u, i);
              }
            }
          }
        }
      }


    }

  private:

    template<typename PHG>
    void updateBlock(const PHG& phg, PartitionID i) {
      if (vertexPQs[i].empty() || phg.partWeight(i) >= context.partition.max_part_weights[i]) {
        if (blockPQ.contains(i)) {
          blockPQ.remove(i);
        }
      } else {
        if (!blockPQ.contains(i)) {
          blockPQ.insert(i, vertexPQs[i].topKey());
        }
      }
    }

    size_t handle(HypernodeID u, PartitionID p) const {
      return size_t(u) * k +  p;
    }

    const Context& context;

    PartitionID k;

    FMStats& runStats;

    FMSharedData& sharedData;

    BlockPriorityQueue blockPQ;

    vec<VertexPriorityQueue> vertexPQs;

    Km1GainComputer gc;
  };


}