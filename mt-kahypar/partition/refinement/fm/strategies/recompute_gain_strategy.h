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

  /*
    TODO
    we're still mashing pq layout and gain strategy into one class.
    --> extract gain strategy into a template parameter for pq layout. problem is that gain_delta strategy needs access to the pq
   */

  class RecomputeGainStrategy {
  public:

    using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;

    static constexpr bool uses_gain_cache = false;

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
    void insertIntoPQ(const PHG& phg, const HypernodeID v) {
      auto [target, gain] = computeBestTargetBlock(phg, v);
      sharedData.targetPart[v] = target;
      pq.insert(v, gain);
      runStats.pushes++;
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateGain(const PHG& phg, const HypernodeID v, const Move& /*move*/) {
      ASSERT(pq.contains(v));
      auto [target, gain] = computeBestTargetBlock(phg, v);
      sharedData.targetPart[v] = target;
      pq.adjustKey(v, gain);
    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool findNextMove(const PHG& phg, Move& m) {
      if (pq.empty()) {
        return false;
      }
      while (true) {
        // TODO we're likely gonna have to kick the retries. they're too expensive and we want to 'just trust the gains'
        const HypernodeID u = pq.top();
        const Gain estimated_gain = pq.topKey();
        auto [to, gain] = computeBestTargetBlock(phg, u);
        if (gain >= estimated_gain) { // accept any gain that is at least as good
          m.node = u; m.to = to; m.from = phg.partID(u);
          m.gain = gain;
          runStats.extractions++;
          pq.deleteTop();
          return true;
        } else {
          runStats.retries++;
          pq.adjustKey(u, gain);
          sharedData.targetPart[u] = to;
        }
      }
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
    void updatePQs(const PHG& /* phg */) {

    }

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void deltaGainUpdates(PHG& , const HyperedgeID , const HyperedgeWeight ,
                          const PartitionID , const HypernodeID ,
                          const PartitionID , const HypernodeID ) {
      // do nothing!
    }

  private:

    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PHG& phg, const HypernodeID u) {
      gc.computeGainsFromScratch(phg, u);

      const HypernodeWeight weight_of_u = phg.nodeWeight(u);
      PartitionID best_target = kInvalidPartition;
      HypernodeWeight best_target_weight = std::numeric_limits<HypernodeWeight>::max();
      Gain best_gain = std::numeric_limits<Gain>::min();
      for (PartitionID target = 0; target < phg.k(); ++target) {
        if (target != from) {
          const HypernodeWeight target_weight = phg.partWeight(target);
          const Gain gain = gc.gains[target];
          if ( (gain > best_gain || (gain == best_gain && target_weight < best_target_weight))
                && target_weight + weight_of_u <= context.partition.max_part_weights[target]) {
            best_target = target;
            best_gain = gain;
            best_target_weight = target_weight;
          }
        }
      }

      return std::make_pair(best_target, best_gain);
    }



  private:
    const Context& context;

    FMStats& runStats;

    FMSharedData& sharedData;

    VertexPriorityQueue pq;

    Km1GainComputer gc;
  };


}