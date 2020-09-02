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
            gain_tmp(sharedData.numParts, 0)
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
    std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PHG& phg,
                                                                   const HypernodeID u) {
      const PartitionID from = phg.partID(u);
      Gain internal_weight = 0;   // weighted affinity with part(u), even if u was moved
      for (HyperedgeID e : phg.incidentEdges(u)) {
        HyperedgeWeight edge_weight = phg.edgeWeight(e);
        if (phg.pinCountInPart(e, from) > 1) {
          internal_weight += edge_weight;
        }

        if constexpr (PHG::supports_connectivity_set) {
          for (PartitionID i : phg.connectivitySet(e)) {
            gain_tmp[i] += edge_weight;
          }
        } else {
          // case for deltaPhg since maintaining connectivity sets is too slow
          for (PartitionID i = 0; i < phg.k(); ++i) {
            if (phg.pinCountInPart(e, i) > 0) {
              gain_tmp[i] += edge_weight;
            }
          }
        }
      }

      const HypernodeWeight weight_of_u = phg.nodeWeight(u);
      PartitionID best_target = kInvalidPartition;
      HypernodeWeight best_target_weight = std::numeric_limits<HypernodeWeight>::max();
      Gain best_gain = std::numeric_limits<Gain>::min();
      for (PartitionID target = 0; target < phg.k(); ++target) {
        if (target != from) {
          const HypernodeWeight target_weight = phg.partWeight(target);
          Gain gain = gain_tmp[target] - internal_weight;
          if ( (gain > best_gain || (gain == best_gain && target_weight < best_target_weight))
                && target_weight + weight_of_u <= context.partition.max_part_weights[target]) {
            best_target = target;
            best_gain = gain;
            best_target_weight = target_weight;
          }
        }
        gain_tmp[target] = 0;
      }

      return std::make_pair(best_target, best_gain);
    }



  private:
    const Context& context;

    FMStats& runStats;

    FMSharedData& sharedData;

    VertexPriorityQueue pq;

    vec<Gain> gain_tmp;
  };


}