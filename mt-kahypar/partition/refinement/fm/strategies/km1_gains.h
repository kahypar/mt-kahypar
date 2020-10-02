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

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {
struct Km1GainComputer {
  Km1GainComputer(PartitionID k) : gains(k, 0) { }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void computeGainsFromScratch(const PHG& phg, const HypernodeID u) {
    for (size_t i = 0; i < gains.size(); ++i) {
      gains[i] = 0;
    }

    const PartitionID from = phg.partID(u);
    Gain internal_weight = 0;   // weighted affinity with part(u) after u would be moved
    for (HyperedgeID e : phg.incidentEdges(u)) {
      HyperedgeWeight edge_weight = phg.edgeWeight(e);
      if (phg.pinCountInPart(e, from) > 1) {
        internal_weight += edge_weight;
      }

      if constexpr (PHG::supports_connectivity_set) {
        for (PartitionID i : phg.connectivitySet(e)) {
          gains[i] += edge_weight;
        }
      } else {
        // case for deltaPhg since maintaining connectivity sets is too slow
        for (size_t i = 0; i < gains.size(); ++i) {
          if (phg.pinCountInPart(e, i) > 0) {
            gains[i] += edge_weight;
          }
        }
      }
    }

    for (size_t i = 0; i < gains.size(); ++i) {
      gains[i] -= internal_weight;
    }
  }


  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PHG& phg,
                                                                 const HypernodeID u,
                                                                 const std::vector<HypernodeWeight>& max_part_weights) {
    computeGainsFromScratch(phg, u);

    const HypernodeWeight weight_of_u = phg.nodeWeight(u);
    const PartitionID from = phg.partID(u);
    PartitionID best_target = kInvalidPartition;
    HypernodeWeight best_target_weight = std::numeric_limits<HypernodeWeight>::max();
    Gain best_gain = std::numeric_limits<Gain>::min();
    for (PartitionID target = 0; target < phg.k(); ++target) {
      if (target != from) {
        const HypernodeWeight target_weight = phg.partWeight(target);
        const Gain gain = gains[target];
        if ( (gain > best_gain || (gain == best_gain && target_weight < best_target_weight))
             && target_weight + weight_of_u <= max_part_weights[target]) {
          best_target = target;
          best_gain = gain;
          best_target_weight = target_weight;
        }
      }
    }

    return std::make_pair(best_target, best_gain);
  }


  vec<Gain> gains;
};
}