/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {
struct Km1GainComputer {
  Km1GainComputer(PartitionID k) : gains(k, 0) { }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void computeGains(const PHG& phg, const HypernodeID u) {
    clear();
    Gain internal_weight = computeGainsPlusInternalWeight(phg, u);
    for (Gain& x : gains) { x -= internal_weight; }
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  Gain computeGainsPlusInternalWeight(const PHG& phg, const HypernodeID u) {
    assert(std::all_of(gains.begin(), gains.end(), [](const Gain& g) { return g == 0; }));
    const PartitionID from = phg.partID(u);
    Gain internal_weight = 0;   // weight that will not be removed from the objective
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
    return internal_weight;
  }

  void clear() {
    std::fill(gains.begin(), gains.end(), 0);
  }

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PHG& phg,
                                                                 const HypernodeID u,
                                                                 const std::vector<HypernodeWeight>& max_part_weights) {
    const HypernodeWeight weight_of_u = phg.nodeWeight(u);
    const PartitionID from = phg.partID(u);
    const Gain internal_weight = computeGainsPlusInternalWeight(phg, u);
    PartitionID best_target = kInvalidPartition;
    HypernodeWeight best_target_weight = std::numeric_limits<HypernodeWeight>::max();
    Gain best_gain = std::numeric_limits<Gain>::min();
    for (PartitionID target = 0; target < int(gains.size()); ++target) {
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
      gains[target] = 0;
    }

    best_gain -= internal_weight;
    return std::make_pair(best_target, best_gain);
  }

  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlockIgnoringBalance(const PartitionedHypergraph& phg,
                                                                                const HypernodeID u) {
    const PartitionID from = phg.partID(u);
    const Gain internal_weight = computeGainsPlusInternalWeight(phg, u);
    PartitionID best_target = kInvalidPartition;
    Gain best_gain = std::numeric_limits<Gain>::min();
    for (PartitionID target = 0; target < int(gains.size()); ++target) {
      if (target != from && gains[target] > best_gain) {
        best_gain = gains[target];
        best_target = target;
      }
      gains[target] = 0;
    }
    best_gain -= internal_weight;
    return std::make_pair(best_target, best_gain);
  }


  vec<Gain> gains;
};

struct TwoWayGainComputer {
  static Gain gainToOtherBlock(const PartitionedHypergraph& phg, const HypernodeID u) {
    Gain gain = 0;
    const PartitionID from = phg.partID(u);
    for (HyperedgeID e : phg.incidentEdges(u)) {
      const auto pcip = phg.pinCountInPart(e, from);
      const auto weight = phg.edgeWeight(e);
      if (pcip == 1) { gain += weight; }
      else if (pcip == phg.edgeSize(e)) { gain -= weight; }
    }
    return gain;
  }
};

}