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