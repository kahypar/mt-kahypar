/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

class CutGainPolicy {

 static constexpr bool enable_heavy_assert = false;
 static PartitionID kInvalidPartition;

 public:
  template<typename HyperGraph>
  static inline Gain calculateGain(const HyperGraph& hypergraph,
                                   const HypernodeID hn,
                                   const PartitionID to) {
    if ( hypergraph.partID(hn) == kInvalidPartition ) {
      return calculateGainForInvalidBlock(hypergraph, hn, to);
    } else {
      return calculateGainForValidBlock(hypergraph, hn, to);
    }
  }

  template<typename HyperGraph>
  static inline Gain calculateGainForInvalidBlock(const HyperGraph& hypergraph,
                                                  const HypernodeID hn,
                                                  const PartitionID to) {
    ASSERT(hypergraph.partID(hn) == kInvalidPartition);
    ASSERT(to != kInvalidPartition);

    Gain gain = 0;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      if ( hypergraph.connectivity(he) == 1 && hypergraph.pinCountInPart(he, to) == 0 ) {
        gain -= hypergraph.edgeWeight(he);
      }
    }
    return gain;
  }

  template<typename HyperGraph>
  static inline Gain calculateGainForValidBlock(const HyperGraph& hypergraph,
                                                const HypernodeID hn,
                                                const PartitionID to) {
    ASSERT(hypergraph.partID(hn) != kInvalidPartition);
    ASSERT(hypergraph.partID(hn) != to);
    ASSERT(to != kInvalidPartition);

    Gain gain = 0;
    const PartitionID from = hypergraph.partID(hn);
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      if ( hypergraph.edgeSize(he) > 1 ) {
        const PartitionID connectivity = hypergraph.connectivity(he);
        const HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
        const HypernodeID pin_count_in_to_part = hypergraph.pinCountInPart(he, to);
        if ( connectivity == 1 &&
             pin_count_in_from_part > 1 ) {
          // In case connectivity is one and there is more than one pin left
          // in from part, we would make the hyperedge cut if move hn to
          // block to
          gain -= hypergraph.edgeWeight(he);
        } else if ( connectivity == 2 &&
                    pin_count_in_from_part == 1 &&
                    pin_count_in_to_part > 0 ) {
          // In case, the connectivity is two and hn is the last pin left
          // of block from in hyperedge he, we would make the hyperedge a
          // non-cut hyperedge by moving hn to block to.
          gain += hypergraph.edgeWeight(he);
        }
      }
    }
    return gain;
  }

  template<typename HyperGraph>
  static inline void deltaGainUpdate(const HyperGraph& hypergraph,
                                     KWayPriorityQueue& pq,
                                     const HypernodeID hn,
                                     const PartitionID from,
                                     const PartitionID to) {
    if ( from == kInvalidPartition ) {
      deltaGainUpdateForInvalidBlock(hypergraph, pq, hn, from, to);
    } else {
      deltaGainUpdateForValidBlock(hypergraph, pq, hn, from, to);
    }

    HEAVY_INITIAL_PARTITIONING_ASSERT([&]() {
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (pin != hn) {
              const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
              for (PartitionID block = 0; block < hypergraph.k(); ++block) {
                if (pq.contains(original_pin_id, block)) {
                  const Gain gain = calculateGain(hypergraph, pin, block);
                  if (pq.key(original_pin_id, block) != gain) {
                    LOG << V(hn);
                    LOG << V(to);
                    LOG << V(he);
                    LOG << V(pin);
                    LOG << V(block);
                    LOG << V(gain);
                    LOG << V(pq.key(original_pin_id, block));
                    return false;
                  }
                }
              }
            }
          }
        }
        return true;
      } (), "Delta Gain Update failed!");
  }

  template<typename HyperGraph>
  static inline void deltaGainUpdateForInvalidBlock(const HyperGraph& hypergraph,
                                                    KWayPriorityQueue& pq,
                                                    const HypernodeID hn,
                                                    const PartitionID,
                                                    const PartitionID to) {
    ASSERT(hypergraph.partID(hn) == to);
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      const HypernodeID pin_count_in_to_part_after = hypergraph.pinCountInPart(he, to);
      const PartitionID connectivity = hypergraph.connectivity(he);
      const HyperedgeWeight he_weight = hypergraph.edgeWeight(he);

      if ( pin_count_in_to_part_after == 1 ) {
        if ( connectivity == 1 ) {
          // Connectivity changed from 0 to 1 => Each move to an other block (except to)
          // would make hyperedge he cut
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
            for ( PartitionID block = 0; block < hypergraph.k(); ++block ) {
              if ( pin != hn && block != to && pq.contains(original_pin_id, block) ) {
                pq.updateKeyBy(original_pin_id, block, -he_weight);
              }
            }
          }
        } else if ( connectivity == 2 ) {
          // Connectivity changed from 1 to 2 => Each move to a block that is not part
          // of the connectivity set of the hyperedge does not increase the cut any more.
          for ( PartitionID block = 0; block < hypergraph.k(); ++block ) {
            if ( block == to || hypergraph.pinCountInPart(he, block) == 0 ) {
              for ( const HypernodeID& pin : hypergraph.pins(he) ) {
                const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
                if ( pin != hn && pq.contains(original_pin_id, block) ) {
                  pq.updateKeyBy(original_pin_id, block, he_weight);
                }
              }
            }
          }
        }
      }
    }
  }

  template<typename HyperGraph>
  static inline void deltaGainUpdateForValidBlock(const HyperGraph& hypergraph,
                                                  KWayPriorityQueue& pq,
                                                  const HypernodeID hn,
                                                  const PartitionID from,
                                                  const PartitionID to) {
    ASSERT(hypergraph.partID(hn) == to);
    ASSERT(from != kInvalidPartition);
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      const HypernodeID he_size = hypergraph.edgeSize(he);

      if ( he_size > 1 ) {
        const HypernodeID pin_count_in_from_part_before = hypergraph.pinCountInPart(he, from) + 1;
        const HypernodeID pin_count_in_to_part_after = hypergraph.pinCountInPart(he, to);
        const HyperedgeWeight he_weight = hypergraph.edgeWeight(he);

        if ( pin_count_in_from_part_before == he_size ) {
          ASSERT(hypergraph.connectivity(he) == 2);
          ASSERT(pin_count_in_to_part_after == 1);
          // In case, the pin count in hyperedge he of block from was equal to the
          // hyperedge size before, we have made the hyperedge a cut hyperedge.
          // All moves to a different block than from do not increase the cut any
          // more. Therefore, we increase the gain of all pins that are contained
          // in PQs different from block from.
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
            for ( PartitionID block = 0; block < hypergraph.k(); ++block ) {
              if ( pin != hn && block != from && pq.contains(original_pin_id, block) ) {
                pq.updateKeyBy(original_pin_id, block, he_weight);
              }
            }
          }
        }

        if ( pin_count_in_to_part_after == he_size - 1 ) {
          // In case, the pin count in hyperedge he of block to is equal to the
          // hyperedge size minus one, we could make the hyperedge a non-cut
          // hyperedge by moving the only pin left in the from block to block to.
          // Therefore, we increase the gain of that pin.
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
            if ( pin != hn && hypergraph.partID(pin) != to && pq.contains(original_pin_id, to) ) {
              pq.updateKeyBy(original_pin_id, to, he_weight);
            }
          }
        }
      }
    }
  }
};

PartitionID CutGainPolicy::kInvalidPartition = -1;


class MaxNetGainPolicy {

 static constexpr bool enable_heavy_assert = false;
 static PartitionID kInvalidPartition;

 public:
  template<typename HyperGraph>
  static inline Gain calculateGain(const HyperGraph& hypergraph,
                                   const HypernodeID hn,
                                   const PartitionID to) {
    ASSERT(to != kInvalidPartition);

    Gain gain = 0;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      if (hypergraph.pinCountInPart(he, to) > 0) {
        gain += hypergraph.edgeWeight(he);
      }
    }
    return gain;
  }

  template<typename HyperGraph>
  static inline void deltaGainUpdate(const HyperGraph& hypergraph,
                                     KWayPriorityQueue& pq,
                                     const HypernodeID hn,
                                     const PartitionID,
                                     const PartitionID to) {
    ASSERT(hypergraph.partID(hn) == to);
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      const HyperedgeWeight he_weight = hypergraph.edgeWeight(he);
      const HypernodeID pins_in_to_part = hypergraph.pinCountInPart(he, to);
      if ( pins_in_to_part == 1 ) {
        // Block to was not part of hyperedge he before
        // => Update gain of all pins in hyperedge to block to
        for ( const HypernodeID pin : hypergraph.pins(he) ) {
          const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
          if ( pq.contains(original_pin_id, to) ) {
            pq.updateKeyBy(original_pin_id, to, he_weight);
          }
        }
      }
    }

    HEAVY_INITIAL_PARTITIONING_ASSERT([&]() {
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (pin != hn) {
              const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
              for (PartitionID block = 0; block < hypergraph.k(); ++block) {
                if (pq.contains(original_pin_id, block)) {
                  const Gain gain = calculateGain(hypergraph, pin, block);
                  if (pq.key(original_pin_id, block) != gain) {
                    LOG << V(hn);
                    LOG << V(to);
                    LOG << V(he);
                    LOG << V(pin);
                    LOG << V(block);
                    LOG << V(gain);
                    LOG << V(pq.key(original_pin_id, block));
                    return false;
                  }
                }
              }
            }
          }
        }
        return true;
      } (), "Delta Gain Update failed!");
  }
};

PartitionID MaxNetGainPolicy::kInvalidPartition = -1;

} // namespace mt_kahypar
