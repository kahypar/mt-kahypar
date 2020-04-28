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

// ! Selects the PQs in a round-robin fashion.
class RoundRobinPQSelectionPolicy {

 static HypernodeID kInvalidHypernode;
 static PartitionID kInvalidPartition;
 static Gain kInvalidGain;

 public:
  template<typename HyperGraph>
  static inline bool pop(const HyperGraph& hypergraph,
                         KWayPriorityQueue& pq,
                         HypernodeID& hn,
                         PartitionID& to,
                         Gain& gain,
                         const bool) {
    ASSERT(to >= kInvalidPartition && to < hypergraph.k());
    hn = kInvalidHypernode;
    gain = kInvalidGain;

    to = (to + 1) % hypergraph.k();
    const PartitionID start_block = to;
    while ( !pq.isEnabled(to) ) {
      to = (to + 1) % hypergraph.k();
      if ( start_block == to ) {
        to = kInvalidPartition;
        return false;
      }
    }

    ASSERT(to != kInvalidPartition && to < hypergraph.k());
    ASSERT(pq.isEnabled(to));
    pq.deleteMaxFromPartition(hn, gain, to);
    ASSERT(hn != kInvalidHypernode);
    return true;
  }

  // As default block we define the block to which all vertices are assigned to
  // before greedy initial partitioning. Experiments have shown that the greedy
  // round robin variant performs best if we leave all vertices unassigned before
  // greedy initial partitioning.
  static inline PartitionID getDefaultBlock() {
    return kInvalidPartition;
  }
};

HypernodeID RoundRobinPQSelectionPolicy::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
PartitionID RoundRobinPQSelectionPolicy::kInvalidPartition = -1;
Gain RoundRobinPQSelectionPolicy::kInvalidGain = std::numeric_limits<Gain>::min();

// ! Selects the PQ which contains the maximum gain move
class GlobalPQSelectionPolicy {

 static HypernodeID kInvalidHypernode;
 static PartitionID kInvalidPartition;
 static Gain kInvalidGain;

 public:
  template<typename HyperGraph>
  static inline bool pop(const HyperGraph&,
                         KWayPriorityQueue& pq,
                         HypernodeID& hn,
                         PartitionID& to,
                         Gain& gain,
                         const bool) {
    hn = kInvalidHypernode;
    to = kInvalidPartition;
    gain = kInvalidGain;

    if ( pq.numNonEmptyParts() > 0 && pq.numEnabledParts() > 0 ) {
      pq.deleteMax(hn, gain, to);
      ASSERT(hn != kInvalidHypernode);
      return true;
    } else {
      return false;
    }
  }

  // As default block we define the block to which all vertices are assigned to
  // before greedy initial partitioning. Experiments have shown that the greedy
  // global variant performs best if we assign all vertices to block 1 before
  // greedy initial partitioning.
  static inline PartitionID getDefaultBlock() {
    return 1;
  }
};

HypernodeID GlobalPQSelectionPolicy::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
PartitionID GlobalPQSelectionPolicy::kInvalidPartition = -1;
Gain GlobalPQSelectionPolicy::kInvalidGain = std::numeric_limits<Gain>::min();

// ! Selects the PQs one by one until they are disabled
class SequentialPQSelectionPolicy {

 static HypernodeID kInvalidHypernode;
 static PartitionID kInvalidPartition;
 static Gain kInvalidGain;

 public:
  template<typename HyperGraph>
  static inline bool pop(const HyperGraph& hypergraph,
                         KWayPriorityQueue& pq,
                         HypernodeID& hn,
                         PartitionID& to,
                         Gain& gain,
                         const bool use_perfect_balanced_as_upper_bound) {
    hn = kInvalidHypernode;
    gain = kInvalidGain;

    if ( use_perfect_balanced_as_upper_bound ) {
      if ( to == kInvalidPartition ) {
        to = 0;
      }

      while ( to < hypergraph.k() && !pq.isEnabled(to) ) {
        ++to;
      }

      if ( to < hypergraph.k() ) {
        ASSERT(pq.size(to) > 0);
        pq.deleteMaxFromPartition(hn, gain, to);
        ASSERT(hn != kInvalidHypernode);
        return true;
      } else {
        return false;
      }
    } else {
      return GlobalPQSelectionPolicy::pop(hypergraph,
        pq, hn, to, gain, use_perfect_balanced_as_upper_bound);
    }
  }

  // As default block we define the block to which all vertices are assigned to
  // before greedy initial partitioning. Experiments have shown that the greedy
  // sequential variant performs best if we assign all vertices to block 1 before
  // greedy initial partitioning.
  static inline PartitionID getDefaultBlock() {
    return 1;
  }
};

HypernodeID SequentialPQSelectionPolicy::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
PartitionID SequentialPQSelectionPolicy::kInvalidPartition = -1;
Gain SequentialPQSelectionPolicy::kInvalidGain = std::numeric_limits<Gain>::min();

} // namespace mt_kahypar
