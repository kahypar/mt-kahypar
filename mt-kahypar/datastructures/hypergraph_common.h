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

#include <cstdint>

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "tests/parallel/topology_mock.h"

#define USE_HARDWARE_MOCK false

namespace mt_kahypar {

#if USE_HARDWARE_MOCK
static constexpr int NUM_NUMA_NODES = 2;
using TopoMock = mt_kahypar::parallel::TopologyMock<NUM_NUMA_NODES>;
using topology_t = mt_kahypar::parallel::topology_t;
using node_t = mt_kahypar::parallel::node_t;
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
#else
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
#endif
using TBBNumaArena = mt_kahypar::parallel::TBBNumaArena<HardwareTopology, false>;

#define UI64(X) static_cast<uint64_t>(X)

using TaskGroupID = size_t;
using RatingType = double;
#if KAHYPAR_USE_64_BIT_IDS
#define ID(X) static_cast<uint64_t>(X)
using HypernodeID = uint64_t;
using HyperedgeID = uint64_t;
#else
#define ID(X) static_cast<uint32_t>(X)
using HypernodeID = uint32_t;
using HyperedgeID = uint32_t;
#endif
using HypernodeWeight = int32_t;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = HyperedgeWeight;

// Graph Types
using NodeID = uint32_t;
using ArcWeight = double;

struct Arc {
  NodeID head;
  ArcWeight weight;

  Arc() :
    head(0),
    weight(0) { }

  Arc(NodeID head, ArcWeight weight) :
    head(head),
    weight(weight) { }
};

// Constant Declarations
static constexpr PartitionID kInvalidPartition = -1;
static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
static constexpr HypernodeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
static constexpr size_t kEdgeHashSeed = 42;

static constexpr HypernodeID invalidNode = std::numeric_limits<HypernodeID>::max();
static constexpr Gain invalidGain = std::numeric_limits<Gain>::min();

struct Move {
  PartitionID from = -1;
  PartitionID to = -1;
  HypernodeID node = invalidNode;
  Gain gain = invalidGain;
};

struct Memento {
  HypernodeID u; // representative
  HypernodeID v; // contraction partner
};

using Batch = parallel::scalable_vector<Memento>;
using BatchVector = parallel::scalable_vector<Batch>;
using VersionedBatchVector = parallel::scalable_vector<BatchVector>;

using MoveID = uint32_t;
using SearchID = uint32_t;

struct NoOpDeltaFunc {
  void operator() (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }
};


/*!
  * This struct is used during multilevel coarsening to efficiently
  * detect parallel hyperedges.
  */
struct ContractedHyperedgeInformation {
  HyperedgeID he = kInvalidHyperedge;
  size_t hash = kEdgeHashSeed;
  size_t size = std::numeric_limits<size_t>::max();
  bool valid = false;
};

struct ParallelHyperedge {
  HyperedgeID removed_hyperedge;
  HyperedgeID representative;
};

// ! Helper function to compute delta for cut-metric after changeNodePart
static HyperedgeWeight cutDelta(const HyperedgeID,
                                const HyperedgeWeight edge_weight,
                                const HypernodeID edge_size,
                                const HypernodeID pin_count_in_from_part_after,
                                const HypernodeID pin_count_in_to_part_after) {
  if ( edge_size > 1 ) {
    if (pin_count_in_to_part_after == edge_size) {
      return -edge_weight;
    } else if (pin_count_in_from_part_after == edge_size - 1 &&
               pin_count_in_to_part_after == 1) {
      return edge_weight;
    }
  }
  return 0;
}

// ! Helper function to compute delta for km1-metric after changeNodePart
static HyperedgeWeight km1Delta(const HyperedgeID,
                                const HyperedgeWeight edge_weight,
                                const HypernodeID,
                                const HypernodeID pin_count_in_from_part_after,
                                const HypernodeID pin_count_in_to_part_after) {
  return (pin_count_in_to_part_after == 1 ? edge_weight : 0) +
         (pin_count_in_from_part_after == 0 ? -edge_weight : 0);
}

} // namespace mt_kahypar
