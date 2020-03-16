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

#include <chrono>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/kway_priority_queue.h"

#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#if KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING
#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/datastructures/numa_hypergraph_factory.h"
#include "mt-kahypar/datastructures/numa_partitioned_hypergraph.h"
#endif

#include "tests/parallel/topology_mock.h"

#define USE_HARDWARE_MOCK false
#define KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING false

namespace mt_kahypar {

#if KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING
static constexpr int NUM_NUMA_NODES = 2;
#else
// In case non-numa-aware partitioning is enabled, we mock
// the system architecture to simulate a system with one
// NUMA node.
static constexpr int NUM_NUMA_NODES = 1;
#endif

#if USE_HARDWARE_MOCK || !KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING
using TopoMock = mt_kahypar::parallel::TopologyMock<NUM_NUMA_NODES>;
using topology_t = mt_kahypar::parallel::topology_t;
using node_t = mt_kahypar::parallel::node_t;
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
#else
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
#endif
using TBBNumaArena = mt_kahypar::parallel::TBBNumaArena<HardwareTopology>;

using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;
using KWayPriorityQueue = kahypar::ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>, true>;
using ThreadLocalKWayPriorityQueue = tbb::enumerable_thread_specific<KWayPriorityQueue>;

#if KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING
  using Hypergraph = ds::NumaHypergraph<ds::StaticHypergraph, HardwareTopology, TBBNumaArena>;
  using HypergraphFactory = ds::NumaHypergraphFactory<
    ds::StaticHypergraph, ds::StaticHypergraphFactory, HardwareTopology, TBBNumaArena>;
  using PartitionedHypergraph = ds::NumaPartitionedHypergraph<
    Hypergraph, HypergraphFactory>;
#else
  using Hypergraph = ds::StaticHypergraph;
  using HypergraphFactory = ds::StaticHypergraphFactory;
  using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, HypergraphFactory>;
#endif
  using Graph = ds::GraphT<Hypergraph>;

struct GlobalTypeTraits {
  using HyperGraph = Hypergraph;
  using HyperGraphFactory = HypergraphFactory;
  using PartitionedHyperGraph = PartitionedHypergraph;
  using TBB = TBBNumaArena;
  using HwTopology = HardwareTopology;
};

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
}  // namespace mt_kahypar
