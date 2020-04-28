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
#include "tests/parallel/topology_mock.h"

#define TRACK_BORDER_VERTICES true

namespace mt_kahypar {

static constexpr int NUM_NUMA_NODES = 1;

using TopoMock = mt_kahypar::parallel::TopologyMock<NUM_NUMA_NODES>;
using topology_t = mt_kahypar::parallel::topology_t;
using node_t = mt_kahypar::parallel::node_t;
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
using TBBNumaArena = mt_kahypar::parallel::TBBNumaArena<HardwareTopology>;

using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;
using KWayPriorityQueue = kahypar::ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>, true>;
using ThreadLocalKWayPriorityQueue = tbb::enumerable_thread_specific<KWayPriorityQueue>;

using Hypergraph = ds::StaticHypergraph;
using HypergraphFactory = ds::StaticHypergraphFactory;
template<bool track_border_vertices = TRACK_BORDER_VERTICES>
using PartitionedHypergraph = ds::PartitionedHypergraph<
  Hypergraph, HypergraphFactory, track_border_vertices>;
using Graph = ds::GraphT<Hypergraph>;

struct GlobalTypeTraits {
  using HyperGraph = Hypergraph;
  using HyperGraphFactory = HypergraphFactory;
  template<bool track_border_vertices = TRACK_BORDER_VERTICES>
  using PartitionedHyperGraph = PartitionedHypergraph<track_border_vertices>;
  using TBB = TBBNumaArena;
  using HwTopology = HardwareTopology;
};

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
}  // namespace mt_kahypar
