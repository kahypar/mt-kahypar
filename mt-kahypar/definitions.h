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

#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"
#include "mt-kahypar/datastructures/streaming_hypergraph.h"


#include "tests/parallel/topology_mock.h"

namespace kahypar {

using TopoMock = kahypar::parallel::TopologyMock<2>;
using topology_t = kahypar::parallel::topology_t;
using node_t = kahypar::parallel::node_t;
using HardwareTopology = kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
// using HardwareTopology = kahypar::parallel::HardwareTopology<>;
using TBBNumaArena = kahypar::parallel::TBBNumaArena<HardwareTopology>;

using HypernodeID = uint32_t;
using HyperedgeID = uint32_t;
using HypernodeWeight = int32_t;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = HyperedgeWeight;

using StreamingHypergraph = kahypar::ds::StreamingHypergraph<HypernodeID, 
                                                             HyperedgeID, 
                                                             HypernodeWeight, 
                                                             HyperedgeWeight, 
                                                             PartitionID,
                                                             HardwareTopology,
                                                             TBBNumaArena>;

} // namespace kahypar