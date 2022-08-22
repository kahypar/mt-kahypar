/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <chrono>
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_initializer.h"

/** NOTE
 * We have different configurations of our partitioner (resulting in several executables)
 * using different hypergraph data structures. The concrete implementation is set via compile
 * definitions in our CMake script and is configured in this file here. If you use a CMake version < 3.19
 * each executable will depend on all header files included in the following, since previous CMake versions do not evaluate
 * #ifdef definitions. This means that all executables are recompiled if you change one of these header files.
 * (even though the target does not depend on it). In later CMake versions (>= 3.19), the header dependencies are
 * extracted from the compiler, which evaluates the following #ifdef defintions correctly.
 * Therefore, we strongly recommend to use a CMake version >= 3.19 to reduce dependencies between
 * the different compile targets.
 */


#ifdef USE_GRAPH_PARTITIONER
#ifdef USE_STRONG_PARTITIONER
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#else
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_graph_factory.h"
#endif
#include "mt-kahypar/datastructures/partitioned_graph.h"
#include "mt-kahypar/datastructures/delta_partitioned_graph.h"
#else
#ifdef USE_STRONG_PARTITIONER
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"
#else
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#endif
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"
#endif

namespace mt_kahypar {

using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
using TBBInitializer = mt_kahypar::parallel::TBBInitializer<HardwareTopology, false>;

#ifdef USE_GRAPH_PARTITIONER
#ifdef USE_STRONG_PARTITIONER
using Hypergraph = ds::DynamicGraph;
using HypergraphFactory = ds::DynamicGraphFactory;
#else
using Hypergraph = ds::StaticGraph;
using HypergraphFactory = ds::StaticGraphFactory;
#endif
using PartitionedHypergraph = ds::PartitionedGraph<Hypergraph, HypergraphFactory>;
using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
#else
#ifdef USE_STRONG_PARTITIONER
using Hypergraph = ds::DynamicHypergraph;
using HypergraphFactory = ds::DynamicHypergraphFactory;
#else
using Hypergraph = ds::StaticHypergraph;
using HypergraphFactory = ds::StaticHypergraphFactory;
#endif
using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, HypergraphFactory>;
using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
#endif
#ifdef USE_GRAPH_PARTITIONER
using PartIdType = CAtomic<PartitionID>;
#else
using PartIdType = PartitionID;
#endif

using ParallelHyperedge = Hypergraph::ParallelHyperedge;

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

}  // namespace mt_kahypar
