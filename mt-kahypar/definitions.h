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

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"

namespace mt_kahypar {

using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;
using KWayPriorityQueue = kahypar::ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>, true>;
using ThreadLocalKWayPriorityQueue = tbb::enumerable_thread_specific<KWayPriorityQueue>;

#ifdef KAHYPAR_USE_N_LEVEL_PARADIGM
using Hypergraph = ds::DynamicHypergraph;
using HypergraphFactory = ds::DynamicHypergraphFactory;
#else
using Hypergraph = ds::StaticHypergraph;
using HypergraphFactory = ds::StaticHypergraphFactory;
#endif
using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, HypergraphFactory>;
using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
using Graph = ds::GraphT<Hypergraph>;

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

}  // namespace mt_kahypar
