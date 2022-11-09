/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <chrono>

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
