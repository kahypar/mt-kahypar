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

#include "include/libmtkahypartypes.h"

#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_graph_factory.h"
#include "mt-kahypar/datastructures/partitioned_graph.h"
#include "mt-kahypar/datastructures/delta_partitioned_graph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"

namespace mt_kahypar {

struct StaticGraphTypeTraits {
  using Hypergraph = ds::StaticGraph;
  using PartitionedHypergraph = ds::PartitionedGraph<Hypergraph>;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};

struct DynamicGraphTypeTraits {
  using Hypergraph = ds::DynamicGraph;
  using PartitionedHypergraph = ds::PartitionedGraph<Hypergraph>;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};

struct StaticHypergraphTypeTraits {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, ds::ConnectivityInfo>;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

struct DynamicHypergraphTypeTraits {
  using Hypergraph = ds::DynamicHypergraph;
  using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, ds::ConnectivityInfo>;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

struct LargeKHypergraphTypeTraits {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, ds::SparseConnectivityInfo>;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

#define INSTANTIATE_FUNC_WITH_HYPERGRAPHS(FUNC) \
  template FUNC(ds::StaticHypergraph);          \
  template FUNC(ds::DynamicHypergraph);         \
  template FUNC(ds::StaticGraph);               \
  template FUNC(ds::DynamicGraph);

}  // namespace mt_kahypar
