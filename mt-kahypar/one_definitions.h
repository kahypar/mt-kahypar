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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

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

using StaticPartitionedGraph = ds::PartitionedGraph<ds::StaticGraph>;
using DynamicPartitionedGraph = ds::PartitionedGraph<ds::DynamicGraph>;
using StaticPartitionedHypergraph = ds::PartitionedHypergraph<ds::StaticHypergraph, ds::ConnectivityInfo>;
using StaticSparsePartitionedHypergraph = ds::PartitionedHypergraph<ds::StaticHypergraph, ds::SparseConnectivityInfo>;
using DynamicPartitionedHypergraph = ds::PartitionedHypergraph<ds::DynamicHypergraph, ds::ConnectivityInfo>;

struct StaticGraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticGraph;
  using PartitionedHypergraph = StaticPartitionedGraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};

struct DynamicGraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::DynamicGraph;
  using PartitionedHypergraph = DynamicPartitionedGraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};

struct StaticHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = StaticPartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

struct DynamicHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::DynamicHypergraph;
  using PartitionedHypergraph = DynamicPartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

struct LargeKHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = StaticSparsePartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

using TypeTraitsList = kahypar::meta::Typelist<StaticGraphTypeTraits,
                                               DynamicGraphTypeTraits,
                                               StaticHypergraphTypeTraits,
                                               DynamicHypergraphTypeTraits,
                                               LargeKHypergraphTypeTraits>;

#define INSTANTIATE_FUNC_WITH_HYPERGRAPHS(FUNC) \
  template FUNC(ds::StaticHypergraph);          \
  template FUNC(ds::DynamicHypergraph);         \
  template FUNC(ds::StaticGraph);               \
  template FUNC(ds::DynamicGraph);

#define INSTANTIATE_FUNC_WITH_PARTITIONED_HG(FUNC)  \
  template FUNC(StaticPartitionedGraph);            \
  template FUNC(DynamicPartitionedGraph);           \
  template FUNC(StaticPartitionedHypergraph);       \
  template FUNC(StaticSparsePartitionedHypergraph); \
  template FUNC(DynamicPartitionedHypergraph);

#define INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS(C) \
  template class C(StaticGraphTypeTraits);          \
  template class C(DynamicGraphTypeTraits);         \
  template class C(StaticHypergraphTypeTraits);     \
  template class C(DynamicHypergraphTypeTraits);    \
  template class C(LargeKHypergraphTypeTraits);

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS(C)       \
  template class C<StaticGraphTypeTraits>;          \
  template class C<DynamicGraphTypeTraits>;         \
  template class C<StaticHypergraphTypeTraits>;     \
  template class C<DynamicHypergraphTypeTraits>;    \
  template class C<LargeKHypergraphTypeTraits>;

#define INSTANTIATE_CLASS_WITH_HYPERGRAPHS(C) \
  template class C<ds::StaticHypergraph>;     \
  template class C<ds::DynamicHypergraph>;    \
  template class C<ds::StaticGraph>;          \
  template class C<ds::DynamicGraph>;         \

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

}  // namespace mt_kahypar
