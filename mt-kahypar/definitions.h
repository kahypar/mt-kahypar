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

#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#endif
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_graph_factory.h"
#include "mt-kahypar/datastructures/partitioned_graph.h"
#include "mt-kahypar/datastructures/delta_partitioned_graph.h"
#endif
#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"
#endif
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"

namespace mt_kahypar {

#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
using StaticPartitionedGraph = ds::PartitionedGraph<ds::StaticGraph>;
#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
using DynamicPartitionedGraph = ds::PartitionedGraph<ds::DynamicGraph>;
#endif
#endif
using StaticPartitionedHypergraph = ds::PartitionedHypergraph<ds::StaticHypergraph, ds::ConnectivityInfo>;
#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
using DynamicPartitionedHypergraph = ds::PartitionedHypergraph<ds::DynamicHypergraph, ds::ConnectivityInfo>;
#endif
#ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
using StaticSparsePartitionedHypergraph = ds::PartitionedHypergraph<ds::StaticHypergraph, ds::SparseConnectivityInfo>;
#endif

#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
struct StaticGraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticGraph;
  using PartitionedHypergraph = StaticPartitionedGraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};

#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
struct DynamicGraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::DynamicGraph;
  using PartitionedHypergraph = DynamicPartitionedGraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedGraph<PartitionedHypergraph>;
};
#endif
#endif

struct StaticHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = StaticPartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};

#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
struct DynamicHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::DynamicHypergraph;
  using PartitionedHypergraph = DynamicPartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};
#endif

#ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
struct LargeKHypergraphTypeTraits : public kahypar::meta::PolicyBase {
  using Hypergraph = ds::StaticHypergraph;
  using PartitionedHypergraph = StaticSparsePartitionedHypergraph;
  using DeltaPartitionedHypergraph = ds::DeltaPartitionedHypergraph<PartitionedHypergraph>;
};
#endif

#ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
#define ENABLE_LARGE_K(X) X
#else
#define ENABLE_LARGE_K(X)
#endif

#ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
#define ENABLE_N_LEVEL(X) X
#else
#define ENABLE_N_LEVEL(X)
#endif

#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
#define ENABLE_GRAPHS(X) X
#else
#define ENABLE_GRAPHS(X)
#endif

#define COMMA ,
using TypeTraitsList = kahypar::meta::Typelist<StaticHypergraphTypeTraits
                                               ENABLE_N_LEVEL(COMMA DynamicHypergraphTypeTraits)
                                               ENABLE_LARGE_K(COMMA LargeKHypergraphTypeTraits)
                                               ENABLE_GRAPHS(COMMA StaticGraphTypeTraits)
                                               ENABLE_GRAPHS(ENABLE_N_LEVEL(COMMA DynamicGraphTypeTraits))>;

#define INSTANTIATE_FUNC_WITH_HYPERGRAPHS(FUNC)                    \
  template FUNC(ds::StaticHypergraph);                             \
  ENABLE_N_LEVEL(template FUNC(ds::DynamicHypergraph);)            \
  ENABLE_GRAPHS(template FUNC(ds::StaticGraph);)                   \
  ENABLE_GRAPHS(ENABLE_N_LEVEL(template FUNC(ds::DynamicGraph);))

#define INSTANTIATE_CLASS_WITH_HYPERGRAPHS(C)                         \
  template class C<ds::StaticHypergraph>;                             \
  ENABLE_N_LEVEL(template class C<ds::DynamicHypergraph>;)            \
  ENABLE_GRAPHS(template class C<ds::StaticGraph>;)                   \
  ENABLE_GRAPHS(ENABLE_N_LEVEL(template class C<ds::DynamicGraph>;))

#define INSTANTIATE_FUNC_WITH_PARTITIONED_HG(FUNC)                       \
  template FUNC(StaticPartitionedHypergraph);                            \
  ENABLE_LARGE_K(template FUNC(StaticSparsePartitionedHypergraph);)      \
  ENABLE_N_LEVEL(template FUNC(DynamicPartitionedHypergraph);)           \
  ENABLE_GRAPHS(template FUNC(StaticPartitionedGraph);)                  \
  ENABLE_GRAPHS(ENABLE_N_LEVEL(template FUNC(DynamicPartitionedGraph);))

#define INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS(C)                          \
  template class C(StaticHypergraphTypeTraits);                              \
  ENABLE_N_LEVEL(template class C(DynamicHypergraphTypeTraits);)             \
  ENABLE_LARGE_K(template class C(LargeKHypergraphTypeTraits);)              \
  ENABLE_GRAPHS(template class C(StaticGraphTypeTraits);)                    \
  ENABLE_GRAPHS(ENABLE_N_LEVEL(template class C(DynamicGraphTypeTraits);))

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS(C)                               \
  template class C<StaticHypergraphTypeTraits>;                             \
  ENABLE_N_LEVEL(template class C<DynamicHypergraphTypeTraits>;)            \
  ENABLE_LARGE_K(template class C<LargeKHypergraphTypeTraits>;)             \
  ENABLE_GRAPHS(template class C<StaticGraphTypeTraits>;)                   \
  ENABLE_LARGE_K(ENABLE_N_LEVEL(template class C<DynamicGraphTypeTraits>;))

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

}  // namespace mt_kahypar
