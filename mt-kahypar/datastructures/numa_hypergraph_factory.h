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

#include "tbb/parallel_for.h"

#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename Factory = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class NumaHypergraphFactory {

  using NumaHyperGraph = NumaHypergraph<Hypergraph, HardwareTopology, TBBNumaArena>;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  struct NumaMapping {
    parallel::scalable_vector<int> vertices_to_numa_node;
    parallel::scalable_vector<int> edges_to_numa_node;
  };

 public:
  static NumaHyperGraph construct(const HypernodeID num_hypernodes,
                                  const HyperedgeID num_hyperedges,
                                  const HyperedgeVector& edge_vector,
                                  const HyperedgeWeight* hyperedge_weight = nullptr,
                                  const HypernodeWeight* hypernode_weight = nullptr) {
    NumaHyperGraph hypergraph;

    // Compute mapping that maps each vertex and edge to a NUMA node
    NumaMapping numa_mapping = computeVertexAndEdgeToNumaNodeMapping(
      num_hypernodes, num_hyperedges, edge_vector);

    // Allocate empty hypergraphs on each NUMA node
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(
      TBBNumaArena::GLOBAL_TASK_GROUP, [&](const int) {
        hypergraph._hypergraphs.emplace_back();
      });

    // Construct hypergraphs on each NUMA node in parallel
    hypergraph._node_mapping.resize(num_hypernodes);
    hypergraph._edge_mapping.resize(num_hyperedges);
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      TBBNumaArena::GLOBAL_TASK_GROUP, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < hypergraph._hypergraphs.size());
        hypergraph._hypergraphs[node] = Factory::construct(
          node, num_hypernodes, num_hyperedges, edge_vector,
          numa_mapping.vertices_to_numa_node, numa_mapping.edges_to_numa_node,
          hypergraph._node_mapping, hypergraph._edge_mapping,
          hyperedge_weight, hypernode_weight);
      });

    // Setup internal stats of numa hypergraph
    for ( const Hypergraph& numa_hypergraph : hypergraph._hypergraphs ) {
      hypergraph._num_hypernodes += numa_hypergraph.initialNumNodes();
      hypergraph._num_hyperedges += numa_hypergraph.initialNumEdges();
      hypergraph._num_pins += numa_hypergraph.initialNumPins();
      hypergraph._total_weight += numa_hypergraph.totalWeight();
    }
    ASSERT(hypergraph.initialNumNodes() == num_hypernodes);
    ASSERT(hypergraph.initialNumEdges() == num_hyperedges);

    // Remap the vertex and edge ids of each numa hypergraph such that
    // the ids encodes the NUMA node the reside on and position inside
    // that hypergraph
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      TBBNumaArena::GLOBAL_TASK_GROUP, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < hypergraph._hypergraphs.size());
        Factory::remapVertexAndEdgeIds(
          hypergraph._hypergraphs[node],
          hypergraph._node_mapping,
          hypergraph._edge_mapping);
      });

    return hypergraph;
  }

 private:
  NumaHypergraphFactory() { }

  // ! Computes for each hyperedge and vertex a mapping to a numa node, to which
  // ! it will be assigned to.
  static NumaMapping computeVertexAndEdgeToNumaNodeMapping(const HypernodeID num_hypernodes,
                                                           const HyperedgeID num_hyperedges,
                                                           const HyperedgeVector& edge_vector) {
    parallel::scalable_vector<int> vertices_to_numa_node(num_hypernodes, 0);
    parallel::scalable_vector<int> edges_to_numa_node(num_hyperedges, 0);
    ASSERT(edge_vector.size() == num_hyperedges);
    const int used_numa_nodes = TBBNumaArena::instance().num_used_numa_nodes();
    ASSERT(used_numa_nodes > 0);
    const size_t num_hyperedges_per_hypergraph = num_hyperedges / used_numa_nodes;

    parallel::scalable_vector<AtomicCounter> num_hn_occurs_as_pin_on_numa_node(used_numa_nodes,
      AtomicCounter(num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0)));
    tbb::parallel_for(0UL, num_hyperedges, [&](const HyperedgeID he) {
      // Equally split the hyperedges across the numa nodes
      const HyperedgeID node = std::min(static_cast<int>(he / num_hyperedges_per_hypergraph), used_numa_nodes - 1);
      edges_to_numa_node[he] = node;
      for ( const HypernodeID& pin : edge_vector[he] ) {
        ASSERT(pin < num_hypernodes);
        ++num_hn_occurs_as_pin_on_numa_node[node][pin];
      }
    });

    // Compute vertex to numa node mapping based on
    // where each vertex occurs most as pin
    tbb::parallel_for(0UL, num_hypernodes, [&](const HypernodeID hn) {
      int best_count = std::numeric_limits<int>::min();
      int best_node = -1;
      for ( int node = 0; node < used_numa_nodes; ++node ) {
        int count = static_cast<int>(num_hn_occurs_as_pin_on_numa_node[node][hn]);
        if ( count > best_count ) {
          best_node = node;
          best_count = count;
        }
      }
      vertices_to_numa_node[hn] = best_node;
    });

    return NumaMapping { std::move(vertices_to_numa_node), std::move(edges_to_numa_node) };
  }


};

} // namespace ds
} // namespace mt_kahypar