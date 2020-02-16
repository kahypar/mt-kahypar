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
#include "mt-kahypar/utils/timer.h"

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
  static NumaHyperGraph construct(const TaskGroupID task_group_id,
                                  const HypernodeID num_hypernodes,
                                  const HyperedgeID num_hyperedges,
                                  const HyperedgeVector& edge_vector,
                                  const HyperedgeWeight* hyperedge_weight = nullptr,
                                  const HypernodeWeight* hypernode_weight = nullptr) {
    NumaHyperGraph hypergraph;

    // Compute mapping that maps each vertex and edge to a NUMA node
    utils::Timer::instance().start_timer("compute_numa_mapping", "Compute Vertex to NUMA Mapping");
    NumaMapping numa_mapping = computeVertexAndEdgeToNumaNodeMapping(
      num_hypernodes, num_hyperedges, edge_vector);
    utils::Timer::instance().stop_timer("compute_numa_mapping");

    // Allocate empty hypergraphs on each NUMA node
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(
      task_group_id, [&](const int) {
        hypergraph._hypergraphs.emplace_back();
      });

    // Construct hypergraphs on each NUMA node in parallel
    utils::Timer::instance().start_timer("construct_numa_hypergraphs", "Construct NUMA hypergraphs");
    hypergraph._node_mapping.resize(num_hypernodes);
    hypergraph._edge_mapping.resize(num_hyperedges);
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < hypergraph._hypergraphs.size());
        hypergraph._hypergraphs[node] = Factory::construct(
          node, num_hypernodes, num_hyperedges, edge_vector,
          numa_mapping.vertices_to_numa_node, numa_mapping.edges_to_numa_node,
          hypergraph._node_mapping, hypergraph._edge_mapping,
          hyperedge_weight, hypernode_weight);
      });
    utils::Timer::instance().stop_timer("construct_numa_hypergraphs");

    // Setup internal stats of numa hypergraph
    for ( const Hypergraph& numa_hypergraph : hypergraph._hypergraphs ) {
      hypergraph._num_hypernodes += numa_hypergraph.initialNumNodes();
      hypergraph._num_hyperedges += numa_hypergraph.initialNumEdges();
      hypergraph._num_pins += numa_hypergraph.initialNumPins();
      hypergraph._total_degree += numa_hypergraph.initialTotalVertexDegree();
      hypergraph._total_weight += numa_hypergraph.totalWeight();
    }
    ASSERT(hypergraph.initialNumNodes() == num_hypernodes);
    ASSERT(hypergraph.initialNumEdges() == num_hyperedges);

    // Remap the vertex and edge ids of each numa hypergraph such that
    // the ids encodes the NUMA node the reside on and position inside
    // that hypergraph
    utils::Timer::instance().start_timer("compute_global_mappings", "Comp. Global HN and HE Mapping");
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < hypergraph._hypergraphs.size());
        Factory::remapVertexAndEdgeIds(
          hypergraph._hypergraphs[node],
          hypergraph._node_mapping,
          hypergraph._edge_mapping);
      });
    utils::Timer::instance().stop_timer("compute_global_mappings");

    return hypergraph;
  }

  static NumaHyperGraph construct(const TaskGroupID task_group_id,
                                  const NumaHyperGraph& hypergraph,
                                  parallel::scalable_vector<int>&& vertices_to_numa_node) {
    ASSERT(hypergraph.initialNumNodes() == vertices_to_numa_node.size());
    NumaHyperGraph redistributed_hypergraph;

    // Compute mapping that maps each vertex and edge to a NUMA node
    utils::Timer::instance().start_timer("compute_numa_mapping", "Compute Edge to NUMA Mapping");
    NumaMapping numa_mapping = computeEdgeToNumaNodeMapping(
      task_group_id, hypergraph, std::move(vertices_to_numa_node));
    utils::Timer::instance().stop_timer("compute_numa_mapping");

    utils::Timer::instance().start_timer("setup_construction_data", "Setup Hypergraph Construction Data");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
    HyperedgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;

    tbb::parallel_invoke([&] {
      hypernode_weight.assign(num_hypernodes, 1);
      hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
        const HypernodeID original_id = hypergraph.originalNodeID(hn);
        ASSERT(original_id < num_hypernodes);
        hypernode_weight[original_id] = hypergraph.nodeWeight(hn);
      });
    }, [&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.assign(num_hyperedges, 1);
      hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
        const HyperedgeID original_id = hypergraph.originalEdgeID(he);
        ASSERT(original_id < num_hyperedges);
        hyperedge_weight[original_id] = hypergraph.edgeWeight(he);
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          edge_vector[original_id].push_back(hypergraph.originalNodeID(pin));
        }
      });
    });
    utils::Timer::instance().stop_timer("setup_construction_data");

    // Allocate empty hypergraphs on each NUMA node
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(
      task_group_id, [&](const int) {
        redistributed_hypergraph._hypergraphs.emplace_back();
      });

    // Construct hypergraphs on each NUMA node in parallel
    utils::Timer::instance().start_timer("construct_numa_hypergraphs", "Construct NUMA hypergraphs");
    redistributed_hypergraph._node_mapping.resize(num_hypernodes);
    redistributed_hypergraph._edge_mapping.resize(num_hyperedges);
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < redistributed_hypergraph._hypergraphs.size());
        redistributed_hypergraph._hypergraphs[node] = Factory::construct(
          node, num_hypernodes, num_hyperedges, edge_vector,
          numa_mapping.vertices_to_numa_node, numa_mapping.edges_to_numa_node,
          redistributed_hypergraph._node_mapping, redistributed_hypergraph._edge_mapping,
          hyperedge_weight.data(), hypernode_weight.data());
      });
    utils::Timer::instance().stop_timer("construct_numa_hypergraphs");

    // Setup internal stats of numa hypergraph
    for ( const Hypergraph& numa_hypergraph : redistributed_hypergraph._hypergraphs ) {
      redistributed_hypergraph._num_hypernodes += numa_hypergraph.initialNumNodes();
      redistributed_hypergraph._num_hyperedges += numa_hypergraph.initialNumEdges();
      redistributed_hypergraph._num_pins += numa_hypergraph.initialNumPins();
      redistributed_hypergraph._total_degree += numa_hypergraph.initialTotalVertexDegree();
      redistributed_hypergraph._total_weight += numa_hypergraph.totalWeight();
    }
    ASSERT(redistributed_hypergraph.initialNumNodes() == num_hypernodes);
    ASSERT(redistributed_hypergraph.initialNumEdges() == num_hyperedges);

    // Remap the vertex and edge ids of each numa hypergraph such that
    // the ids encodes the NUMA node the reside on and position inside
    // that hypergraph
    utils::Timer::instance().start_timer("compute_global_mappings", "Comp. Global HN and HE Mapping");
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        ASSERT(static_cast<size_t>(node) < redistributed_hypergraph._hypergraphs.size());
        Factory::remapVertexAndEdgeIds(
          redistributed_hypergraph._hypergraphs[node],
          redistributed_hypergraph._node_mapping,
          redistributed_hypergraph._edge_mapping);
      });
    utils::Timer::instance().stop_timer("compute_global_mappings");

    return redistributed_hypergraph;
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

  // ! Computes for each hyperedge and vertex a mapping to a numa node, to which
  // ! it will be assigned to.
  static NumaMapping computeEdgeToNumaNodeMapping(const TaskGroupID task_group_id,
                                                  const NumaHyperGraph& hypergraph,
                                                  parallel::scalable_vector<int>&& vertices_to_numa_node) {
    const int used_numa_nodes = TBBNumaArena::instance().num_used_numa_nodes();
    parallel::scalable_vector<int> edges_to_numa_node(hypergraph.initialNumEdges(), 0);
    hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID& he) {
      parallel::scalable_vector<size_t> num_he_occurs_as_incident_net_on_numa_node(used_numa_nodes, 0);
      for ( const HypernodeID& pin : hypergraph.incidentEdges(he) ) {
        const int node = common::get_numa_node_of_vertex(pin);
        ++num_he_occurs_as_incident_net_on_numa_node[node];
      }

      int best_count = std::numeric_limits<int>::min();
      int best_node = -1;
      for ( int node = 0; node < used_numa_nodes; ++node ) {
        const int count = static_cast<int>(num_he_occurs_as_incident_net_on_numa_node[node]);
        if ( count > best_count ) {
          best_node = node;
          best_count = count;
        }
      }
      ASSERT(hypergraph.originalEdgeID(he) < edges_to_numa_node.size());
      edges_to_numa_node[hypergraph.originalEdgeID(he)] = best_node;
    });

    return NumaMapping { std::move(vertices_to_numa_node), std::move(edges_to_numa_node) };
  }


};

} // namespace ds
} // namespace mt_kahypar