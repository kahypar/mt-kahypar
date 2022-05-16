/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "static_graph_factory.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {
  void StaticGraphFactory::sort_incident_edges(StaticGraph& graph) {
    parallel::scalable_vector<HyperedgeID> edge_ids_of_node;
    edge_ids_of_node.resize(graph._edges.size());
    // sort incident edges of each node, so their ordering is independent of scheduling
    // (and the same as a typical sequential implementation)
    tbb::parallel_for(ID(0), graph._num_nodes, [&](HypernodeID u) {
      const HyperedgeID start = graph.node(u).firstEntry();
      const HyperedgeID end = graph.node(u + 1).firstEntry();
      for (HyperedgeID id = 0; id < end - start; ++id) {
        edge_ids_of_node[start + id] = id;
      }
      std::sort(edge_ids_of_node.begin() + start, edge_ids_of_node.begin() + end, [&](HyperedgeID& a, HyperedgeID& b) {
        return graph.edge(start + a).target() < graph.edge(start + b).target();
      });

      // apply permutation
      // (yes, this applies the permutation defined by edge_ids_of_node, don't think about it)
      for (size_t i = 0; i < end - start; ++i) {
        HyperedgeID target = edge_ids_of_node[start + i];
        while (target < i) {
          target = edge_ids_of_node[start + target];
        }
        std::swap(graph._edges[start + i], graph._edges[start + target]);
        std::swap(graph._unique_edge_ids[start + i], graph._unique_edge_ids[start + target]);
      }
    });
  }

  StaticGraph StaticGraphFactory::construct(
          const HypernodeID num_nodes,
          const HyperedgeID num_edges,
          const HyperedgeVector& edge_vector,
          const HyperedgeWeight* edge_weight,
          const HypernodeWeight* node_weight,
          const bool stable_construction_of_incident_edges) {
    ASSERT(edge_vector.size() == num_edges);

    EdgeVector edges;
    edges.reserve(num_edges);
    for (const auto& e : edge_vector) {
      if (e.size() != 2) {
        ERROR("Using graph data structure; but the input hypergraph is not a graph.");
      }
      edges.push_back({e[0], e[1]});
    }
    return construct_from_graph_edges(num_nodes, num_edges, edges,
                                      edge_weight, node_weight,
                                      stable_construction_of_incident_edges);
  }

  StaticGraph StaticGraphFactory::construct_from_graph_edges(
          const HypernodeID num_nodes,
          const HyperedgeID num_edges,
          const EdgeVector& edge_vector,
          const HyperedgeWeight* edge_weight,
          const HypernodeWeight* node_weight,
          const bool stable_construction_of_incident_edges) {
    StaticGraph graph;
    graph._num_nodes = num_nodes;
    graph._num_edges = 2 * num_edges;
    graph._nodes.resize(num_nodes + 1);
    graph._edges.resize(2 * num_edges);
    graph._unique_edge_ids.resize(2 * num_edges);

    ASSERT(edge_vector.size() == num_edges);

    // Compute degree for each vertex
    utils::Timer::instance().start_timer("compute_ds_sizes", "Precompute DS Size", true);
    ThreadLocalCounter local_degree_per_vertex(num_nodes);
    ThreadLocalCounter local_incident_weight_per_vertex(edge_weight ? num_nodes : 0);
    tbb::parallel_for(ID(0), num_edges, [&](const size_t pos) {
      Counter& num_degree_per_vertex = local_degree_per_vertex.local();
      Counter& num_incident_weight_per_vertex = local_incident_weight_per_vertex.local();
      const HypernodeID pins[2] = {edge_vector[pos].first, edge_vector[pos].second};
      for (const HypernodeID& pin : pins) {
        ASSERT(pin < num_nodes, V(pin) << V(num_nodes));
        ++num_degree_per_vertex[pin];
        if (edge_weight) {
          num_incident_weight_per_vertex[pin] += edge_weight[pos];
        }
      }
    });

    // We sum up the degree per vertex only thread local. To obtain the
    // global degree, we iterate over each thread local counter and sum it up.
    Counter num_degree_per_vertex(num_nodes, 0);
    for (Counter& c : local_degree_per_vertex) {
      tbb::parallel_for(ID(0), num_nodes, [&](const size_t pos) {
        num_degree_per_vertex[pos] += c[pos];
      });
    }
    Counter num_incident_weight_per_vertex(edge_weight ? num_nodes : 0, 0);
    if (edge_weight) {
      for (Counter& c : local_incident_weight_per_vertex) {
        tbb::parallel_for(ID(0), num_nodes, [&](const size_t pos) {
          num_incident_weight_per_vertex[pos] += c[pos];
        });
      }
    }
    utils::Timer::instance().stop_timer("compute_ds_sizes");

    // Compute prefix sum over the degrees. The prefix sum is used than
    // as start position for each node in the edge array.
    utils::Timer::instance().start_timer("compute_prefix_sums", "Compute Prefix Sums", true);
    parallel::TBBPrefixSum<size_t> degree_prefix_sum(num_degree_per_vertex);
    tbb::parallel_scan(tbb::blocked_range<size_t>( 0UL, UI64(num_nodes)), degree_prefix_sum);
    utils::Timer::instance().stop_timer("compute_prefix_sums");

    utils::Timer::instance().start_timer("setup_hypergraph", "Setup hypergraph", true);
    ASSERT(degree_prefix_sum.total_sum() == 2 * num_edges);

    AtomicCounter incident_edges_position(num_nodes,
                                         parallel::IntegralAtomicWrapper<size_t>(0));

    auto setup_edges = [&] {
      tbb::parallel_for(ID(0), num_edges, [&](const size_t pos) {
        const HypernodeID pin0 = edge_vector[pos].first;
        const HyperedgeID incident_edges_pos0 = degree_prefix_sum[pin0] + incident_edges_position[pin0]++;
        ASSERT(incident_edges_pos0 < graph._edges.size());
        StaticGraph::Edge& edge0 = graph._edges[incident_edges_pos0];
        const HypernodeID pin1 = edge_vector[pos].second;
        const HyperedgeID incident_edges_pos1 = degree_prefix_sum[pin1] + incident_edges_position[pin1]++;
        ASSERT(incident_edges_pos1 < graph._edges.size());
        StaticGraph::Edge& edge1 = graph._edges[incident_edges_pos1];

        edge0.setTarget(pin1);
        edge0.setSource(pin0);
        edge1.setTarget(pin0);
        edge1.setSource(pin1);

        graph._unique_edge_ids[incident_edges_pos0] = pos;
        graph._unique_edge_ids[incident_edges_pos1] = pos;

        if (edge_weight) {
          edge0.setWeight(edge_weight[pos]);
          edge1.setWeight(edge_weight[pos]);
        }
      });
    };

    auto setup_nodes = [&] {
      tbb::parallel_for(ID(0), num_nodes, [&](const size_t pos) {
        StaticGraph::Node& node = graph._nodes[pos];
        node.enable();
        node.setFirstEntry(degree_prefix_sum[pos]);
        node.setIncidentWeight(edge_weight ? num_incident_weight_per_vertex[pos] :
                               degree_prefix_sum[pos + 1] - degree_prefix_sum[pos]);
        if ( node_weight ) {
          node.setWeight(node_weight[pos]);
        }
      });
    };

    auto init_communities = [&] {
      graph._community_ids.resize(num_nodes, 0);
    };

    tbb::parallel_invoke(setup_edges, setup_nodes, init_communities);

    // Add Sentinel
    graph._nodes.back() = StaticGraph::Node(graph._edges.size());
    if (stable_construction_of_incident_edges) {
      sort_incident_edges(graph);
    }
    graph.computeAndSetTotalNodeWeight(parallel_tag_t());
    utils::Timer::instance().stop_timer("setup_hypergraph");
    return graph;
  }
}