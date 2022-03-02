/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "dynamic_graph_factory.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {
  DynamicGraph DynamicGraphFactory::construct(
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

  DynamicGraph DynamicGraphFactory::construct_from_graph_edges(
          const HypernodeID num_nodes,
          const HyperedgeID num_edges,
          const EdgeVector& edge_vector,
          const HyperedgeWeight* edge_weight,
          const HypernodeWeight* node_weight,
          // TODO(maas): should we support stable edge construction?
          const bool stable_construction_of_incident_edges) {
    ASSERT(edge_vector.size() == num_edges);
    DynamicGraph graph;
    graph._num_edges = 2 * num_edges;

    // TODO: calculate required id range
    tbb::parallel_invoke([&] {
      graph._nodes.resize(num_nodes + 1);
      tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID n) {
        // setup nodes
        DynamicGraph::Hypernode& node = graph._nodes[n];
        node.enable();
        if ( node_weight ) {
          node.setWeight(node_weight[n]);
        }
      });
      // Compute total weight of graph
      graph.updateTotalWeight(parallel_tag_t());
    }, [&] {
      graph._adjacency_array = DynamicAdjacencyArray(num_nodes, edge_vector, edge_weight);
      if (stable_construction_of_incident_edges) {
        graph._adjacency_array.sortIncidentEdges();
      }
    }, [&] {
      graph._acquired_nodes.assign(
        num_nodes, parallel::IntegralAtomicWrapper<bool>(false));
    }, [&] {
      graph._contraction_tree.initialize(num_nodes);
    });
    return graph;
  }
}
