/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <functional>

#include <boost/range/irange.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
namespace ds {
class Graph {

 public:
  using ArcWeight = double;

  struct Arc {
    NodeID head;
    ArcWeight weight;
    Arc(NodeID head, ArcWeight weight) :
      head(head),
      weight(weight) { }
  };

  using AdjList = parallel::scalable_vector<Arc>;

  static constexpr size_t coarseGrainsize = 20000;

  explicit Graph(const size_t num_nodes) :
    _num_arcs(0),
    _total_volume(0.0),
    _adj_list(),
    _edge_weights() {
    tbb::parallel_invoke([&] {
      _adj_list.resize(num_nodes);
    }, [&] {
      _edge_weights.assign(num_nodes, 0);
    });
  }

  ~Graph() {
    tbb::parallel_invoke([&] {
      parallel::parallel_free(_adj_list);
    }, [&] {
      parallel::free(_edge_weights);
    });
  }

  size_t numNodes() const {
    return _adj_list.size();
  }

  void setNumNodes(const size_t n) {
    _adj_list.resize(n);
    _edge_weights.resize(n);
  }

  size_t numArcs() const {
    return _num_arcs;
  }

  void setNumArcs(const size_t numArcs) {
    _num_arcs = numArcs;
  }

  auto nodes() const {
    return boost::irange<NodeID>(0, static_cast<NodeID>(numNodes()));
  }

  void addDirectedEdge(const NodeID u, const NodeID v, const ArcWeight weight) {
    addArc(u, v, weight);
  }

  void addUndirectedEdge(const NodeID u, const NodeID v, const ArcWeight weight) {
    addArc(u, v, weight);
    addArc(v, u, weight);
  }

  void addSelfLoop(const NodeID u, const ArcWeight weight) {
    _edge_weights[u] += weight;
  }

  const AdjList& arcsOf(const NodeID u) const {
    return _adj_list[u];
  }

  AdjList& arcsOf(const NodeID u) {
    return _adj_list[u];
  }

  size_t degree(const NodeID u) const {
    return _adj_list[u].size();
  }

  ArcWeight nodeVolume(const NodeID u) const {
    return _edge_weights[u];
  }

  void setNodeVolume(const NodeID u, const ArcWeight vol) {
    _edge_weights[u] = vol;
  }

  ArcWeight totalVolume() const {
    return _total_volume;
  }

  void setTotalVolume(const ArcWeight totalVolume) {
    _total_volume = totalVolume;
  }

  auto nodesParallelCoarseChunking() const {
    return tbb::blocked_range<NodeID>(0, static_cast<NodeID>(numNodes()), coarseGrainsize);
  }

  void finalize() {
    _total_volume = 0.0;
    _num_arcs = 0;
    tbb::enumerable_thread_specific<ArcWeight> local_total_volume(0.0);
    tbb::enumerable_thread_specific<size_t> local_num_arcs(0);
    tbb::parallel_for(0U, static_cast<NodeID>(numNodes()), [&](const NodeID u) {
      local_num_arcs.local() += degree(u);
      local_total_volume.local() += computeNodeVolume(u);
    });
    _total_volume = local_total_volume.combine(std::plus<ArcWeight>());
    _num_arcs = local_num_arcs.combine(std::plus<size_t>());
  }

 private:
  void addArc(const NodeID u, const NodeID v, const ArcWeight weight) {
    assert(u != v);
    _adj_list[u].emplace_back(v, weight);
  }

  ArcWeight computeNodeVolume(const NodeID u) {
    for (Arc& arc : arcsOf(u))
      _edge_weights[u] += arc.weight;
    return _edge_weights[u];
  }

  size_t _num_arcs;
  ArcWeight _total_volume;
  parallel::scalable_vector<AdjList> _adj_list;
  parallel::scalable_vector<ArcWeight> _edge_weights;
};

class AdjListStarExpansion {
 private:
  using ArcWeight = Graph::ArcWeight;

  AdjListStarExpansion() { }

 public:
  template <typename HyperGraph>
  static Graph constructGraph(const HyperGraph& hg,
                              const Context& context) {
    utils::Timer::instance().start_timer("construct_graph", "Construct Graph");
    Graph graph(hg.initialNumNodes() + hg.initialNumEdges());

    switch (context.preprocessing.community_detection.edge_weight_function) {
      case LouvainEdgeWeight::degree:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID hn) -> Graph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he)) * static_cast<ArcWeight>(hg.nodeDegree(hn)) / static_cast<ArcWeight>(hg.edgeSize(he));
          });
        break;
      case LouvainEdgeWeight::non_uniform:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID) -> Graph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he)) / static_cast<ArcWeight>(hg.edgeSize(he));
          });
        break;
      case LouvainEdgeWeight::uniform:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID) -> Graph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he));
          });
        break;
      case LouvainEdgeWeight::hybrid:
        ERROR("Only uniform/non-uniform/degree edge weight is allowed at graph construction.");
      default:
        ERROR("Unknown edge weight for bipartite graph.");
    }
    utils::Timer::instance().stop_timer("construct_graph");

    return graph;
  }

  template <typename HyperGraph>
  static void restrictClusteringToHypernodes(const HyperGraph& hg, Clustering& C) {
    C.resize(hg.initialNumNodes());
    C.shrink_to_fit();
    C.compactify();
  }

 private:
  // TaggedInteger would be great here
  template <typename HyperGraph>
  static NodeID mapHypernode(const HyperGraph& hg, const HypernodeID u) {
    return hg.originalNodeID(u);
  }

  template <typename HyperGraph>
  static NodeID mapHyperedge(const HyperGraph& hg, const HyperedgeID e) {
    return hg.initialNumNodes() + hg.originalEdgeID(e);
  }

  template <typename HyperGraph>
  static HypernodeID mapNodeToHypernode(const HyperGraph& hg, const NodeID u) {
    return hg.globalNodeID(u);
  }

  template <typename HyperGraph>
  static HyperedgeID mapNodeToHyperedge(const HyperGraph& hg, const NodeID u) {
    return hg.globalEdgeID(u - hg.initialNumNodes());
  }

  template <class HyperGraph, class EdgeWeightFunction>
  static void fill(Graph& graph,
                   const HyperGraph& hg,
                   EdgeWeightFunction ewf) {
    tbb::parallel_invoke(
      [&]() {
          // WARNING! This function does not skip deactivated nodes because KaHyPar exposes pairs of iterators, not a range type that implements .empty() as required by parallel_for
          tbb::parallel_for(NodeID(0), NodeID(hg.initialNumNodes()), [&](const HypernodeID original_id) {
            const HypernodeID hn = hg.globalNodeID(original_id);
            const NodeID graph_hn = mapHypernode(hg, hn);
            for (const HyperedgeID& he : hg.incidentEdges(hn))
              graph.addDirectedEdge(graph_hn, mapHyperedge(hg, he), ewf(hg, he, hn));
          });
        },

      [&]() {
          tbb::parallel_for(HyperedgeID(0), HyperedgeID(hg.initialNumEdges()), [&](const HyperedgeID original_id) {
            const HyperedgeID he = hg.globalEdgeID(original_id);
            const NodeID graph_he = mapHyperedge(hg, he);
            for (const HypernodeID& hn : hg.pins(he))
              graph.addDirectedEdge(graph_he, mapHypernode(hg, hn), ewf(hg, he, hn));
          });
        }
      );

    graph.finalize();
  }
};

}  // namespace ds
}  // namespace mt_kahypar
