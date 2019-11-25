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
#include <tbb/parallel_reduce.h>

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
namespace ds {
class AdjListGraph {
 public:
  using ArcWeight = double;

  struct Arc {
    NodeID head;
    ArcWeight weight;
    Arc(NodeID head, ArcWeight weight) :
      head(head),
      weight(weight) { }
  };

  using AdjList = std::vector<Arc>;

  static constexpr size_t coarseGrainsize = 20000;

  explicit AdjListGraph(const size_t numNodes) :
    _numArcs(0),
    _totalVolume(0.0),
    _adj(numNodes),
    _volume(numNodes, 0) {  /* , selfLoopWeight(numNodes, 0) */
  }

  void addArc(const NodeID u, const NodeID v, const ArcWeight weight) {
    assert(u != v);
    _adj[u].emplace_back(v, weight);
  }

  void addHalfEdge(const NodeID u, const NodeID v, const ArcWeight weight) {
    addArc(u, v, weight);
  }

  void addEdge(const NodeID u, const NodeID v, const ArcWeight weight) {
    addArc(u, v, weight);
    addArc(v, u, weight);
  }

  void addSelfLoop(const NodeID u, const ArcWeight weight) {
    _volume[u] += weight;        // since the contraction will also
  }

  // ArcWeight getSelfLoopWeight(const NodeID u) const { return selfLoopWeight[u]; }

  const AdjList& arcsOf(const NodeID u) const {
    return _adj[u];
  }

  AdjList& arcsOf(const NodeID u) {
    return _adj[u];
  }

  size_t degree(const NodeID u) const {
    return _adj[u].size();
  }

  ArcWeight nodeVolume(const NodeID u) const {
    return _volume[u];
  }

  ArcWeight computeNodeVolume(const NodeID u) {
    for (Arc& arc : arcsOf(u))
      _volume[u] += arc.weight;
    return _volume[u];
  }

  void setNodeVolume(const NodeID u, const ArcWeight vol) {
    _volume[u] = vol;
  }

  ArcWeight totalVolume() const {
    return _totalVolume;
  }

  void setTotalVolume(const ArcWeight totalVolume) {
    _totalVolume = totalVolume;
  }

  size_t numNodes() const {
    return _adj.size();
  }

  void setNumNodes(const size_t n) {
    _adj.resize(n);
    _volume.resize(n);
  }

  auto nodes() const {
    return boost::irange<NodeID>(0, static_cast<NodeID>(numNodes()));
  }

  size_t numArcs() const {
    return _numArcs;
  }

  void setNumArcs(const size_t numArcs) {
    _numArcs = numArcs;
  }

  auto nodesParallelCoarseChunking() const {
    return tbb::blocked_range<NodeID>(0, static_cast<NodeID>(numNodes()), coarseGrainsize);
  }

  void finalize(bool parallel = false) {
    if (parallel) {
      // Do in one pass
      LOG << "Finalizing AdjListGraph in parallel is currently deactivated. Abort";
      std::abort();

#if 0
      using PassType = std::pair<ArcWeight, size_t>;
      std::tie(_totalVolume, _numArcs) = tbb::parallel_reduce(
        /* indices */
        nodes(),
        // nodesParallelCoarseChunking(),
        /* initial */
        std::make_pair(ArcWeight(0.0), size_t(0)),
        /* accumulate */
        [&](const tbb::blocked_range<NodeID>& r, PassType current) -> PassType {
            for (NodeID u = r.begin(), last = r.end(); u < last; ++u) {             // why does that not boil down to a range-based for-loop???
              current.first += computeNodeVolume(u);
              current.second += degree(u);
            }
            return current;
          },
        /* join */
        [](const PassType& l, const PassType& r) -> PassType {
            return std::make_pair(l.first + r.first, l.second + r.second);
          }
        );
#endif
    } else {
      _totalVolume = 0.0;
      _numArcs = 0;
      for (const NodeID u : nodes()) {
        _numArcs += degree(u);
        _totalVolume += computeNodeVolume(u);
      }
    }
  }

 private:
  size_t _numArcs;
  ArcWeight _totalVolume;
  std::vector<AdjList> _adj;
  std::vector<ArcWeight> _volume;
  // std::vector<ArcWeight> selfLoopWeight;
};

// Ben Style template everything and do ID mapping on-the-fly in case copying stuff becomes a bottleneck
class AdjListStarExpansion {
 private:
  using ArcWeight = AdjListGraph::ArcWeight;

  AdjListStarExpansion() { }

 public:
  template< typename HyperGraph >
  static AdjListGraph contructGraph(const HyperGraph& hg, const Context& context) {
    AdjListGraph graph(hg.initialNumNodes() + hg.initialNumEdges());
    bool isGraph = true;
    for (const HyperedgeID he : hg.edges()) {
      if (hg.edgeSize(he) > 2) {
        isGraph = false;
        break;
      }
    }

    if (isGraph) {
      graph.setNumNodes(hg.initialNumNodes());
      for (HypernodeID u : hg.nodes()) {
        for (HyperedgeID e : hg.incidentEdges(u)) {
          for (HypernodeID p : hg.pins(e)) {
            if (p != u)
              graph.addHalfEdge(mapHypernode(hg, u), mapHypernode(hg, p), hg.edgeWeight(e));
          }
        }
      }
      graph.finalize();
      return graph;
    }

    // This is literally disgusting
    switch (context.preprocessing.community_detection.edge_weight_function) {
      case LouvainEdgeWeight::degree:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID hn) -> AdjListGraph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he)) * static_cast<ArcWeight>(hg.nodeDegree(hn)) / static_cast<ArcWeight>(hg.edgeSize(he));
          });
        break;
      case LouvainEdgeWeight::non_uniform:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID) -> AdjListGraph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he)) / static_cast<ArcWeight>(hg.edgeSize(he));
          });
        break;
      case LouvainEdgeWeight::uniform:
        fill(graph, hg, [&](const HyperGraph& hg, const HyperedgeID he, const HypernodeID) -> AdjListGraph::ArcWeight {
            return static_cast<ArcWeight>(hg.edgeWeight(he));
          });
        break;
      case LouvainEdgeWeight::hybrid:
        LOG << "Only uniform/non-uniform/degree edge weight is allowed at graph construction.";
        std::exit(-1);
      default:
        LOG << "Unknown edge weight for bipartite graph.";
        std::exit(-1);
    }

    return graph;
  }

  template< typename HyperGraph >
  static void restrictClusteringToHypernodes(const HyperGraph& hg, Clustering& C) {
    C.resize(hg.initialNumNodes());
    C.shrink_to_fit();
  }

 private:

  // TaggedInteger would be great here

  template< typename HyperGraph >
  static NodeID mapHyperedge(const HyperGraph& hg, const HyperedgeID e) {
    return hg.initialNumNodes() + hg.originalEdgeID(e);
  }

  template< typename HyperGraph >
  static NodeID mapHypernode(const HyperGraph& hg, const HypernodeID u) {
    return hg.originalNodeID(u);
  }

  template< typename HyperGraph >
  static HypernodeID mapNodeToHypernode(const HyperGraph& hg, const NodeID u) {
    return hg.globalNodeID(u);
  }

  template< typename HyperGraph >
  static HyperedgeID mapNodeToHyperedge(const HyperGraph& hg, const NodeID u) {
    return hg.globalEdgeID(u - hg.initialNumNodes());
  }

  template <class HyperGraph, class EdgeWeightFunction>
  static void fill(AdjListGraph& graph,
                   const HyperGraph& hg,
                   EdgeWeightFunction ewf,
                   bool parallel = false) {
    if (parallel) {
      tbb::parallel_invoke(
        [&]() {
            // WARNING! This function does not skip deactivated nodes because KaHyPar exposes pairs of iterators, not a range type that implements .empty() as required by parallel_for
            tbb::parallel_for(NodeID(0), NodeID(hg.initialNumNodes()), [&](const HypernodeID hn) {
              const NodeID graph_hn = mapHypernode(hg, hn);
              for (HyperedgeID he : hg.incidentEdges(hn))
                graph.addHalfEdge(graph_hn, mapHyperedge(hg, he), ewf(hg, he, hn));
            });
          },

        [&]() {
            tbb::parallel_for(HyperedgeID(0), HyperedgeID(hg.initialNumEdges()), [&](const HyperedgeID he) {
              const NodeID graph_he = mapHyperedge(hg, he);
              for (HypernodeID hn : hg.pins(he))
                graph.addHalfEdge(graph_he, mapHypernode(hg, hn), ewf(hg, he, hn));
            });
          }
        );
    } else {
      for (HypernodeID hn : hg.nodes())
        for (HyperedgeID he : hg.incidentEdges(hn))
          graph.addEdge(mapHypernode(hg, hn), mapHyperedge(hg, he), ewf(hg, he, hn));
    }

    graph.finalize(parallel);
  }
};

class CSRGraph {
  // TODO. refactor commons of AdjListGraph into GraphBase
  // since contraction is super fast, we might even screw the parallelization there if the benefit to local moving justifies it
};
}  // namespace ds
}  // namespace mt_kahypar
