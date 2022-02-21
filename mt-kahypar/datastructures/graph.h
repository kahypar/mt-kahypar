/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include <functional>
#include <cmath>
#include <boost/range/irange.hpp>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/range.h"


namespace mt_kahypar {
namespace ds {

/*!
 * CSR Graph Data Structure
 */
class Graph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  struct TmpGraphBuffer {
    explicit TmpGraphBuffer(const size_t num_nodes,
                            const size_t num_arcs) :
      tmp_indices("Preprocessing", "tmp_indices", num_nodes + 1),
      tmp_pos("Preprocessing", "tmp_pos", num_nodes),
      tmp_node_volumes("Preprocessing", "tmp_node_volumes", num_nodes),
      tmp_arcs("Preprocessing", "tmp_arcs", num_arcs),
      valid_arcs("Preprocessing", "valid_arcs", num_arcs) { }

    ds::Array<parallel::IntegralAtomicWrapper<size_t>> tmp_indices;
    ds::Array<parallel::IntegralAtomicWrapper<size_t>> tmp_pos;
    ds::Array<parallel::AtomicWrapper<ArcWeight>> tmp_node_volumes;
    ds::Array<Arc> tmp_arcs;
    ds::Array<size_t> valid_arcs;
  };

 public:
  using AdjacenceIterator = typename ds::Array<Arc>::const_iterator;

 public:
  Graph(Hypergraph& hypergraph, const LouvainEdgeWeight edge_weight_type);
  Graph(Graph&& other);
  Graph& operator= (Graph&& other);
  ~Graph();

  // ! Number of nodes in the graph
  size_t numNodes() const {
    return _num_nodes;
  }

  // ! Number of arcs in the graph
  size_t numArcs() const {
    return _num_arcs;
  }

  // ! Iterator over all nodes of the graph
  auto nodes() const {
    return boost::irange<NodeID>(0, static_cast<NodeID>(numNodes()));
  }

  // ! Iterator over all adjacent vertices of u
  // ! If 'n' is set, then only an iterator over the first n elements is returned
  IteratorRange<AdjacenceIterator> arcsOf(const NodeID u,
                                          const size_t n = std::numeric_limits<size_t>::max()) const {
    ASSERT(u < _num_nodes);
    const size_t start = _indices[u];
    size_t end = _indices[u + 1];
    if ( n < ( end - start ) ) {
      end = start + n;
    }
    return IteratorRange<AdjacenceIterator>(
      _arcs.cbegin() + start, _arcs.cbegin() + end);
  }

  // ! Degree of vertex u
  size_t degree(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _indices[u + 1] - _indices[u];
  }

  // ! Maximum degree of a vertex
  size_t max_degree() const {
    return _max_degree;
  }

  // ! Total Volume of the graph
  ArcWeight totalVolume() const {
    return _total_volume;
  }

  // ! Node volume of vertex u
  ArcWeight nodeVolume(const NodeID u) const {
    ASSERT(u < _num_nodes);
    return _node_volumes[u];
  }

  // ! Projects the clustering of the (likely bipartite star-expansion) graph to the hypergraph
  void restrictClusteringToHypernodes(const Hypergraph& hg, ds::Clustering& C) const {
    C.resize(hg.initialNumNodes());
  }

  bool canBeUsed(const bool verbose = true) const;

  /*!
   * Contracts the graph based on the community structure passed as argument.
   * In the first step the community ids are compactified (via parallel prefix sum)
   * which also determines the node ids in the coarse graph. Afterwards, we create
   * a temporary graph which contains all arcs that will not form a selfloop in the
   * coarse graph. Finally, the weights of each multiedge in that temporary graph
   * are aggregated and the result is written to the final contracted graph.
   */
  Graph contract(Clustering& communities, bool low_memory);

  Graph contract_low_memory(Clustering& communities);

  void allocateContractionBuffers() {
    _tmp_graph_buffer = new TmpGraphBuffer(_num_nodes, _num_arcs);
  }

 private:
  Graph();

  /*!
   * Constructs a graph from a given hypergraph.
   */
  template<typename F>
  void construct(const Hypergraph& hypergraph, const F& edge_weight_func);

  template<typename F>
  void constructBipartiteGraph(const Hypergraph& hypergraph, F& edge_weight_func);

  template<typename F>
  void constructGraph(const Hypergraph& hypergraph,
                      const F& edge_weight_func);

  ArcWeight computeNodeVolume(const NodeID u) {
    ASSERT(u < _num_nodes);
    ArcWeight x = 0.0;
    for (const Arc& arc : arcsOf(u)) {
      x += arc.weight;
    }
    _node_volumes[u] = x;
    return x;
  }

  // ! Number of nodes
  size_t _num_nodes;
  // ! Number of arcs
  size_t _num_arcs;
  // ! Total volume of the graph (= sum of arc weights)
  ArcWeight _total_volume;
  // ! Maximum degree of a node
  size_t _max_degree;

  // ! Index Vector
  ds::Array<size_t> _indices;
  // ! Arcs
  ds::Array<Arc> _arcs;
  // ! Node Volumes (= sum of arc weights for each node)
  ds::Array<ArcWeight> _node_volumes;
  // ! Data that is reused throughout the louvain method
  // ! to construct and contract a graph and to prevent expensive allocations
  TmpGraphBuffer* _tmp_graph_buffer;
};

}  // namespace ds

// expose
using Graph = ds::Graph;

}  // namespace mt_kahypar
