/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
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

#pragma once

#include <boost/range/irange.hpp>

#include "tbb/parallel_for.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

class SeparatedNodes {

  static constexpr bool enable_heavy_assert = false;

  class Node {
   public:
    Node() :
      original_node(kInvalidHypernode),
      begin(0),
      weight(1) { }

    explicit Node(const HypernodeID original_node, const HyperedgeID begin, const HypernodeWeight weight) :
      original_node(original_node),
      begin(begin),
      weight(weight) { }

    // ! ID of the original node in the graph
    HypernodeID original_node;
    // ! Index of the first element in _inward_edges
    HyperedgeID begin;
    // ! Node weight
    HypernodeWeight weight;
  };

 public:
  class Edge {
   public:
    Edge() :
      target(0),
      weight(1) { }

    explicit Edge(HypernodeID target, HyperedgeWeight weight) :
      target(target),
      weight(weight) { }

    // ! index of target node
    HypernodeID target;
    // ! hyperedge weight
    HyperedgeWeight weight;
  };

  // ! Iterator to iterate over the nodes
  using HypernodeIterator = boost::range_detail::integer_iterator<HyperedgeID>;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = const Edge*;

  explicit SeparatedNodes() :
    _num_nodes(0),
    _num_graph_nodes(0),
    _num_edges(0),
    _total_weight(0),
    _nodes{ Node(kInvalidHypernode, 0, 0) },
    _outward_incident_weight(),
    _graph_nodes_begin(),
    _inward_edges(),
    _outward_edges(),
    _batch_indices_and_weights(),
    _savepoints() {
      _batch_indices_and_weights.assign(2, {0, 0});
    }

  explicit SeparatedNodes(HypernodeID num_graph_nodes) :
    _num_nodes(0),
    _num_graph_nodes(num_graph_nodes),
    _num_edges(0),
    _total_weight(0),
    _nodes{ Node(kInvalidHypernode, 0, 0) },
    _outward_incident_weight(num_graph_nodes),
    _graph_nodes_begin(),
    _inward_edges(),
    _outward_edges(),
    _batch_indices_and_weights(),
    _savepoints() {
      _batch_indices_and_weights.assign(2, {0, 0});
    }

  SeparatedNodes(const SeparatedNodes&) = delete;
  SeparatedNodes & operator= (const SeparatedNodes &) = delete;

  SeparatedNodes(SeparatedNodes&& other) :
    _num_nodes(other._num_nodes),
    _num_graph_nodes(other._num_graph_nodes),
    _num_edges(other._num_edges),
    _total_weight(other._total_weight),
    _nodes(std::move(other._nodes)),
    _outward_incident_weight(std::move(other._outward_incident_weight)),
    _graph_nodes_begin(std::move(other._graph_nodes_begin)),
    _inward_edges(std::move(other._inward_edges)),
    _outward_edges(std::move(other._outward_edges)),
    _batch_indices_and_weights(std::move(other._batch_indices_and_weights)),
    _savepoints(std::move(other._savepoints)) { }

  SeparatedNodes & operator= (SeparatedNodes&& other) {
    _num_nodes = other._num_nodes;
    _num_graph_nodes = other._num_graph_nodes;
    _num_edges = other._num_edges;
    _total_weight = other._total_weight;
    _nodes = std::move(other._nodes);
    _outward_incident_weight = std::move(other._outward_incident_weight);
    _graph_nodes_begin = std::move(other._graph_nodes_begin);
    _inward_edges = std::move(other._inward_edges);
    _outward_edges = std::move(other._outward_edges);
    _batch_indices_and_weights = std::move(other._batch_indices_and_weights);
    _savepoints = std::move(other._savepoints);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Number of hypernodes
  HypernodeID numNodes() const {
    return _num_nodes;
  }

  HypernodeID numGraphNodes() const {
    return _num_graph_nodes;
  }

  // ! Number of edges
  HypernodeID numEdges() const {
    return _num_edges;
  }

  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ####################### Iterators #######################

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return IteratorRange<HypernodeIterator>(
      boost::range_detail::integer_iterator<HypernodeID>(0),
      boost::range_detail::integer_iterator<HypernodeID>(_num_nodes));
  }

  // ! Edges leading from the graph node outward to a separated node
  IteratorRange<IncidenceIterator> outwardEdges(const HypernodeID graph_node) const {
    ASSERT(graph_node < _num_graph_nodes, "Graph node" << graph_node << "does not exist");
    ASSERT(!_graph_nodes_begin.empty(), "Graph nodes not initialized!");
    return IteratorRange<IncidenceIterator>(
      _outward_edges.data() + _graph_nodes_begin[graph_node],
      _outward_edges.data() + _graph_nodes_begin[graph_node + 1]);
  }

  // ! Edges leading from the separated node to a graph node
  IteratorRange<IncidenceIterator> inwardEdges(const HypernodeID separated_node) const {
    return IteratorRange<IncidenceIterator>(
      _inward_edges.data() + node(separated_node).begin,
      _inward_edges.data() + node(separated_node + 1).begin);
  }

    // ####################### Node Information #######################

  // ! Original ID of a vertex
  HypernodeID originalHypernodeID(const HypernodeID u) const {
    return node(u).original_node;
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return node(u).weight;
  }

  HyperedgeID inwardDegree(const HypernodeID separated_node) const {
    return node(separated_node + 1).begin - node(separated_node).begin;
  }

  HyperedgeID outwardDegree(const HypernodeID u) const {
    ASSERT(u < _num_graph_nodes, "Graph node" << u << "does not exist");
    ASSERT(!_graph_nodes_begin.empty(), "Graph nodes not initialized!");
    return _graph_nodes_begin[u + 1] - _graph_nodes_begin[u];
  }

  HyperedgeWeight outwardIncidentWeight(const HypernodeID u) const {
    ASSERT(u < _num_graph_nodes, "Graph node" << u << "does not exist");
    return _outward_incident_weight[u].load();
  }

  // ####################### Batches and Savepoints #######################

  HypernodeID currentBatchIndex() const {
    return _batch_indices_and_weights[_batch_indices_and_weights.size() - 2].first;
  }

  // ! Returns the index of the current batch, which is also the number
  // ! of nodes left after the pop operation.
  HypernodeID popBatch();

  void cleanBatchState();

  void setSavepoint();

  void restoreSavepoint();

  // ####################### Contract / Uncontract #######################

  /*!
   * Adds the given nodes. Each node is specified by its weight and its
   * starting index in the edge vector.
   */
  void addNodes(const vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>>& nodes, const vec<Edge>& edges);

  /*!
   * Contracts a given community structure. Note that the mapping must contain
   * the node id of the node in the coarsened graph (i.e. as returned by StaticGraph::contract).
   *
   * If kInvalidHypernode is used as community id, the node is considered removed.
   */
  void contract(const vec<HypernodeID>& communities, const HypernodeID& num_coarsened_graph_nodes);

  void initializeOutwardEdges();

  // ####################### Extract Block #######################

  // ! extracts one block as new separated nodes structure
  SeparatedNodes extract(PartitionID block, const vec<HypernodeID>& graph_node_mapping,
                         const vec<CAtomic<PartitionID>>& part_ids) const;

  // ####################### Initialization / Reset Functions #######################

  // ! Copy static hypergraph in parallel
  SeparatedNodes copy(parallel_tag_t) const;

  // ! Copy static hypergraph sequential
  SeparatedNodes copy() const;

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

 private:
  static_assert(std::is_trivially_copyable<Node>::value, "Node is not trivially copyable");
  static_assert(std::is_trivially_copyable<Edge>::value, "Hyperedge is not trivially copyable");

  // ####################### Node Information #######################

  // ! Accessor for node-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Node& node(const HypernodeID u) const {
    ASSERT(u <= _num_nodes, "Node" << u << "does not exist");
    return _nodes[u];
  }

  // ! Accessor for node-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Node& node(const HypernodeID u) {
    ASSERT(u <= _num_nodes, "Node" << u << "does not exist");
    return _nodes[u];
  }

  // ! Number of separated nodes
  HypernodeID _num_nodes;
  // ! Number of graph nodes
  HypernodeID _num_graph_nodes;
  // ! Number of edges (note: we have an outward and inward edge for each)
  HyperedgeID _num_edges;
  // ! Total weight of the nodes
  HypernodeWeight _total_weight;

  // ! Nodes
  vec<Node> _nodes;
  vec<parallel::IntegralAtomicWrapper<HyperedgeWeight>> _outward_incident_weight;
  Array<HyperedgeID> _graph_nodes_begin;

  // ! Edges
  vec<Edge> _inward_edges;
  Array<Edge> _outward_edges;

  // ! Batches
  // - allow to restore nodes of a previous state (but not edges)
  vec<std::pair<HyperedgeID, HypernodeWeight>> _batch_indices_and_weights;
  // ! Savepoints
  // - allow to restore edges of a previous state, must be set manually
  // edges, edge_indices, num_graph_nodes
  vec<std::tuple<vec<Edge>, vec<HyperedgeID>, HypernodeID>> _savepoints;
};

} // namespace ds
} // namespace mt_kahypar
