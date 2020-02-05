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

#include <type_traits>

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// Forward
template <typename Hypergraph,
          typename Factory,
          typename HardwareTopology,
          typename TBBNumaArena>
class NumaHypergraphFactory;

template <typename Hypergraph = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class NumaHypergraph {

  static_assert(!Hypergraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = true;
  static constexpr bool is_partitioned = false;

  explicit NumaHypergraph() :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_degree(0),
    _total_weight(0),
    _hypergraphs(),
    _node_mapping(),
    _edge_mapping(),
    _community_node_mapping() { }

  NumaHypergraph(const NumaHypergraph&) = delete;
  NumaHypergraph & operator= (const NumaHypergraph &) = delete;

  NumaHypergraph(NumaHypergraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_pins(other._num_pins),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _hypergraphs(std::move(other._hypergraphs)),
    _node_mapping(std::move(other._node_mapping)),
    _edge_mapping(std::move(other._edge_mapping)),
    _community_node_mapping(std::move(other._community_node_mapping)) { }

  NumaHypergraph & operator= (NumaHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _hypergraphs = std::move(other._hypergraphs);
    _node_mapping = std::move(other._node_mapping);
    _edge_mapping = std::move(other._edge_mapping);
    _community_node_mapping = std::move(other._community_node_mapping);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Number of NUMA hypergraphs
  size_t numNumaHypergraphs() const {
    return _hypergraphs.size();
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Initial number of hypernodes on numa node
  HypernodeID initialNumNodes(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].initialNumNodes();
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    HypernodeID num_removed_hypernodes = 0;
    for ( const Hypergraph& hypergraph : _hypergraphs ) {
      num_removed_hypernodes += hypergraph.numRemovedHypernodes();
    }
    return num_removed_hypernodes;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  // ! Initial number of hyperedges on numa node
  HyperedgeID initialNumEdges(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].initialNumEdges();
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Initial number of pins on numa node
  HypernodeID initialNumPins(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].initialNumPins();
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _total_degree;
  }

  // ! Initial sum of the degree of all vertices on numa node
  HypernodeID initialTotalVertexDegree(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].initialTotalVertexDegree();
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Recomputes the total weight of the hypergraph (in parallel)
  void updateTotalWeight(const TaskGroupID task_group_id) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].updateTotalWeight();
        });

    _total_weight = 0;
    for ( Hypergraph& hypergraph : _hypergraphs ) {
      _total_weight += hypergraph.totalWeight();
    }
  }

  // ! Recomputes the total weight of the hypergraph (sequential)
  void updateTotalWeight() {
    _total_weight = 0;
    for ( Hypergraph& hypergraph : _hypergraphs ) {
      hypergraph.updateTotalWeight();
      _total_weight += hypergraph.totalWeight();
    }
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID task_group_id, const F& f) {
    static_cast<const NumaHypergraph&>(*this).doParallelForAllNodes(task_group_id, f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID task_group_id, const F& f) const {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        _hypergraphs[node].doParallelForAllNodes(task_group_id, f);
      });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) {
    static_cast<const NumaHypergraph&>(*this).doParallelForAllEdges(task_group_id, f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) const {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        _hypergraphs[node].doParallelForAllEdges(task_group_id, f);
      });
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  ConcatenatedRange<IteratorRange<HypernodeIterator>> nodes() const {
    ASSERT(!_hypergraphs.empty());
    ConcatenatedRange<IteratorRange<HypernodeIterator>> iterator;
    for (const Hypergraph& hg : _hypergraphs) {
      iterator.concat(hg.nodes());
    }
    return iterator;
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HypernodeIterator> nodes(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  ConcatenatedRange<IteratorRange<HyperedgeIterator>> edges() const {
    ASSERT(!_hypergraphs.empty());
    ConcatenatedRange<IteratorRange<HyperedgeIterator>> iterator;
    for (const Hypergraph& hg : _hypergraphs) {
      iterator.concat(hg.edges());
    }
    return iterator;
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HyperedgeIterator> edges(const int node) const {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node].edges();
  }

  /*!
   * Illustration for different incidentEdges iterators:
   *
   * Structure of incident nets for a vertex u during community coarsening:
   *
   * | <-- single-pin community hyperedges --> | <-- multi-pin community hyperedges --> | <-- invalid hyperedges --> |
   *
   *                                             <--   multiPinIncidentEdges(u,c)   -->
   *
   *  <--------------------------- activeIncidentEdges(u,c) -------------------------->
   *
   *  <---------------------------------------------- incidentEdges(u) --------------------------------------------->
   */

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).incidentEdges(u);
  }

  // ! Returns a range to loop over all VALID hyperedges of hypernode u.
  IteratorRange<IncidentNetsIterator> multiPinIncidentEdges(const HypernodeID u, const PartitionID community_id) const {
    return hypergraph_of_vertex(u).multiPinIncidentEdges(u, community_id);
  }

  // TODO function name should reflect its purpose
  // ! Returns a range to loop over the set of all VALID incident hyperedges of hypernode u that are not single-pin community hyperedges.
  IteratorRange<IncidentNetsIterator> activeIncidentEdges(const HypernodeID u, const PartitionID community_id) const {
    return hypergraph_of_vertex(u).activeIncidentEdges(u, community_id);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    return hypergraph_of_edge(e).pins(e);
  }

  // ! Returns a range to loop over the pins of hyperedge e that belong to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).pins(e, community_id);
  }

  // ! Returns a range to loop over the set of communities contained in hyperedge e.
  // IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! Can be used to map the global vertex ids to a consecutive range
  // ! of nodes between [0,|V|).
  HypernodeID originalNodeID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).originalNodeID(u);
  }

  // ! Reverse operation of originalNodeID(u)
  HypernodeID globalNodeID(const HypernodeID u) const {
    ASSERT(u < _node_mapping.size());
    return _node_mapping[u];
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    hypergraph_of_vertex(u).setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeDegree(u);
  }

  // ! Returns, if the corresponding vertex is high degree vertex
  bool isHighDegreeVertex(const HypernodeID u) const {
    return hypergraph_of_vertex(u).isHighDegreeVertex(u);
  }

  // ! Marks all vertices with a degree greater the threshold
  // ! as high degree vertex
  void markAllHighDegreeVertices(const TaskGroupID task_group_id,
                                 const HypernodeID high_degree_threshold) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].markAllHighDegreeVertices(task_group_id, high_degree_threshold);
        });
  }

  // ! Number of invalid incident nets
  HyperedgeID numInvalidIncidentNets(const HypernodeID u) const {
    return hypergraph_of_vertex(u).numInvalidIncidentNets(u);
  }

  // ! Contraction index of the vertex in the contraction hierarchy
  HypernodeID contractionIndex(const HypernodeID u) const {
    return hypergraph_of_vertex(u).contractionIndex(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return hypergraph_of_vertex(u).nodeIsEnabled(u);
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).disableHypernode(u);
  }

  // ! Removes a hypernode (must be enabled before)
  void removeHypernode(const HypernodeID u) {
    hypergraph_of_vertex(u).removeHypernode(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Returns for a edge of the hypergraph its original edge id
  // ! Can be used to map the global edge ids to a consecutive range
  // ! of edges between [0,|E|).
  HypernodeID originalEdgeID(const HyperedgeID e) const {
    return hypergraph_of_edge(e).originalEdgeID(e);
  }

  // ! Reverse operation of originalEdgeID(e)
  HypernodeID globalEdgeID(const HyperedgeID e) const {
    ASSERT(e < _edge_mapping.size());
    return _edge_mapping[e];
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeWeight(e);
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    hypergraph_of_edge(e).setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeSize(e);
  }

  // ! Hash value defined over the pins of a hyperedge
  size_t edgeHash(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeHash(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return hypergraph_of_edge(e).edgeIsEnabled(e);
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hypergraph_of_edge(e).disableHyperedge(e);
  }

  // ####################### Community Hyperedge Information #######################

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeWeight(e, community_id);
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    hypergraph_of_edge(e).setEdgeWeight(e, community_id, weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeSize(e, community_id);
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    return hypergraph_of_edge(e).edgeHash(e, community_id);
  }

  // ####################### Community Information #######################

  // ! Number of communities
  PartitionID numCommunities() const {
    return _hypergraphs[0].numCommunities();
  }

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).communityID(u);
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID u, const PartitionID community_id) {
    hypergraph_of_vertex(u).setCommunityID(u, community_id);
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  // HypernodeID communityNodeId(const HypernodeID u) const

  // ! Number of hypernodes in community
  HypernodeID numCommunityHypernodes(const PartitionID community) const {
    HypernodeID num_community_hypernodes = 0;
    for ( const Hypergraph& hypergraph : _hypergraphs ) {
      num_community_hypernodes += hypergraph.numCommunityHypernodes(community);
    }
    return num_community_hypernodes;
  }

  // ! Number of pins in community
  HypernodeID numCommunityPins(const PartitionID community) const {
    HypernodeID num_community_pins = 0;
    for ( const Hypergraph& hypergraph : _hypergraphs ) {
      num_community_pins += hypergraph.numCommunityPins(community);
    }
    return num_community_pins;
  }

  // ! Total degree of community
  HyperedgeID communityDegree(const PartitionID community) const {
    HypernodeID total_community_degree = 0;
    for ( const Hypergraph& hypergraph : _hypergraphs ) {
      total_community_degree += hypergraph.communityDegree(community);
    }
    return total_community_degree;
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    return hypergraph_of_vertex(e).numCommunitiesInHyperedge(e);
  }

  // ! Numa node to which community is assigned to
  PartitionID communityNumaNode(const PartitionID community_id) const {
    ASSERT(static_cast<size_t>(community_id) < _community_node_mapping.size());
    return _community_node_mapping[community_id];
  }

  // ! Sets the community to numa node mapping
  void setCommunityNodeMapping(std::vector<PartitionID>&& community_node_mapping) {
    ASSERT(community_node_mapping.size() == static_cast<size_t>(numCommunities()));
    _community_node_mapping = std::move(community_node_mapping);
  }

  // ####################### Initialization / Reset Functions #######################

  /*!
   * Initializes community-related information after all vertices are assigned to a community.
   * This includes:
   *  1.) Number of Communities
   *  2.) Number of Vertices per Community
   *  3.) Number of Pins per Community
   *  4.) For each hypernode v of community C, we compute a unique id within
   *      that community in the range [0, |C|)
   */
  void initializeCommunities(const TaskGroupID task_group_id) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].initializeCommunities(task_group_id, _hypergraphs);
        });
  }

  /*!
  * Initializes community hyperedges.
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to a range of consecutive pins with
  *       same community in that hyperedge
  */
  void initializeCommunityHyperedges(const TaskGroupID task_group_id) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].initializeCommunityHyperedges(task_group_id, _hypergraphs);
        });
  }

 private:
  template <typename HyperGraph,
            typename Factory,
            typename HwTopology,
            typename TBB>
  friend class NumaHypergraphFactory;

  // ####################### Helper Functions #######################

  const Hypergraph& hypergraph_of_vertex(const HypernodeID u) const {
    int node = common::get_numa_node_of_vertex(u);
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node];
  }

  Hypergraph& hypergraph_of_vertex(const HypernodeID u) {
    return const_cast<Hypergraph&>(static_cast<const NumaHypergraph&>(*this).hypergraph_of_vertex(u));
  }

  const Hypergraph& hypergraph_of_edge(const HyperedgeID e) const {
    int node = common::get_numa_node_of_edge(e);
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node];
  }

  Hypergraph& hypergraph_of_edge(const HyperedgeID e) {
    return const_cast<Hypergraph&>(static_cast<const NumaHypergraph&>(*this).hypergraph_of_edge(e));
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! NUMA Hypergraphs
  parallel::scalable_vector<Hypergraph> _hypergraphs;
  // ! Mapping from original node id to its hypergraph node id
  parallel::scalable_vector<HypernodeID> _node_mapping;
  // ! Mapping from original edge id to its hypergraph edge id
  parallel::scalable_vector<HyperedgeID> _edge_mapping;
  // ! Mapping from community id to its NUMA node which they are assigned to
  parallel::scalable_vector<PartitionID> _community_node_mapping;
};

} // namespace ds
} // namespace mt_kahypar