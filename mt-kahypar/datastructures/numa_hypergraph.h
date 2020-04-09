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
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// Forward
template <typename Hypergraph,
          typename Factory,
          typename HardwareTopology,
          typename TBBNumaArena>
class NumaHypergraphFactory;

template <typename HyperGraph = Mandatory,
          typename HwTopology = Mandatory,
          typename TBB = Mandatory>
class NumaHypergraph {

  static_assert(!HyperGraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");
  static_assert(!HyperGraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

 static constexpr bool enable_heavy_assert = false;

  // During contractions we temporary memcpy all incident nets of a collapsed
  // vertex to consecutive range in a temporary incident nets structure.
  // Afterwards, we sort that range and remove duplicates. However, it turned
  // out that this become a major sequential bottleneck in presence of high
  // degree vertices. Therefore, all vertices with temporary degree greater
  // than this threshold are contracted with a special procedure.
  static constexpr HyperedgeID HIGH_DEGREE_CONTRACTION_THRESHOLD = ID(500000);

 public:
  // ! Type Traits
  using Hypergraph = HyperGraph;
  using HardwareTopology = HwTopology;
  using TBBNumaArena = TBB;

  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = true;
  static constexpr bool is_partitioned = false;

  static_assert(sizeof(HypernodeID) == 8, "Hypernode ID must be 8 byte");
  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(sizeof(HyperedgeID) == 8, "Hyperedge ID must be 8 byte");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
  // ! Iterator to iterate over the set of communities contained in a hyperedge
  using CommunityIterator = typename Hypergraph::CommunityIterator;
  // ! Buffer used for contractions
  using TmpContractionBuffer = typename Hypergraph::TmpContractionBuffer;

  explicit NumaHypergraph() :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_removed_hyperedges(0),
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
    _num_removed_hyperedges(other._num_removed_hyperedges),
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
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _hypergraphs = std::move(other._hypergraphs);
    _node_mapping = std::move(other._node_mapping);
    _edge_mapping = std::move(other._edge_mapping);
    _community_node_mapping = std::move(other._community_node_mapping);
    return *this;
  }

  ~NumaHypergraph() {
    freeInternalData();
  }

  // ####################### General Hypergraph Stats #######################

  // ! Number of NUMA hypergraphs
  size_t numNumaHypergraphs() const {
    return _hypergraphs.size();
  }

  Hypergraph& numaHypergraph(const int node) {
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node];
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

  // ! Number of removed hyperedges
  HyperedgeID numRemovedHyperedges() const {
    return _num_removed_hyperedges;
  }

  // ! Set the number of removed hyperedges
  void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
    _num_removed_hyperedges = num_removed_hyperedges;
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

  // ! Returns a range to loop over all active multi-pin hyperedges of hypernode u.
  IteratorRange<IncidentNetsIterator> multiPinIncidentEdges(const HypernodeID u, const PartitionID community_id) const {
    return hypergraph_of_vertex(u).multiPinIncidentEdges(u, community_id);
  }

  // ! Returns a range to loop over the set of all active incident
  // ! hyperedges of hypernode u that are not single-pin community hyperedges.
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
  IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {
    return hypergraph_of_edge(e).communities(e);
  }

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! Can be used to map the global vertex ids to a consecutive range
  // ! of nodes between [0,|V|).
  HypernodeID originalNodeID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).originalNodeID(u);
  }

  // ! Reverse operation of originalNodeID(u)
  HypernodeID globalNodeID(const HypernodeID u) const {
    ASSERT(u < _node_mapping.size(),
           V(u) << V(_node_mapping.size()) << "u is not in the correct range for original node IDs");
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
  HypernodeID communityNodeId(const HypernodeID u) const {
    return hypergraph_of_vertex(u).communityNodeId(u);
  }

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

  bool hasCommunityNodeMapping() const {
    return _community_node_mapping.size() > 0;
  }

  // ! Numa node to which community is assigned to
  PartitionID communityNumaNode(const PartitionID community_id) const {
    ASSERT(static_cast<size_t>(community_id) < _community_node_mapping.size());
    return _community_node_mapping[community_id];
  }

  // ! Sets the community to numa node mapping
  void setCommunityNodeMapping(parallel::scalable_vector<PartitionID>&& community_node_mapping) {
    _community_node_mapping = std::move(community_node_mapping);
  }

  // ! Returns a copy of community to numa node mapping
  parallel::scalable_vector<PartitionID> communityNodeMapping() const {
    return _community_node_mapping;
  }

  // ####################### Contract / Uncontract #######################

  Memento contract(const HypernodeID, const HypernodeID) {
    ERROR("contract(u,v) is not supported in numa hypergraph");
    return Memento();
  }

  Memento contract(const HypernodeID, const HypernodeID, const PartitionID) {
    ERROR("contract(u,v,c) is not supported in numa hypergraph");
    return Memento();
  }

  /*!
   * Contracts a given community structure. All vertices with the same label
   * are collapsed into the same vertex. The resulting single-pin and parallel
   * hyperedges are removed from the contracted graph. The function returns
   * the contracted hypergraph and a mapping which specifies a mapping from
   * community label (given in 'communities') to a vertex in the coarse hypergraph.
   *
   * \param communities Community structure that should be contracted
   * \param task_group_id Task Group ID
   */
  NumaHypergraph contract(parallel::scalable_vector<HypernodeID>& communities,
                          const TaskGroupID task_group_id) const {
    using Hypernode = typename Hypergraph::Hypernode;
    using Hyperedge = typename Hypergraph::Hyperedge;
    ASSERT(communities.size() == _num_hypernodes);
    const int num_numa_nodes = _hypergraphs.size();
    ASSERT(TBBNumaArena::instance().num_used_numa_nodes() == num_numa_nodes);

    parallel::scalable_vector<TmpContractionBuffer*> tmp_contraction_buffer;
    for ( int node = 0; node < num_numa_nodes; ++node ) {
      ASSERT(_hypergraphs[node]._tmp_contraction_buffer);
      tmp_contraction_buffer.push_back(_hypergraphs[node]._tmp_contraction_buffer);

      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_hypernodes) <=
        tmp_contraction_buffer.back()->tmp_hypernodes.size());
      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_hypernodes) <=
        tmp_contraction_buffer.back()->tmp_num_incident_nets.size());
      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_hypernodes) <=
        tmp_contraction_buffer.back()->hn_weights.size());
      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_hyperedges) <=
        tmp_contraction_buffer.back()->tmp_hyperedges.size());
      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_pins) <=
        tmp_contraction_buffer.back()->tmp_incidence_array.size());
      ASSERT(static_cast<size_t>(_hypergraphs[node]._num_hyperedges) <=
        tmp_contraction_buffer.back()->valid_hyperedges.size());
    }

    // #################### STAGE 1 ####################
    // Remapping of vertex ids
    utils::Timer::instance().start_timer("preprocess_contractions", "Preprocess Contractions");
    parallel::scalable_vector<HypernodeID> mapping(_num_hypernodes, kInvalidHypernode);
    parallel::scalable_vector<HypernodeID> num_numa_hypernodes_prefix_sum(num_numa_nodes + 1, 0);
    for ( HypernodeID id = 0; id < _num_hypernodes; ++id ) {
      const HypernodeID hn = globalNodeID(id);
      if ( nodeIsEnabled(hn) ) {
        HypernodeID community = communities[id];
        if ( mapping[community] == kInvalidHypernode ) {
          const int node = common::get_numa_node_of_vertex(hn);
          // Setup mapping from community id to a global node id in the
          // contracted hypergraph. Note, a global node id encodes the
          // position in the hypernode vector (in its corresponding streaming
          // hypergraph) and the numa node id.
          mapping[community] =  common::get_global_vertex_id(
            node, num_numa_hypernodes_prefix_sum[node + 1]++);
          communities[id] = mapping[community];
        } else {
          communities[id] = mapping[community];
        }
      } else {
        communities[id] = kInvalidHypernode;
      }
    }
    for ( int node = 1; node <= num_numa_nodes; ++node ) {
      num_numa_hypernodes_prefix_sum[node] += num_numa_hypernodes_prefix_sum[node - 1];
    }

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
      return communities[originalNodeID(hn)];
    };

    // Mapping from a coarse vertex id of the coarse hypergraph to its
    // original id in the coarse hypergraph
    auto map_to_original_id_in_coarse_hypergraph = [&](const HypernodeID coarse_hn) {
      const int node = common::get_numa_node_of_vertex(coarse_hn);
      return num_numa_hypernodes_prefix_sum[node] +
        common::get_local_position_of_vertex(coarse_hn);
    };

    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
        const HypernodeID num_hypernodes_on_numa_node =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        tbb::parallel_invoke([&] {
          tmp_contraction_buffer[node]->hn_weights.assign(num_hypernodes_on_numa_node,
            parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
        }, [&] {
          tmp_contraction_buffer[node]->tmp_hypernodes.assign(num_hypernodes_on_numa_node, Hypernode(true));
        }, [&] {
          tmp_contraction_buffer[node]->tmp_num_incident_nets.assign(num_hypernodes_on_numa_node,
            parallel::IntegralAtomicWrapper<size_t>(0));
        });
      });

    doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
      const int node = common::get_numa_node_of_vertex(coarse_hn);
      const HypernodeID local_id = common::get_local_position_of_vertex(coarse_hn);
      ASSERT(static_cast<size_t>(node) < tmp_contraction_buffer.size());
      ASSERT(local_id < (num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node]));
      // Weight vector is atomic => thread-safe
      tmp_contraction_buffer[node]->hn_weights[local_id] += nodeWeight(hn);
      // In case community detection is enabled all vertices matched to one vertex
      // in the contracted hypergraph belong to same community. Otherwise, all communities
      // are default assigned to community 0
      tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setCommunityID(communityID(hn));
      // Aggregate upper bound for number of incident nets of the contracted vertex
      tmp_contraction_buffer[node]->tmp_num_incident_nets[local_id] += nodeDegree(hn);
    });
    utils::Timer::instance().stop_timer("preprocess_contractions");

    // #################### STAGE 2 ####################
    // In this step hyperedges and incident nets of vertices are contracted inside the temporary
    // buffers. The vertex ids of pins are already remapped to the vertex ids in the coarse
    // graph and duplicates are removed. Also nets that become single-pin hyperedges are marked
    // as invalid. All incident nets of vertices that are collapsed into one vertex in the coarse
    // graph are also aggregate in a consecutive memory range and duplicates are removed. Note
    // that parallel and single-pin hyperedges are not removed from the incident nets (will be done
    // in a postprocessing step).
    utils::Timer::instance().start_timer("contract_incidence_structure", "Contract Incidence Structures");
    parallel::scalable_vector<Vector<HyperedgeID>> tmp_incident_nets(num_numa_nodes);
    ConcurrentBucketMap<HyperedgeHash> hyperedge_hash_map;
    hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);
    tbb::parallel_invoke([&] {
      // Contract Hyperedges
      utils::Timer::instance().start_timer("contract_hyperedges", "Contract Hyperedges", true);
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
        auto& tmp_hyperedges = tmp_contraction_buffer[node]->tmp_hyperedges;
        auto& tmp_incidence_array = tmp_contraction_buffer[node]->tmp_incidence_array;
        auto& valid_hyperedges = tmp_contraction_buffer[node]->valid_hyperedges;
        tbb::parallel_for(ID(0), _hypergraphs[node].initialNumEdges(), [&, node](const HyperedgeID& local_id) {
          const HyperedgeID he = common::get_global_edge_id(node, local_id);
          if ( edgeIsEnabled(he) ) {
            // Copy hyperedge and pins to temporary buffer
            const Hyperedge& e = _hypergraphs[node]._hyperedges[local_id];
            ASSERT(static_cast<size_t>(local_id) < tmp_hyperedges.size());
            ASSERT(e.firstEntry() <= tmp_incidence_array.size());
            tmp_hyperedges[local_id] = e;
            valid_hyperedges[local_id] = 1;

            // Map pins to vertex ids in coarse graph
            const size_t incidence_array_start = e.firstEntry();
            const size_t incidence_array_end = e.firstInvalidEntry();
            for ( size_t pos = incidence_array_start; pos < incidence_array_end; ++pos ) {
              const HypernodeID pin = _hypergraphs[node]._incidence_array[pos];
              ASSERT(pos < tmp_incidence_array.size());
              tmp_incidence_array[pos] = map_to_coarse_hypergraph(pin);
            }

            // Remove duplicates and disabled vertices
            auto first_entry_it = tmp_incidence_array.begin() + incidence_array_start;
            std::sort(first_entry_it, tmp_incidence_array.begin() + incidence_array_end);
            auto first_invalid_entry_it = std::unique(first_entry_it,
              tmp_incidence_array.begin() + incidence_array_end);
            while ( first_entry_it != first_invalid_entry_it && *(first_invalid_entry_it - 1) == kInvalidHypernode ) {
              --first_invalid_entry_it;
            }

            // Update size of hyperedge in temporary hyperedge buffer
            const size_t contracted_size = std::distance(
              tmp_incidence_array.begin() + incidence_array_start, first_invalid_entry_it);
            tmp_hyperedges[local_id].setSize(contracted_size);

            if ( contracted_size > 1 ) {
              // Compute hash of contracted hyperedge
              size_t he_hash = kEdgeHashSeed;
              for ( size_t pos = incidence_array_start; pos < incidence_array_start + contracted_size; ++pos ) {
                he_hash += kahypar::math::hash(tmp_incidence_array[pos]);
              }
              hyperedge_hash_map.insert(he_hash,
                HyperedgeHash { he, he_hash, contracted_size, true });
            } else {
              // Hyperedge becomes a single-pin hyperedge
              valid_hyperedges[local_id] = 0;
              tmp_hyperedges[local_id].disable();
            }
          } else {
            valid_hyperedges[local_id] = 0;
          }
        });
      });
      utils::Timer::instance().stop_timer("contract_hyperedges");
    }, [&] {
      // Contract Incident Nets
      utils::Timer::instance().start_timer("tmp_contract_incident_nets", "Tmp Contract Incident Nets", true);

      // Compute start position the incident nets of a coarse vertex in the
      // temporary incident nets array with a parallel prefix sum
      parallel::scalable_vector<Vector<parallel::IntegralAtomicWrapper<size_t>>> tmp_incident_nets_pos(num_numa_nodes);
      parallel::scalable_vector<parallel::TBBPrefixSum<
        parallel::IntegralAtomicWrapper<size_t>>> tmp_incident_nets_prefix_sum;
      for ( int node = 0; node < num_numa_nodes; ++node ) {
        tmp_incident_nets_prefix_sum.emplace_back(tmp_contraction_buffer[node]->tmp_num_incident_nets);
      }
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
        const HypernodeID num_hypernodes_on_numa_node =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        tbb::parallel_invoke([&] {
          tbb::parallel_scan(tbb::blocked_range<size_t>(
            0UL, UI64(num_hypernodes_on_numa_node)), tmp_incident_nets_prefix_sum[node]);
          const HyperedgeID total_degree_on_numa_node = tmp_incident_nets_prefix_sum[node].total_sum();
          tmp_incident_nets[node].assign(total_degree_on_numa_node, 0);
        }, [&] {
          tmp_incident_nets_pos[node].assign(num_hypernodes_on_numa_node,
            parallel::IntegralAtomicWrapper<size_t>(0));
        });
      });

      // Write the incident nets of each contracted vertex to the temporary incident net array
      doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
        const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
        const int coarse_node = common::get_numa_node_of_vertex(coarse_hn);
        const HypernodeID coarse_local_id = common::get_local_position_of_vertex(coarse_hn);
        const int node = common::get_numa_node_of_vertex(hn);
        const HypernodeID local_id = common::get_local_position_of_vertex(hn);
        const HyperedgeID node_degree = nodeDegree(hn);
        ASSERT(coarse_node < num_numa_nodes);
        ASSERT(coarse_local_id < tmp_incident_nets_pos[coarse_node].size());
        size_t incident_nets_pos = tmp_incident_nets_prefix_sum[coarse_node][coarse_local_id] +
          tmp_incident_nets_pos[coarse_node][coarse_local_id].fetch_add(node_degree);
        ASSERT(incident_nets_pos + node_degree <= tmp_incident_nets_prefix_sum[coarse_node][coarse_local_id + 1]);
        memcpy(tmp_incident_nets[coarse_node].data() + incident_nets_pos,
               _hypergraphs[node]._incident_nets.data() +
               _hypergraphs[node]._hypernodes[local_id].firstEntry(),
               sizeof(HyperedgeID) * node_degree);
      });

      // Setup temporary hypernodes
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
        std::mutex high_degree_vertex_mutex;
        parallel::scalable_vector<HypernodeID> high_degree_vertices;
        const HypernodeID num_hypernodes_on_numa_node =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        tbb::parallel_for(ID(0), num_hypernodes_on_numa_node, [&, node](const HypernodeID& local_id) {
          // Remove duplicates
          const size_t incident_nets_start = tmp_incident_nets_prefix_sum[node][local_id];
          const size_t incident_nets_end = tmp_incident_nets_prefix_sum[node][local_id + 1];
          const size_t tmp_degree = incident_nets_end - incident_nets_start;
          const HypernodeID coarse_hn = common::get_global_vertex_id(node, local_id);
          const HypernodeID original_coarse_id = map_to_original_id_in_coarse_hypergraph(coarse_hn);
          if ( tmp_degree <= HIGH_DEGREE_CONTRACTION_THRESHOLD ) {
            std::sort(tmp_incident_nets[node].begin() + incident_nets_start,
                      tmp_incident_nets[node].begin() + incident_nets_end);
            auto first_invalid_entry_it = std::unique(tmp_incident_nets[node].begin() + incident_nets_start,
                                                      tmp_incident_nets[node].begin() + incident_nets_end);


            // Setup pointers to temporary incident nets
            const size_t contracted_size = std::distance(
              tmp_incident_nets[node].begin() + incident_nets_start, first_invalid_entry_it);
            tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setSize(contracted_size);
          } else {
            std::lock_guard<std::mutex> lock(high_degree_vertex_mutex);
            high_degree_vertices.push_back(local_id);
          }
          tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setWeight(
            tmp_contraction_buffer[node]->hn_weights[local_id]);
          tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setFirstEntry(incident_nets_start);
          tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setOriginalNodeID(original_coarse_id);
        });

        if ( !high_degree_vertices.empty() ) {
          // High degree vertices are treated special, because sorting and afterwards
          // removing duplicates can become a major sequential bottleneck. Therefore,
          // we distribute the incident nets of a high degree vertex into our concurrent
          // bucket map. As a result all equal incident nets reside in the same bucket
          // afterwards. In a second step, we process each bucket in parallel and apply
          // for each bucket the duplicate removal procedure from above.
          ConcurrentBucketMap<HyperedgeID> duplicate_incident_nets_map;
          for ( const HypernodeID& local_id : high_degree_vertices ) {
            const size_t incident_nets_start = tmp_incident_nets_prefix_sum[node][local_id];
            const size_t incident_nets_end = tmp_incident_nets_prefix_sum[node][local_id + 1];
            const size_t tmp_degree = incident_nets_end - incident_nets_start;

            // Insert incident nets into concurrent bucket map
            duplicate_incident_nets_map.reserve_for_estimated_number_of_insertions(tmp_degree);
            tbb::parallel_for(incident_nets_start, incident_nets_end, [&](const size_t pos) {
              HyperedgeID he = tmp_incident_nets[node][pos];
              duplicate_incident_nets_map.insert(he, std::move(he));
            });

            // Process each bucket in parallel and remove duplicates
            std::atomic<size_t> incident_nets_pos(incident_nets_start);
            tbb::parallel_for(0UL, duplicate_incident_nets_map.numBuckets(), [&](const size_t bucket) {
              auto& incident_net_bucket = duplicate_incident_nets_map.getBucket(bucket);
              std::sort(incident_net_bucket.begin(), incident_net_bucket.end());
              auto first_invalid_entry_it = std::unique(incident_net_bucket.begin(), incident_net_bucket.end());
              const size_t bucket_degree = std::distance(incident_net_bucket.begin(), first_invalid_entry_it);
              const size_t tmp_incident_nets_pos = incident_nets_pos.fetch_add(bucket_degree);
              memcpy(tmp_incident_nets[node].data() + tmp_incident_nets_pos,
                    incident_net_bucket.data(), sizeof(HyperedgeID) * bucket_degree);
              duplicate_incident_nets_map.clear(bucket);
            });

            // Update number of incident nets of high degree vertex
            const size_t contracted_size = incident_nets_pos.load() - incident_nets_start;
            tmp_contraction_buffer[node]->tmp_hypernodes[local_id].setSize(contracted_size);
          }
        }
      });
      utils::Timer::instance().stop_timer("tmp_contract_incident_nets");
    });
    utils::Timer::instance().stop_timer("contract_incidence_structure");

    // #################### STAGE 3 ####################
    // In the step before we aggregated hyperedges within a bucket data structure.
    // Hyperedges with the same hash/footprint are stored inside the same bucket.
    // We iterate now in parallel over each bucket and sort each bucket
    // after its hash. A bucket is processed by one thread and parallel
    // hyperedges are detected by comparing the pins of hyperedges with
    // the same hash.
    utils::Timer::instance().start_timer("remove_parallel_hyperedges", "Remove Parallel Hyperedges");

    // Helper function that checks if two hyperedges are parallel
    // Note, pins inside the hyperedges are sorted.
    auto check_if_hyperedges_are_parallel = [&](const HyperedgeID lhs,
                                                const HyperedgeID rhs) {
      const int node_lhs = common::get_numa_node_of_edge(lhs);
      const int node_rhs = common::get_numa_node_of_edge(rhs);
      const HyperedgeID local_id_lhs = common::get_local_position_of_edge(lhs);
      const HyperedgeID local_id_rhs = common::get_local_position_of_edge(rhs);
      const Hyperedge& lhs_he = tmp_contraction_buffer[node_lhs]->tmp_hyperedges[local_id_lhs];
      const Hyperedge& rhs_he = tmp_contraction_buffer[node_rhs]->tmp_hyperedges[local_id_rhs];
      auto& tmp_incidence_array_lhs = tmp_contraction_buffer[node_lhs]->tmp_incidence_array;
      auto& tmp_incidence_array_rhs = tmp_contraction_buffer[node_rhs]->tmp_incidence_array;
      if ( lhs_he.size() == rhs_he.size() ) {
        const size_t lhs_start = lhs_he.firstEntry();
        const size_t rhs_start = rhs_he.firstEntry();
        for ( size_t i = 0; i < lhs_he.size(); ++i ) {
          const size_t lhs_pos = lhs_start + i;
          const size_t rhs_pos = rhs_start + i;
          if ( tmp_incidence_array_lhs[lhs_pos] != tmp_incidence_array_rhs[rhs_pos] ) {
            return false;
          }
        }
        return true;
      } else {
        return false;
      }
    };

    tbb::parallel_for(0UL, hyperedge_hash_map.numBuckets(), [&](const size_t bucket) {
      auto& hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
      std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
        [&](const HyperedgeHash& lhs, const HyperedgeHash& rhs) {
          return lhs.hash < rhs.hash || (lhs.hash == rhs.hash && lhs.size < rhs.size);
        });

      // Parallel Hyperedge Detection
      for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
        HyperedgeHash& contracted_he_lhs = hyperedge_bucket[i];
        if ( contracted_he_lhs.valid ) {
          const HyperedgeID lhs_he = contracted_he_lhs.he;
          const int node_lhs = common::get_numa_node_of_edge(lhs_he);
          const HyperedgeID local_id_lhs = common::get_local_position_of_edge(lhs_he);
          HyperedgeWeight lhs_weight = tmp_contraction_buffer[node_lhs]->tmp_hyperedges[local_id_lhs].weight();
          for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
            HyperedgeHash& contracted_he_rhs = hyperedge_bucket[j];
            const HyperedgeID rhs_he = contracted_he_rhs.he;
            const int node_rhs = common::get_numa_node_of_edge(rhs_he);
            const HyperedgeID local_id_rhs = common::get_local_position_of_edge(rhs_he);
            if ( contracted_he_rhs.valid &&
                 contracted_he_lhs.hash == contracted_he_rhs.hash &&
                 check_if_hyperedges_are_parallel(lhs_he, rhs_he) ) {
                // Hyperedges are parallel
                lhs_weight += tmp_contraction_buffer[node_rhs]->tmp_hyperedges[local_id_rhs].weight();
                contracted_he_rhs.valid = false;
                tmp_contraction_buffer[node_rhs]->valid_hyperedges[local_id_rhs] = false;
            } else if ( contracted_he_lhs.hash != contracted_he_rhs.hash ) {
              // In case, hash of both are not equal we go to the next hyperedge
              // because we compared it with all hyperedges that had an equal hash
              break;
            }
          }
          tmp_contraction_buffer[node_lhs]->tmp_hyperedges[local_id_lhs].setWeight(lhs_weight);
        }
      }
      hyperedge_hash_map.free(bucket);
    });
    utils::Timer::instance().stop_timer("remove_parallel_hyperedges");

    // #################### STAGE 4 ####################
    // Coarsened hypergraph is constructed here by writting data from temporary
    // buffers to corresponding members in coarsened hypergraph. For the
    // incidence array, we compute a prefix sum over the hyperedge sizes in
    // the contracted hypergraph which determines the start position of the pins
    // of each net in the incidence array. Furthermore, we postprocess the incident
    // nets of each vertex by removing invalid hyperedges and remapping hyperedge ids.
    // Incident nets are also written to the incident nets array with the help of a prefix
    // sum over the node degrees.
    utils::Timer::instance().start_timer("contract_hypergraph", "Contract Hypergraph");

    NumaHypergraph hypergraph;
    // Allocate empty hypergraphs on each NUMA node
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        hypergraph._hypergraphs.emplace_back();
        hypergraph._hypergraphs.back()._node = node;
      });

    // Compute number of hyperedges in coarse graph (those flagged as valid)
    parallel::scalable_vector<parallel::TBBPrefixSum<size_t>> he_mapping;
    parallel::scalable_vector<HyperedgeID> num_numa_hyperedges_prefix_sum(num_numa_nodes + 1, 0);
    for ( int node = 0; node < num_numa_nodes; ++node ) {
      he_mapping.emplace_back(tmp_contraction_buffer[node]->valid_hyperedges);
    }
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
      tbb::parallel_invoke([&] {
        tbb::parallel_scan(tbb::blocked_range<size_t>(
          0UL, UI64(_hypergraphs[node]._num_hyperedges)), he_mapping[node]);
      }, [&] {
        const HypernodeID num_hypernodes_on_numa_node =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        hypergraph._hypergraphs[node]._hypernodes.resize(num_hypernodes_on_numa_node);
      });
    });
    for ( int node = 1; node <= num_numa_nodes; ++node ) {
      num_numa_hyperedges_prefix_sum[node] =
        num_numa_hyperedges_prefix_sum[node - 1] + he_mapping[node - 1].total_sum();
    }
    hypergraph._num_hypernodes = num_numa_hypernodes_prefix_sum[num_numa_nodes];
    hypergraph._num_hyperedges = num_numa_hyperedges_prefix_sum[num_numa_nodes];

    // Mapping from a numa node and local edge id of the coarse hypergraph to its
    // original edge id in the coarse hypergraph
    auto map_to_original_edge_id_in_coarse_hypergraph = [&](const int node, const HyperedgeID local_id) {
      return num_numa_hyperedges_prefix_sum[node] + local_id;
    };

    hypergraph._node_mapping.resize(hypergraph._num_hypernodes);
    hypergraph._edge_mapping.resize(hypergraph._num_hyperedges);
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
      tbb::parallel_invoke([&, node] {
        utils::Timer::instance().start_timer("setup_hyperedges", "Setup Hyperedges", true);
        auto& tmp_hyperedges = tmp_contraction_buffer[node]->tmp_hyperedges;
        auto& tmp_incidence_array = tmp_contraction_buffer[node]->tmp_incidence_array;

        utils::Timer::instance().start_timer("compute_he_pointer", "Compute HE Pointer", true);
        // Compute start position of each hyperedge in incidence array
        const HyperedgeID num_hyperedges_on_numa_node =
          num_numa_hyperedges_prefix_sum[node + 1] - num_numa_hyperedges_prefix_sum[node];
        parallel::scalable_vector<size_t> he_sizes;
        parallel::TBBPrefixSum<size_t> num_pins_prefix_sum(he_sizes);
        tbb::parallel_invoke([&] {
          he_sizes.assign(_hypergraphs[node]._num_hyperedges, 0);
          tbb::parallel_for(ID(0), _hypergraphs[node]._num_hyperedges, [&](const HyperedgeID& local_id) {
            if ( he_mapping[node].value(local_id) /* valid hyperedge contained in coarse graph */ ) {
              he_sizes[local_id] = tmp_hyperedges[local_id].size();
            }
          });

          tbb::parallel_scan(tbb::blocked_range<size_t>(
            0UL, UI64(_hypergraphs[node]._num_hyperedges)), num_pins_prefix_sum);
          const HypernodeID num_pins_on_numa_node = num_pins_prefix_sum.total_sum();
          hypergraph._hypergraphs[node]._num_hyperedges = num_hyperedges_on_numa_node;
          hypergraph._hypergraphs[node]._num_pins = num_pins_on_numa_node;
          hypergraph._hypergraphs[node]._incidence_array.resize(num_pins_on_numa_node);
        }, [&] {
          hypergraph._hypergraphs[node]._hyperedges.resize(num_hyperedges_on_numa_node);
        });
        utils::Timer::instance().stop_timer("compute_he_pointer");

        utils::Timer::instance().start_timer("setup_incidence_array", "Setup Incidence Array", true);
        // Write hyperedges from temporary buffers to incidence array
        tbb::parallel_for(ID(0), _hypergraphs[node]._num_hyperedges, [&](const HyperedgeID& local_id) {
          if ( he_mapping[node].value(local_id) /* hyperedge is valid */ ) {
            const size_t he_pos = he_mapping[node][local_id];
            const HyperedgeID original_he_id = map_to_original_edge_id_in_coarse_hypergraph(node, he_pos);
            const size_t incidence_array_start = num_pins_prefix_sum[local_id];
            Hyperedge& e = hypergraph._hypergraphs[node]._hyperedges[he_pos];
            e = std::move(tmp_hyperedges[local_id]);
            const size_t tmp_incidence_array_start = e.firstEntry();
            std::memcpy(hypergraph._hypergraphs[node]._incidence_array.data() + incidence_array_start,
                        tmp_incidence_array.data() + tmp_incidence_array_start,
                        sizeof(HypernodeID) * e.size());
            e.setFirstEntry(incidence_array_start);
            e.setOriginalEdgeID(original_he_id);
            ASSERT(original_he_id < hypergraph._edge_mapping.size());
            hypergraph._edge_mapping[original_he_id] = common::get_global_edge_id(node, he_pos);
          }
        });
        utils::Timer::instance().stop_timer("setup_incidence_array");
        utils::Timer::instance().stop_timer("setup_hyperedges");
      }, [&, node] {
        utils::Timer::instance().start_timer("setup_hypernodes", "Setup Hypernodes", true);
        auto& tmp_hypernodes = tmp_contraction_buffer[node]->tmp_hypernodes;
        auto& tmp_numa_incident_nets = tmp_incident_nets[node];

        utils::Timer::instance().start_timer("compute_num_incident_nets", "Compute Num Incident Nets", true);
        // Remap hyperedge ids in temporary incident nets to hyperedge ids of the
        // coarse hypergraph and remove singple-pin/parallel hyperedges.
        const HypernodeID num_hypernodes_on_numa_node =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        parallel::scalable_vector<size_t> incident_nets_sizes(num_hypernodes_on_numa_node, 0);
        tbb::parallel_for(ID(0), num_hypernodes_on_numa_node, [&](const HypernodeID& local_id) {
          const size_t incident_nets_start =  tmp_hypernodes[local_id].firstEntry();
          size_t incident_nets_end = tmp_hypernodes[local_id].firstInvalidEntry();
          for ( size_t pos = incident_nets_start; pos < incident_nets_end; ++pos ) {
            const HyperedgeID he = tmp_numa_incident_nets[pos];
            const int node = common::get_numa_node_of_edge(he);
            const HyperedgeID local_he_id = common::get_local_position_of_edge(he);
            if ( he_mapping[node].value(local_he_id) /* hyperedge is valid */ ) {
              tmp_numa_incident_nets[pos] = common::get_global_edge_id(node, he_mapping[node][local_he_id]);
            } else {
              std::swap(tmp_numa_incident_nets[pos--], tmp_numa_incident_nets[--incident_nets_end]);
            }
          }
          const size_t incident_nets_size = incident_nets_end - incident_nets_start;
          tmp_hypernodes[local_id].setSize(incident_nets_size);
          incident_nets_sizes[local_id] = incident_nets_size;
        });

        // Compute start position of the incident nets for each vertex inside
        // the coarsened incident net array
        parallel::TBBPrefixSum<size_t> num_incident_nets_prefix_sum(incident_nets_sizes);
        tbb::parallel_scan(tbb::blocked_range<size_t>(
          0UL, UI64(num_hypernodes_on_numa_node)), num_incident_nets_prefix_sum);
        const size_t total_degree_on_numa_node = num_incident_nets_prefix_sum.total_sum();
        hypergraph._hypergraphs[node]._num_hypernodes = num_hypernodes_on_numa_node;
        hypergraph._hypergraphs[node]._total_degree = total_degree_on_numa_node;
        hypergraph._hypergraphs[node]._incident_nets.resize(total_degree_on_numa_node);
        utils::Timer::instance().stop_timer("compute_num_incident_nets");

        utils::Timer::instance().start_timer("setup_incident_nets", "Setup Incidenct Nets", true);
        // Write incident nets from temporary buffer to incident nets array
        tbb::parallel_for(ID(0), num_hypernodes_on_numa_node, [&](const HypernodeID& local_id) {
          const size_t incident_nets_start = num_incident_nets_prefix_sum[local_id];
          Hypernode& hn = hypergraph._hypergraphs[node]._hypernodes[local_id];
          hn = std::move(tmp_hypernodes[local_id]);
          const size_t tmp_incident_nets_start = hn.firstEntry();
          std::memcpy(hypergraph._hypergraphs[node]._incident_nets.data() + incident_nets_start,
                      tmp_numa_incident_nets.data() + tmp_incident_nets_start,
                      sizeof(HyperedgeID) * hn.size());
          hn.setFirstEntry(incident_nets_start);
          const HypernodeID original_hn_id = hn.originalNodeID();
          ASSERT(original_hn_id < hypergraph._node_mapping.size());
          hypergraph._node_mapping[original_hn_id] = common::get_global_vertex_id(node, local_id);
        });
        utils::Timer::instance().stop_timer("setup_incident_nets");
        utils::Timer::instance().stop_timer("setup_hypernodes");
      });
    });
    utils::Timer::instance().stop_timer("contract_hypergraph");

    // Compute stats of NUMA hypergraph
    for ( int node = 0; node < num_numa_nodes; ++node ) {
      hypergraph._num_pins += hypergraph.initialNumPins(node);
      hypergraph._total_degree += hypergraph.initialTotalVertexDegree(node);
    }

    // Initialize Communities and Update Total Weight
    utils::Timer::instance().start_timer("setup_communities", "Setup Communities");
    tbb::parallel_invoke([&] {
      const Hypergraph& hg = _hypergraphs[0];
      if ( hg._community_support.isInitialized() ) {
        hypergraph.initializeCommunities(task_group_id);
        if ( hg._community_support.areCommunityHyperedgesInitialized() ) {
          hypergraph.initializeCommunityHyperedges(task_group_id);
        }
      }
    }, [&] {
      hypergraph.updateTotalWeight(task_group_id);
      hypergraph._community_node_mapping = _community_node_mapping;
    });
    utils::Timer::instance().stop_timer("setup_communities");

    // Pass contraction buffer to contracted hypergraph
    for ( int node = 0; node < num_numa_nodes; ++node ) {
      hypergraph._hypergraphs[node]._tmp_contraction_buffer = _hypergraphs[node]._tmp_contraction_buffer;
    }
    return hypergraph;
  }

  void uncontract(const Memento&, parallel::scalable_vector<HyperedgeID>&) {
    ERROR("uncontract(memento,parallel_he) is not supported in numa hypergraph");
  }

  void uncontract(const std::vector<Memento>&,
                  parallel::scalable_vector<HyperedgeID>&,
                  const kahypar::ds::FastResetFlagArray<>&,
                  const bool) {
    ERROR("uncontract(...) is not supported in numa hypergraph");
  }

  void restoreDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("restoreDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in numa hypergraph");
  }

  parallel::scalable_vector<HyperedgeID> findDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("findDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in numa hypergraph");
    return parallel::scalable_vector<HyperedgeID>();
  }

  // ####################### Remove / Restore Hyperedges #######################

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).removeIncidentEdgeFromHypernode(he, pin);
    }
    ++_num_removed_hyperedges;
    hypergraph_of_edge(he).disableHyperedge(he);
  }

  void removeSinglePinCommunityEdge(const HyperedgeID, const PartitionID) {
    ERROR("removeSinglePinCommunityEdge(e,c) is not supported in numa hypergraph");
  }

  void removeParallelEdge(const HyperedgeID, const PartitionID) {
    ERROR("removeParallelEdge(e,c) is not supported in numa hypergraph");
  }

  // ! Restores an hyperedge of a certain size.
  void restoreEdge(const HyperedgeID he, const size_t,
                   const HyperedgeID representative = kInvalidHyperedge) {
    unused(representative);
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "already enabled");
    hypergraph_of_edge(he).enableHyperedge(he);
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_vertex(pin).insertIncidentEdgeToHypernode(he, pin);
    }
  }

  // ! Restores a single-pin hyperedge
  void restoreSinglePinHyperedge(const HyperedgeID he) {
    restoreEdge(he, 1);
  }

  void restoreParallelHyperedge(const HyperedgeID,
                                const Memento&,
                                parallel::scalable_vector<HyperedgeID>&,
                                const kahypar::ds::FastResetFlagArray<>* batch_hypernodes = nullptr) {
    unused(batch_hypernodes);
    ERROR("restoreParallelHyperedge(...) is not supported in numa hypergraph");
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

    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].finalizeCommunityNodeIds(_hypergraphs, task_group_id);
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

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   *
   * In case, template parameter Hypergraph is the dynamic hypergraph:
   * The pins of the original hyperedge are sorted in decreasing order of their
   * contraction index. The contraction index of a vertex v is defined as the index
   * of the contraction (u,v) in the contraction history, where v occurs as contraction
   * partner. This is done to fullfil the invariants required by the uncontraction method.
   *
   * Note this function have to be called after parallel community coarsening such
   * that uncontractions can be performed correctly.
   */
  void removeCommunityHyperedges(const TaskGroupID task_group_id) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].removeCommunityHyperedges(task_group_id, {}, _hypergraphs);
        });
  }

  void buildContractionHierarchy(const std::vector<Memento>&) {
    ERROR("buildContractionHierarchy(mementos) is not supported in numa hypergraph");
  }

  void invalidateDisabledHyperedgesFromIncidentNets(const TaskGroupID) {
    ERROR("invalidateDisabledHyperedgesFromIncidentNets(id) is not supported in numa hypergraph");
  }

  // ####################### Copy #######################

  // ! Copy numa hypergraph in parallel
  // ! TODO(heuer): in case dynamic hypergraph is used, vertex and edge ids must
  // ! be also compactified
  NumaHypergraph copy(const TaskGroupID task_group_id) {
    NumaHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    tbb::parallel_invoke([&] {
      hypergraph._hypergraphs.resize(_hypergraphs.size());
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
            hypergraph._hypergraphs[node] = _hypergraphs[node].copy(task_group_id);
          });
    }, [&] {
      hypergraph._node_mapping.resize(_node_mapping.size());
      memcpy(hypergraph._node_mapping.data(), _node_mapping.data(),
        sizeof(HypernodeID) * _node_mapping.size());
    }, [&] {
      hypergraph._edge_mapping.resize(_edge_mapping.size());
      memcpy(hypergraph._edge_mapping.data(), _edge_mapping.data(),
        sizeof(HyperedgeID) * _edge_mapping.size());
    }, [&] {
      hypergraph._community_node_mapping.resize(_community_node_mapping.size());
      memcpy(hypergraph._community_node_mapping.data(), _community_node_mapping.data(),
        sizeof(PartitionID) * _community_node_mapping.size());
    });

    return hypergraph;
  }

  // ! Copy numa hypergraph sequential
  // ! TODO(heuer): in case dynamic hypergraph is used, vertex and edge ids must
  // ! be also compactified
  NumaHypergraph copy() {
    NumaHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    hypergraph._hypergraphs.resize(_hypergraphs.size());
    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      // Note, in case we copy sequential the hypergraphs (which are originally
      // located one on each NUMA node) are copied all to the NUMA node where
      // the calling thread is running.
      hypergraph._hypergraphs[node] = _hypergraphs[node].copy();
    }

    hypergraph._node_mapping.resize(_node_mapping.size());
    memcpy(hypergraph._node_mapping.data(), _node_mapping.data(),
      sizeof(HypernodeID) * _node_mapping.size());
    hypergraph._edge_mapping.resize(_edge_mapping.size());
    memcpy(hypergraph._edge_mapping.data(), _edge_mapping.data(),
      sizeof(HyperedgeID) * _edge_mapping.size());
    hypergraph._community_node_mapping.resize(_community_node_mapping.size());
    memcpy(hypergraph._community_node_mapping.data(), _community_node_mapping.data(),
      sizeof(PartitionID) * _community_node_mapping.size());

    return hypergraph;
  }

  // Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      tbb::parallel_invoke([&] {
        tbb::parallel_for(0UL, _hypergraphs.size(), [&](const size_t i) {
          _hypergraphs[i].freeInternalData();
        }, tbb::static_partitioner());
      }, [&] {
        parallel::parallel_free(_node_mapping,
          _edge_mapping, _community_node_mapping);
      });
    }
    _num_hypernodes = 0;
    _num_hyperedges = 0;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      utils::MemoryTreeNode* numa_hypergraph_node =parent->addChild(
        "NUMA Hypergraph " + std::to_string(node));
      _hypergraphs[node].memoryConsumption(numa_hypergraph_node);
    }

    parent->addChild("Global Node Mapping",
      sizeof(HypernodeID) * _node_mapping.size());
    parent->addChild("Global Edge Mapping",
      sizeof(HyperedgeID) * _edge_mapping.size());
    parent->addChild("Community NUMA Node Mapping",
      sizeof(PartitionID) * _community_node_mapping.size());
  }

 private:
  template <typename Hypgraph,
            typename Factory,
            typename HwTopo,
            typename TBBArena>
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
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
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