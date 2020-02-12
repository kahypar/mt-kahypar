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
#include "mt-kahypar/datastructures/streaming_map.h"
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

 public:
  // ! Type Traits
  using Hypergraph = HyperGraph;
  using HardwareTopology = HwTopology;
  using TBBNumaArena = TBB;

  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = true;
  static constexpr bool is_partitioned = false;

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
  std::pair<NumaHypergraph, parallel::scalable_vector<HypernodeID>> contract(
    const parallel::scalable_vector<HypernodeID>& communities,
    const TaskGroupID task_group_id) const {
    ASSERT(communities.size() == _num_hypernodes);
    const int num_numa_nodes = _hypergraphs.size();
    ASSERT(TBBNumaArena::instance().num_used_numa_nodes() == num_numa_nodes);

    // #################### STAGE 1 ####################
    // Remapping of vertex ids
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
        }
      }
    }
    for ( int node = 1; node <= num_numa_nodes; ++node ) {
      num_numa_hypernodes_prefix_sum[node] += num_numa_hypernodes_prefix_sum[node - 1];
    }
    const HypernodeID num_hypernodes = num_numa_hypernodes_prefix_sum[num_numa_nodes];

    parallel::scalable_vector<HypernodeID> hn_to_numa_node;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
    parallel::scalable_vector<PartitionID> community_ids;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<uint8_t>> is_high_degree_vertex;
    tbb::parallel_invoke([&] {
      hn_to_numa_node.resize(num_hypernodes, 0);
    }, [&] {
      hn_weights.assign(num_hypernodes,
        parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
    }, [&] {
      community_ids.assign(num_hypernodes, 0);
    }, [&] {
      is_high_degree_vertex.assign(num_hypernodes,
        parallel::IntegralAtomicWrapper<uint8_t>(false));
    });

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
      return mapping[communities[originalNodeID(hn)]];
    };

    // Mapping from a coarse vertex id of the coarse hypergraph to its
    // original id in the coarse hypergraph
    auto map_to_original_id_in_coarse_hypergraph = [&](const HypernodeID coarse_hn) {
      const int node = common::get_numa_node_of_vertex(coarse_hn);
      return num_numa_hypernodes_prefix_sum[node] +
        common::get_local_position_of_vertex(coarse_hn);
    };

    doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
      const int node = common::get_numa_node_of_vertex(coarse_hn);
      const HypernodeID original_id = map_to_original_id_in_coarse_hypergraph(coarse_hn);
      ASSERT(original_id < num_hypernodes);
      hn_to_numa_node[original_id] = node;
      // Weight vector is atomic => thread-safe
      hn_weights[original_id] += nodeWeight(hn);
      // In case community detection is enabled all vertices matched to one vertex
      // in the contracted hypergraph belong to same community. Otherwise, all communities
      // are default assigned to community 0
      community_ids[original_id] = communityID(hn);
      // Vector is atomic => thread-safe
      is_high_degree_vertex[original_id].fetch_or(isHighDegreeVertex(hn));
    });

    // #################### STAGE 2 ####################
    // We iterate over all hyperedges in parallel and remap their ids
    // to the ones determined in the step before. Furthermore, duplicates
    // and disabled hyperedges are removed. The hyperedges are then inserted
    // into a streaming map with their hash as key. All hyperedges with the same
    // hash are then present in the same bucket of the streaming map, which
    // makes it possible to detect parallel hyperedges in parallel.
    StreamingMap<size_t, ContractedHyperedge> hash_to_hyperedge;
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].contractHyperedges(
              task_group_id, hash_to_hyperedge, map_to_coarse_hypergraph);
        });

    using HyperedgeMap = parallel::scalable_vector<parallel::scalable_vector<ContractedHyperedge>>;
    HyperedgeMap hyperedge_buckets(hash_to_hyperedge.size());
    hash_to_hyperedge.copy(hyperedge_buckets, [&](const size_t key) {
      return key % hash_to_hyperedge.size();
    });

    // #################### STAGE 3 ####################
    // We iterate now in parallel over each bucket and sort each bucket
    // after its hash. A bucket is processed by one thread and parallel
    // hyperedges are detected by comparing the pins of hyperedges with
    // the same hash.

    // Helper function that checks if two hyperedges are parallel
    // Note, pins inside the hyperedges are sorted.
    auto check_if_hyperedges_are_parallel = [](const parallel::scalable_vector<HypernodeID>& lhs,
                                                const parallel::scalable_vector<HypernodeID>& rhs) {
      HEAVY_COARSENING_ASSERT(std::is_sorted(lhs.cbegin(), lhs.cend()));
      HEAVY_COARSENING_ASSERT(std::is_sorted(rhs.cbegin(), rhs.cend()));
      if ( lhs.size() == rhs.size() ) {
        for ( size_t i = 0; i < lhs.size(); ++i ) {
          if ( lhs[i] != rhs[i] ) {
            return false;
          }
        }
        return true;
      } else {
        return false;
      }
    };

    // Stores the prefix sum over the number of hyperedges in each bucket
    parallel::scalable_vector<HyperedgeID> num_hyperedges_prefix_sum(hyperedge_buckets.size() + 1, 0);
    // Stores the prefix sum over the number of hyperedges in each bucket for each numa node
    parallel::scalable_vector<parallel::scalable_vector<HyperedgeID>> num_numa_hyperedges_prefix_sum(
      hyperedge_buckets.size() + 1, parallel::scalable_vector<HyperedgeID>(num_numa_nodes, 0));
    // Stores the prefix sum over the number of pins in each bucket for each numa node
    parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> num_numa_pins_prefix_sum(
      hyperedge_buckets.size() + 1, parallel::scalable_vector<HypernodeID>(num_numa_nodes, 0));
    // Stores the prefix sum over the number of incident nets for each numa node and vertex
    parallel::scalable_vector<parallel::scalable_vector<
      parallel::IntegralAtomicWrapper<HyperedgeID>>> num_numa_incident_nets(num_numa_nodes);
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        const HypernodeID num_numa_hypernodes =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        num_numa_incident_nets[node].resize(num_numa_hypernodes);
      });

    tbb::parallel_for(0UL, hyperedge_buckets.size(), [&](const size_t bucket) {
      parallel::scalable_vector<ContractedHyperedge>& hyperedge_bucket = hyperedge_buckets[bucket];
      std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
        [&](const ContractedHyperedge& lhs, const ContractedHyperedge& rhs) {
          return lhs.hash < rhs.hash || (lhs.hash == rhs.hash && lhs.hyperedge.size() < rhs.hyperedge.size());
        });

      // Parallel Hyperedge Detection
      for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
        ContractedHyperedge& contracted_he_lhs = hyperedge_bucket[i];
        if ( !contracted_he_lhs.is_parallel ) {
          ASSERT(contracted_he_lhs.node < num_numa_nodes);
          const int node = contracted_he_lhs.node;
          // Determine position for each hyperedge and its pin in the hyperedge vector
          // and incidence array of its corresponding streaming hypergraph
          contracted_he_lhs.he_idx = num_numa_hyperedges_prefix_sum[bucket + 1][node]++;
          contracted_he_lhs.pin_idx = num_numa_pins_prefix_sum[bucket + 1][node];
          num_numa_pins_prefix_sum[bucket + 1][node] += contracted_he_lhs.hyperedge.size();
          ++num_hyperedges_prefix_sum[bucket + 1];
          // Aggregate the number of incident nets of each vertex
          for ( const HypernodeID& pin : contracted_he_lhs.hyperedge ) {
            const int pin_node = common::get_numa_node_of_vertex(pin);
            const HypernodeID local_id = common::get_local_position_of_vertex(pin);
            ASSERT(pin_node < num_numa_nodes);
            ASSERT(local_id < num_numa_incident_nets[pin_node].size());
            ++num_numa_incident_nets[pin_node][local_id];
          }

          for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
            ContractedHyperedge& contracted_he_rhs = hyperedge_bucket[j];
            if ( !contracted_he_rhs.is_parallel &&
                 contracted_he_lhs.hash == contracted_he_rhs.hash &&
                 check_if_hyperedges_are_parallel(
                   contracted_he_lhs.hyperedge, contracted_he_rhs.hyperedge) ) {
                // Hyperedges are parallel
                contracted_he_lhs.weight += contracted_he_rhs.weight;
                contracted_he_rhs.is_parallel = true;
            } else if ( contracted_he_lhs.hash != contracted_he_rhs.hash ) {
              // In case, hash of both are not equal we go to the next hyperedge
              // because we compared it with all hyperedges that had an equal hash
              break;
            }
          }
        }
      }
    });

    for ( size_t i = 1; i <= hyperedge_buckets.size(); ++i ) {
      num_hyperedges_prefix_sum[i] += num_hyperedges_prefix_sum[i - 1];
      for ( int node = 0; node < num_numa_nodes; ++node ) {
        num_numa_hyperedges_prefix_sum[i][node] += num_numa_hyperedges_prefix_sum[i - 1][node];
        num_numa_pins_prefix_sum[i][node] += num_numa_pins_prefix_sum[i - 1][node];
      }
    }

    // #################### STAGE 4 ####################
    // Initialize hypergraph
    NumaHypergraph hypergraph;

    // Prefix sum over the number of incident nets on each NUMA node
    parallel::scalable_vector<parallel::TBBPrefixSum<
      parallel::IntegralAtomicWrapper<HyperedgeID>>> num_numa_incident_nets_prefix_sum;
    // Allocate empty hypergraphs on each NUMA node
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        hypergraph._hypergraphs.emplace_back();
        num_numa_incident_nets_prefix_sum.emplace_back(num_numa_incident_nets[node]);
      });


    // Compute prefix sum over the number of incident nets on each NUMA node
    // => Used to compute positions of the incident nets of a vertex in each
    // NUMA hypergraph
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        const HypernodeID num_numa_hypernodes =
          num_numa_hypernodes_prefix_sum[node + 1] - num_numa_hypernodes_prefix_sum[node];
        tbb::parallel_scan(tbb::blocked_range<HypernodeID>(
          0UL, num_numa_hypernodes), num_numa_incident_nets_prefix_sum[node]);
      });

    // Setup stats of each hypergraph on each NUMA node in parallel
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        Hypergraph& hg = hypergraph._hypergraphs[node];
        hg._node = node;
        hg._num_hypernodes = num_numa_hypernodes_prefix_sum[node + 1] -
          num_numa_hypernodes_prefix_sum[node];
        hg._num_hyperedges = num_numa_hyperedges_prefix_sum.back()[node];
        hg._num_pins = num_numa_pins_prefix_sum.back()[node];
        hg._total_degree = num_numa_incident_nets_prefix_sum[node].total_sum();
        hg._incident_nets.resize(hg._total_degree);
      });

    // Construct hypergraphs on each NUMA node in parallel
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HyperedgeID>> incident_nets_pos(
      num_hypernodes, parallel::IntegralAtomicWrapper<HyperedgeID>(0));
    hypergraph._node_mapping.resize(num_hypernodes);
    hypergraph._edge_mapping.resize(num_hyperedges_prefix_sum.back());
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(
      task_group_id, [&](const int node) {
        Hypergraph& hg = hypergraph._hypergraphs[node];

        tbb::parallel_invoke([&] {
          // Setup hypernodes
          using Hypernode = typename Hypergraph::Hypernode;
          hg._hypernodes.resize(hg._num_hypernodes);
          tbb::parallel_for(0UL, hg._num_hypernodes, [&](const HypernodeID id) {
            const size_t incident_nets_pos = num_numa_incident_nets_prefix_sum[node][id];
            const size_t incident_nets_size = id == 0 ?
              num_numa_incident_nets_prefix_sum[node][id + 1] :
              num_numa_incident_nets_prefix_sum[node][id + 1] -
              num_numa_incident_nets_prefix_sum[node][id];
            ASSERT(incident_nets_pos + incident_nets_size <= hg._total_degree);
            const HypernodeID hn = common::get_global_vertex_id(node, id);
            const HypernodeID original_id = map_to_original_id_in_coarse_hypergraph(hn);
            ASSERT(original_id < num_hypernodes);
            hypergraph._node_mapping[original_id] = hn;
            Hypernode& hn_obj = hg._hypernodes[id];
            hn_obj.enable();
            hn_obj.setFirstEntry(incident_nets_pos);
            hn_obj.setSize(incident_nets_size);
            hn_obj.setOriginalNodeID(original_id);
            hn_obj.setWeight(hn_weights[original_id]);
            hn_obj.setCommunityID(community_ids[original_id]);
            if ( is_high_degree_vertex[original_id] ) {
              hn_obj.markAsHighDegreeVertex();
            }
          });
        }, [&] {
          // Setup hyperedges, incidence and incident nets array
          using Hyperedge = typename Hypergraph::Hyperedge;
          hg._hyperedges.resize(hg._num_hyperedges);
          hg._incidence_array.resize(hg._num_pins);

          tbb::parallel_for(0UL, hyperedge_buckets.size(), [&](const size_t bucket) {
            HyperedgeID original_he_id = num_hyperedges_prefix_sum[bucket];
            for ( ContractedHyperedge contracted_he : hyperedge_buckets[bucket] ) {
              if ( !contracted_he.is_parallel ) {
                if ( contracted_he.node == node ) {
                  const HyperedgeID he = contracted_he.he_idx +
                    num_numa_hyperedges_prefix_sum[bucket][node];
                  const size_t incidence_array_pos = contracted_he.pin_idx +
                    num_numa_pins_prefix_sum[bucket][node];
                  ASSERT(he < hg._num_hyperedges);
                  ASSERT(incidence_array_pos + contracted_he.hyperedge.size() <=
                        hg._incidence_array.size());

                  // Setup hyperedge
                  ASSERT(original_he_id < num_hyperedges_prefix_sum.back());
                  Hyperedge& he_obj = hg._hyperedges[he];
                  he_obj.enable();
                  he_obj.setFirstEntry(incidence_array_pos);
                  he_obj.setSize(contracted_he.hyperedge.size());
                  he_obj.setOriginalEdgeID(original_he_id);
                  he_obj.setWeight(contracted_he.weight);
                  he_obj.hash() = contracted_he.hash;

                  // Copy content of hyperedge to incidence array
                  memcpy(hg._incidence_array.data() + incidence_array_pos,
                        contracted_he.hyperedge.data(), sizeof(HypernodeID) * he_obj.size());

                  // Add hyperedge as incident net to all pins
                  const HyperedgeID global_he = common::get_global_edge_id(node, he);
                  hypergraph._edge_mapping[original_he_id] = global_he;
                  for ( const HypernodeID& pin : hg.pins(global_he) ) {
                    const int pin_node = common::get_numa_node_of_vertex(pin);
                    const HypernodeID local_id = common::get_local_position_of_vertex(pin);
                    const HypernodeID original_pin_id = map_to_original_id_in_coarse_hypergraph(pin);
                    ASSERT(pin_node < num_numa_nodes);
                    ASSERT(local_id < hypergraph._hypergraphs[pin_node]._num_hypernodes);
                    ASSERT(original_pin_id < incident_nets_pos.size());
                    const size_t incident_nets_position =
                      num_numa_incident_nets_prefix_sum[pin_node][local_id] +
                      incident_nets_pos[original_pin_id]++;
                    ASSERT(incident_nets_position <
                          hypergraph._hypergraphs[pin_node]._incident_nets.size());
                    hypergraph._hypergraphs[pin_node]._incident_nets[incident_nets_position] = global_he;
                  }
                }
                ++original_he_id;
              }
            }
          });
        });
      });

    // Compute stats of NUMA hypergraph
    for ( int node = 0; node < num_numa_nodes; ++node ) {
      hypergraph._num_hypernodes += hypergraph.initialNumNodes(node);
      hypergraph._num_hyperedges += hypergraph.initialNumEdges(node);
      hypergraph._num_pins += hypergraph.initialNumPins(node);
      hypergraph._total_degree += hypergraph.initialTotalVertexDegree(node);
    }

    // Initialize Communities and Update Total Weight
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

    return std::make_pair(std::move(hypergraph), std::move(mapping));
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
          _hypergraphs[node].finalizeCommunityNodeIds(_hypergraphs);
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