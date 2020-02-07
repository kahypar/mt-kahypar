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

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/community_support.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// Forward
class StaticHypergraphFactory;

class StaticHypergraph {

  static constexpr bool enable_heavy_assert = false;

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HypernodeID;

    Hypernode() :
      _begin(0),
      _size(0),
      _original_id(kInvalidHypernode),
      _weight(1),
      _community_id(0),
      _is_high_degree_vertex(false),
      _valid(false) { }

    bool isDisabled() const {
      return _valid == false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incident_nets to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    HypernodeID originalNodeID() const {
      return _original_id;
    }

    void setOriginalNodeID(const HypernodeID original_id) {
      _original_id = original_id;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    PartitionID communityID() const {
      ASSERT(!isDisabled());
      return _community_id;
    }

    void setCommunityID(const PartitionID community_id) {
      ASSERT(!isDisabled());
      _community_id = community_id;
    }

    bool isHighDegreeVertex() const {
      return _is_high_degree_vertex;
    }

    void markAsHighDegreeVertex() {
      _is_high_degree_vertex = true;
    }

   private:
    // ! Index of the first element in _incident_nets
    size_t _begin;
    // ! Number of incident nets
    size_t _size;
    // ! Original id of hypernode
    HypernodeID _original_id;
    // ! Hypernode weight
    HyperedgeWeight _weight;
    // ! Community id
    PartitionID _community_id;
    // ! Indicates, wheather this vertex is a high degree node or not
    bool _is_high_degree_vertex;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /**
   * Represents a hyperedge of the hypergraph and contains all information
   * associated with a net (except connectivity information).
   */
  class Hyperedge {
   public:
    using IDType = HyperedgeID;

    Hyperedge() :
      _begin(0),
      _size(0),
      _original_id(kInvalidHyperedge),
      _weight(1),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
    // ! when iterating over the set of all nodes/edges.
    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    bool isDisabled() const {
      return _valid == false;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incidence_array to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    HyperedgeID originalEdgeID() const {
      return _original_id;
    }

    void setOriginalEdgeID(const HyperedgeID original_id) {
      _original_id = original_id;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    size_t& hash() {
      return _hash;
    }

    size_t hash() const {
      return _hash;
    }

    bool operator== (const Hyperedge& rhs) const {
      return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
    }

    bool operator!= (const Hyperedge& rhs) const {
      return _begin != rhs._begin || _size != rhs._size || _weight != rhs._weight;
    }

   private:
    // ! Index of the first element in _incidence_array
    size_t _begin;
    // ! Number of pins
    size_t _size;
    // ! Original id of hyperedge
    HyperedgeID _original_id;
    // ! hyperedge weight
    HyperedgeWeight _weight;
    // ! Hash of pins
    size_t _hash;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /*!
   * Iterator for HypergraphElements (Hypernodes/Hyperedges)
   *
   * The iterator is used in for-each loops over all hypernodes/hyperedges.
   * In order to support iteration over coarsened hypergraphs, this iterator
   * skips over HypergraphElements marked as invalid.
   * Iterating over the set of vertices \f$V\f$ therefore is linear in the
   * size \f$|V|\f$ of the original hypergraph - even if it has been coarsened
   * to much smaller size. The same also holds true for for-each loops over
   * the set of hyperedges.
   *
   * In order to be as generic as possible, the iterator does not expose the
   * internal Hypernode/Hyperedge representations. Instead only handles to
   * the respective elements are returned, i.e. the IDs of the corresponding
   * hypernodes/hyperedges.
   *
   */
  template <typename ElementType>
  class HypergraphElementIterator :
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                         typename ElementType::IDType,   // value_type
                         std::ptrdiff_t,   // difference_type
                         const typename ElementType::IDType*,   // pointer
                         typename ElementType::IDType> {   // reference
   public:
    using IDType = typename ElementType::IDType;

    /*!
     * Construct a HypergraphElementIterator
     * See GenericHypergraph::nodes() or GenericHypergraph::edges() for usage.
     *
     * If start_element is invalid, the iterator advances to the first valid
     * element.
     *
     * \param start_element A pointer to the starting position
     * \param id The index of the element the pointer points to
     * \param max_id The maximum index allowed
     */
    HypergraphElementIterator(const ElementType* start_element, IDType id, IDType max_id) :
      _id(id),
      _max_id(max_id),
      _element(start_element) {
      if (_id != _max_id && _element->isDisabled()) {
        operator++ ();
      }
    }

    // ! Returns the id of the element the iterator currently points to.
    IDType operator* () const {
      return _id;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    HypergraphElementIterator & operator++ () {
      ASSERT(_id < _max_id);
      do {
        ++_id;
        ++_element;
      } while (_id < _max_id && _element->isDisabled());
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    HypergraphElementIterator operator++ (int) {
      HypergraphElementIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const HypergraphElementIterator& rhs) {
      return _id != rhs._id;
    }

    bool operator== (const HypergraphElementIterator& rhs) {
      return _id == rhs._id;
    }

   private:
    // Handle to the HypergraphElement the iterator currently points to
    IDType _id = 0;
    // Maximum allowed index
    IDType _max_id = 0;
    // HypergraphElement the iterator currently points to
    const ElementType* _element = nullptr;
  };

  static_assert(std::is_trivially_copyable<Hypernode>::value, "Hypernode is not trivially copyable");
  static_assert(std::is_trivially_copyable<Hyperedge>::value, "Hyperedge is not trivially copyable");

  using IncidenceArray = parallel::scalable_vector<HypernodeID>;
  using IncidentNets = parallel::scalable_vector<HyperedgeID>;

 public:
  static constexpr bool is_static_hypergraph = true;
  static constexpr bool is_numa_aware = false;
  static constexpr bool is_partitioned = false;

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename IncidenceArray::const_iterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename IncidentNets::const_iterator;
  // ! Iterator to iterate over the set of communities contained in a hyperedge
  using CommunityIterator = typename CommunitySupport<StaticHypergraph>::CommunityIterator;

  explicit StaticHypergraph() :
    _node(0),
    _num_hypernodes(0),
    _num_removed_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_degree(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
    _community_support() { }

  StaticHypergraph(const StaticHypergraph&) = delete;
  StaticHypergraph & operator= (const StaticHypergraph &) = delete;

  StaticHypergraph(StaticHypergraph&& other) :
    _node(other._node),
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_pins(other._num_pins),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _community_support(std::move(other._community_support)) { }

  StaticHypergraph & operator= (StaticHypergraph&& other) {
    _node = other._node;
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _hypernodes = std::move(other._hypernodes);
    _incident_nets = std::move(other._incident_nets);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
    _community_support = std::move(other._community_support);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Number of NUMA hypergraphs
  size_t numNumaHypergraphs() const {
    return 1UL;
  }

  int numaNode() const {
    return _node;
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Initial number of hypernodes on numa node
  HypernodeID initialNumNodes(const int) const {
    return _num_hypernodes;
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _num_removed_hypernodes;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  // ! Initial number of hyperedges on numa node
  HyperedgeID initialNumEdges(const int) const {
    return _num_hyperedges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Initial number of pins on numa node
  HypernodeID initialNumPins(const int) const {
    return _num_pins;
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _total_degree;
  }

  // ! Initial sum of the degree of all vertices on numa node
  HypernodeID initialTotalVertexDegree(const int) const {
    return _total_degree;
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Recomputes the total weight of the hypergraph (parallel)
  void updateTotalWeight(const TaskGroupID) {
    _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 0,
      [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          weight += this->_hypernodes[hn].weight();
        }
        return weight;
      }, std::plus<HypernodeWeight>());
  }

  // ! Recomputes the total weight of the hypergraph (sequential)
  void updateTotalWeight() {
    _total_weight = 0;
    for ( const HypernodeID& hn : nodes() ) {
      _total_weight += nodeWeight(hn);
    }
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID task_group_id, const F& f) {
    static_cast<const StaticHypergraph&>(*this).doParallelForAllNodes(task_group_id, f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID, const F& f) const {
    tbb::parallel_for(0UL, _num_hypernodes, [&](const HypernodeID& id) {
      const HypernodeID hn = common::get_global_vertex_id(_node, id);
      if ( nodeIsEnabled(hn) ) {
        f(hn);
      }
    });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) {
    static_cast<const StaticHypergraph&>(*this).doParallelForAllEdges(task_group_id, f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID, const F& f) const {
    tbb::parallel_for(0UL, _num_hyperedges, [&](const HyperedgeID& id) {
      const HyperedgeID he = common::get_global_edge_id(_node, id);
      if ( edgeIsEnabled(he) ) {
        f(he);
      }
    });
  }

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    const HypernodeID start = common::get_global_vertex_id(_node, 0);
    const HypernodeID end = common::get_global_vertex_id(_node, _num_hypernodes);
    return IteratorRange<HypernodeIterator>(
      HypernodeIterator(_hypernodes.data(), start, end),
      HypernodeIterator(_hypernodes.data() + _num_hypernodes, end, end));
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HypernodeIterator> nodes(const int) const {
    return nodes();
  }

  // ! Returns a range of the active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    const HyperedgeID start = common::get_global_edge_id(_node, 0);
    const HyperedgeID end = common::get_global_edge_id(_node, _num_hyperedges);
    return IteratorRange<HyperedgeIterator>(
      HyperedgeIterator(_hyperedges.data(), start, end),
      HyperedgeIterator(_hyperedges.data() + _num_hyperedges, end, end));
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HyperedgeIterator> edges(const int) const {
    return edges();
  }


  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    const Hypernode& hn = hypernode(u);
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets.cbegin() + hn.firstEntry(),
      _incident_nets.cbegin() + hn.firstInvalidEntry());
  }

  // ! Returns a range to loop over all active multi-pin hyperedges of hypernode u.
  IteratorRange<IncidentNetsIterator> multiPinIncidentEdges(const HypernodeID, const PartitionID) const {
    ERROR("multiPinIncidentEdges(u,c) is not supported in static hypergraph");
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets.cend(), _incident_nets.cend());
  }

  // ! Returns a range to loop over the set of all active incident
  // ! hyperedges of hypernode u that are not single-pin community hyperedges.
  IteratorRange<IncidentNetsIterator> activeIncidentEdges(const HypernodeID, const PartitionID) const {
    ERROR("activeIncidentEdges(u,c) is not supported in static hypergraph");
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets.cend(), _incident_nets.cend());
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const Hyperedge& he = hyperedge(e);
    return IteratorRange<IncidenceIterator>(
      _incidence_array.cbegin() + he.firstEntry(),
      _incidence_array.cbegin() + he.firstInvalidEntry());
  }

  // ! Returns a range to loop over the pins of hyperedge e that belong to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.pins(*this, e, community_id);
  }

  // ! Returns a range to loop over the set of communities contained in hyperedge e.
  IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {
    return _community_support.communities(e);
  }

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! Can be used to map the global vertex ids to a consecutive range
  // ! of nodes between [0,|V|).
  HypernodeID originalNodeID(const HypernodeID u) const {
    return hypernode(u).originalNodeID();
  }

  // ! Reverse operation of originalNodeID(u)
  HypernodeID globalNodeID(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes);
    return u;
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).size();
  }

  // ! Returns, if the corresponding vertex is high degree vertex
  bool isHighDegreeVertex(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).isHighDegreeVertex();
  }
  // ! Marks all vertices with a degree greater the threshold
  // ! as high degree vertex
  void markAllHighDegreeVertices(const TaskGroupID task_group_id,
                                 const HypernodeID high_degree_threshold) {
    doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      if ( nodeDegree(hn) >= high_degree_threshold ) {
        hypernode(hn).markAsHighDegreeVertex();
      }
    });
  }

  // ! Number of invalid incident nets
  HyperedgeID numInvalidIncidentNets(const HypernodeID) const {
    ERROR("numInvalidIncidentNets(u) is not supported in static hypergraph");
    return kInvalidHyperedge;
  }

  // ! Contraction index of the vertex in the contraction hierarchy
  HypernodeID contractionIndex(const HypernodeID) const {
    ERROR("contractionIndex(u) is not supported in static hypergraph");
    return kInvalidHypernode;
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  // ! Disables a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  // ! Removes a hypernode (must be enabled before)
  void removeHypernode(const HypernodeID u) {
    hypernode(u).disable();
    ++_num_removed_hypernodes;
  }

  // ####################### Hyperedge Information #######################

  // ! Returns for a edge of the hypergraph its original edge id
  // ! Can be used to map the global edge ids to a consecutive range
  // ! of edges between [0,|E|).
  HyperedgeID originalEdgeID(const HyperedgeID e) const {
    return hyperedge(e).originalEdgeID();
  }

  // ! Reverse operation of originalEdgeID(e)
  HypernodeID globalEdgeID(const HyperedgeID e) const {
    ASSERT(e < _num_hyperedges);
    return e;
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).weight();
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).setWeight(weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).size();
  }

  // ! Hash value defined over the pins of a hyperedge
  size_t edgeHash(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).hash();
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return !hyperedge(e).isDisabled();
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hyperedge(e).enable();
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
  }

  // ####################### Community Hyperedge Information #######################

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeWeight(e, community_id);
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    _community_support.setEdgeWeight(e, community_id, weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeSize(e, community_id);
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeHash(e, community_id);
  }

  // ####################### Community Information #######################

  // ! Number of communities
  PartitionID numCommunities() const {
    return _community_support.numCommunities();
  }

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityID();
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID u, const PartitionID community_id) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setCommunityID(community_id);
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  HypernodeID communityNodeId(const HypernodeID u) const {
    return _community_support.communityNodeId(u);
  }

  // ! Number of hypernodes in community
  HypernodeID numCommunityHypernodes(const PartitionID community) const {
    return _community_support.numCommunityHypernodes(community);
  }

  // ! Number of pins in community
  HypernodeID numCommunityPins(const PartitionID community) const {
    return _community_support.numCommunityPins(community);
  }

  // ! Total degree of community
  HyperedgeID communityDegree(const PartitionID community) const {
    return _community_support.communityDegree(community);
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    return _community_support.numCommunitiesInHyperedge(e);
  }

  // ! Numa node to which community is assigned to
  PartitionID communityNumaNode(const PartitionID) const {
    return 0;
  }

  // ! Sets the community to numa node mapping
  void setCommunityNodeMapping(std::vector<PartitionID>&&) { }

  // ####################### Contract / Uncontract #######################

  Memento contract(const HypernodeID, const HypernodeID) {
    ERROR("contract(u,v) is not supported in static hypergraph");
    return Memento();
  }

  Memento contract(const HypernodeID, const HypernodeID, const PartitionID) {
    ERROR("contract(u,v,c) is not supported in static hypergraph");
    return Memento();
  }

  std::pair<StaticHypergraph, parallel::scalable_vector<HypernodeID>> contract(
    const parallel::scalable_vector<HypernodeID>&,
    const TaskGroupID) const {
    ERROR("contract(communities,id) is not supported in static hypergraph");
    return std::make_pair(StaticHypergraph(), parallel::scalable_vector<HypernodeID>());
  }

  void uncontract(const Memento&, parallel::scalable_vector<HyperedgeID>&) {
    ERROR("uncontract(memento,parallel_he) is not supported in static hypergraph");
  }

  void uncontract(const std::vector<Memento>&,
                  parallel::scalable_vector<HyperedgeID>&,
                  const kahypar::ds::FastResetFlagArray<>&,
                  const bool) {
    ERROR("uncontract(...) is not supported in static hypergraph");
  }

  void restoreDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("restoreDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in static hypergraph");
  }

  parallel::scalable_vector<HyperedgeID> findDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("findDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in static hypergraph");
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
      removeIncidentEdgeFromHypernode(he, pin);
    }
    disableHyperedge(he);
  }

  void removeSinglePinCommunityEdge(const HyperedgeID, const PartitionID) {
    ERROR("removeSinglePinCommunityEdge(e,c) is not supported in static hypergraph");
  }

  void removeParallelEdge(const HyperedgeID, const PartitionID) {
    ERROR("removeParallelEdge(e,c) is not supported in static hypergraph");
  }

  // ! Restores an hyperedge of a certain size.
  void restoreEdge(const HyperedgeID he, const size_t,
                   const HyperedgeID representative = kInvalidHyperedge) {
    unused(representative);
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "already enabled");
    enableHyperedge(he);
    for ( const HypernodeID& pin : pins(he) ) {
      insertIncidentEdgeToHypernode(he, pin);
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
    ERROR("restoreParallelHyperedge(...) is not supported in static hypergraph");
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
  void initializeCommunities(const TaskGroupID,
                             const parallel::scalable_vector<StaticHypergraph>& hypergraphs = {}) {
    _community_support.initialize(*this, hypergraphs);
  }

  /*!
  * Initializes community hyperedges.
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to a range of consecutive pins with
  *       same community in that hyperedge
  */
  void initializeCommunityHyperedges(const TaskGroupID,
                                     const parallel::scalable_vector<StaticHypergraph>& hypergraphs = {}) {
    _community_support.initializeCommunityHyperedges(*this, hypergraphs);
  }

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   */
  void removeCommunityHyperedges(const TaskGroupID,
                                 const parallel::scalable_vector<HypernodeID>& contraction_index = {},
                                 const parallel::scalable_vector<StaticHypergraph>& hypergraphs = {}) {
    _community_support.removeCommunityHyperedges(contraction_index, hypergraphs);
  }

  void buildContractionHierarchy(const std::vector<Memento>&) {
    ERROR("buildContractionHierarchy(mementos) is not supported in static hypergraph");
  }

  void invalidateDisabledHyperedgesFromIncidentNets(const TaskGroupID) {
    ERROR("invalidateDisabledHyperedgesFromIncidentNets(id) is not supported in static hypergraph");
  }

  // ####################### Copy #######################

  // ! Copy static hypergraph in parallel
  StaticHypergraph copy(const TaskGroupID task_group_id) {
    StaticHypergraph hypergraph;

    hypergraph._node = _node;
    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    tbb::parallel_invoke([&] {
      hypergraph._hypernodes.resize(_hypernodes.size());
      memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
        sizeof(Hypernode) * _hypernodes.size());
    }, [&] {
      hypergraph._incident_nets.resize(_incident_nets.size());
      memcpy(hypergraph._incident_nets.data(), _incident_nets.data(),
        sizeof(HyperedgeID) * _incident_nets.size());
    }, [&] {
      hypergraph._hyperedges.resize(_hyperedges.size());
      memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
        sizeof(Hyperedge) * _hyperedges.size());
    }, [&] {
      hypergraph._incidence_array.resize(_incidence_array.size());
      memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
        sizeof(HypernodeID) * _incidence_array.size());
    }, [&] {
      hypergraph._community_support = _community_support.copy(task_group_id);
    });

    return hypergraph;
  }

  // ! Copy static hypergraph sequential
  StaticHypergraph copy() {
    StaticHypergraph hypergraph;

    hypergraph._node = _node;
    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
      sizeof(Hypernode) * _hypernodes.size());
    hypergraph._incident_nets.resize(_incident_nets.size());
    memcpy(hypergraph._incident_nets.data(), _incident_nets.data(),
      sizeof(HyperedgeID) * _incident_nets.size());

    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
      sizeof(Hyperedge) * _hyperedges.size());
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
      sizeof(HypernodeID) * _incidence_array.size());

    hypergraph._community_support = _community_support.copy();

    return hypergraph;
  }

 private:
  friend class StaticHypergraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;
  template <typename Hypergraph,
            typename HardwareTopology,
            typename TBBNumaArena>
  friend class NumaHypergraph;

  // ####################### Hypernode Information #######################

  // ! Accessor for hypernode-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    const HyperedgeID local_pos = common::get_local_position_of_vertex(u);
    ASSERT(_node == common::get_numa_node_of_vertex(u),
      "Hypernode" << u << "is not part of hypergraph on numa node" << _node);
    ASSERT(local_pos <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[local_pos];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const StaticHypergraph&>(*this).hypernode(u));
  }

  // ####################### Hyperedge Information #######################

  // ! Accessor for hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
    const HyperedgeID local_pos = common::get_local_position_of_edge(e);
    ASSERT(_node == common::get_numa_node_of_edge(e),
      "Hyperedge" << e << "is not part of hypergraph on numa node" << _node);
    ASSERT(local_pos <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[local_pos];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const StaticHypergraph&>(*this).hyperedge(e));
  }

  // ####################### Initialization / Reset Functions #######################

  void finalizeCommunityNodeIds(const parallel::scalable_vector<StaticHypergraph>& hypergraphs) {
    _community_support.finalizeCommunityNodeIds(*this, hypergraphs);
  }

  // ####################### Remove / Restore Hyperedges #######################

  // ! Removes hyperedge e from the incident nets of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID e,
                                                                       const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");

    Hypernode& hn = hypernode(u);
    size_t incident_nets_pos = hn.firstEntry();
    for ( ; incident_nets_pos < hn.firstInvalidEntry(); ++incident_nets_pos ) {
      if ( _incident_nets[incident_nets_pos] == e ) {
        break;
      }
    }
    ASSERT(incident_nets_pos < hn.firstInvalidEntry());
    swap(_incident_nets[incident_nets_pos], _incident_nets[hn.firstInvalidEntry() - 1]);
    hn.setSize(hn.size() - 1);
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID e,
                                                                     const HypernodeID u) {
    using std::swap;
    Hypernode& hn = hypernode(u);
    ASSERT(!hn.isDisabled(), "Hypernode" << u << "is disabled");
    HEAVY_REFINEMENT_ASSERT(std::count(_incident_nets.cbegin() + hn.firstEntry(),
                                       _incident_nets.cend() + hn.firstInvalidEntry(), e) == 0,
                        "HN" << u << "is already connected to HE" << e);
    const size_t incident_nets_start = hn.firstInvalidEntry();
    const size_t incident_nets_end = hypernode(u + 1).firstEntry();
    size_t incident_nets_pos = incident_nets_start;
    for ( ; incident_nets_pos < incident_nets_end; ++incident_nets_pos ) {
      if ( _incident_nets[incident_nets_pos] == e ) {
        break;
      }
    }
    ASSERT(incident_nets_pos < incident_nets_end);
    swap(_incident_nets[incident_nets_start], _incident_nets[incident_nets_pos]);
    hn.setSize(hn.size() + 1);
  }

  // ! NUMA node of hypergraph
  int _node;
  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of removed hypernodes
  HypernodeID _num_removed_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! Hypernodes
  parallel::scalable_vector<Hypernode> _hypernodes;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Hyperedges
  parallel::scalable_vector<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;

  // ! Community Information and Stats
  CommunitySupport<StaticHypergraph> _community_support;
};

} // namespace ds
} // namespace mt_kahypar