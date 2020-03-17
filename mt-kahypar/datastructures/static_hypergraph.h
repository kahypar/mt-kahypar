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

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_scan.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/community_support.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

// Forward
class StaticHypergraphFactory;

class StaticHypergraph {

  static constexpr bool enable_heavy_assert = false;

  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

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
      _valid(false) { }

    // Sentinel Constructor
    Hypernode(const size_t begin) :
      _begin(begin),
      _size(0),
      _original_id(kInvalidHypernode),
      _weight(1),
      _community_id(0),
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

    // Sentinel Constructor
    Hyperedge(const size_t begin) :
      _begin(begin),
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

  ~StaticHypergraph() {
    freeInternalData();
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
    _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
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
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& id) {
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
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& id) {
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

  bool hasCommunityNodeMapping() const {
    return false;
  }

  // ! Numa node to which community is assigned to
  PartitionID communityNumaNode(const PartitionID) const {
    return 0;
  }

  // ! Sets the community to numa node mapping
  void setCommunityNodeMapping(parallel::scalable_vector<PartitionID>&&) { }

  // ! Returns a copy of community to numa node mapping
  parallel::scalable_vector<PartitionID> communityNodeMapping() const {
    return { 0 };
  }

  // ####################### Contract / Uncontract #######################

  Memento contract(const HypernodeID, const HypernodeID) {
    ERROR("contract(u,v) is not supported in static hypergraph");
    return Memento();
  }

  Memento contract(const HypernodeID, const HypernodeID, const PartitionID) {
    ERROR("contract(u,v,c) is not supported in static hypergraph");
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
  std::pair<StaticHypergraph, parallel::scalable_vector<HypernodeID>> contract(
    const parallel::scalable_vector<HypernodeID>& communities,
    const TaskGroupID task_group_id) const {
    ASSERT(communities.size() == _num_hypernodes);

    // #################### STAGE 1 ####################
    // Remapping of vertex ids
    utils::Timer::instance().start_timer("compute_cluster_mapping", "Compute Cluster Mapping");
    parallel::scalable_vector<HypernodeID> mapping(_num_hypernodes, kInvalidHypernode);
    HypernodeID num_hypernodes = 0;
    for ( const HypernodeID& hn : nodes() ) {
      ASSERT(hn < _num_hypernodes);
      HypernodeID community = communities[hn];
      if ( mapping[community] == kInvalidHypernode ) {
        // Setup mapping from community id to a vertex id
        // in the contracted hypergraph
        mapping[community] = num_hypernodes++;
      }
    }

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
      ASSERT(hn < communities.size());
      return mapping[communities[hn]];
    };

    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
    parallel::scalable_vector<PartitionID> community_ids;
    tbb::parallel_invoke([&] {
      hn_weights.assign(num_hypernodes, parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
    }, [&] {
      community_ids.assign(num_hypernodes, 0);
    });

    doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
      const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
      ASSERT(coarse_hn < num_hypernodes);
      // Weight vector is atomic => thread-safe
      hn_weights[coarse_hn] += nodeWeight(hn);
      // In case community detection is enabled all vertices matched to one vertex
      // in the contracted hypergraph belong to same community. Otherwise, all communities
      // are default assigned to community 0
      community_ids[coarse_hn] = communityID(hn);
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");


    // #################### STAGE 2 ####################
    // We iterate over all hyperedges in parallel and remap their ids
    // to the ones determined in the step before. Furthermore, duplicates
    // and disabled hyperedges are removed. The hyperedges are then inserted
    // into a streaming map with their hash as key. All hyperedges with the same
    // hash are then present in the same bucket of the streaming map, which
    // makes it possible to detect parallel hyperedges in parallel.
    utils::Timer::instance().start_timer("contracting_hyperedges", "Contracting Hyperedges");
    StreamingMap<size_t, ContractedHyperedge> hash_to_hyperedge;
    contractHyperedges(task_group_id, hash_to_hyperedge, map_to_coarse_hypergraph);

    using HyperedgeMap = parallel::scalable_vector<parallel::scalable_vector<ContractedHyperedge>>;
    HyperedgeMap hyperedge_buckets(hash_to_hyperedge.size());
    hash_to_hyperedge.copy(hyperedge_buckets, [&](const size_t key) {
      return key % hash_to_hyperedge.size();
    });
    utils::Timer::instance().stop_timer("contracting_hyperedges");

    // #################### STAGE 3 ####################
    // We iterate now in parallel over each bucket and sort each bucket
    // after its hash. A bucket is processed by one thread and parallel
    // hyperedges are detected by comparing the pins of hyperedges with
    // the same hash.

    utils::Timer::instance().start_timer("remove_parallel_hyperedges", "Remove Parallel Hyperedges");
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
    // Stores the prefix sum over the number of pins in each bucket
    parallel::scalable_vector<HyperedgeID> num_pins_prefix_sum(hyperedge_buckets.size() + 1, 0);
    // For each node we aggregate the number of incident nets
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> num_incident_nets(
      num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0));

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

          // Determine position for each hyperedge and its pin
          // in the hyperedge vector and incidence array
          contracted_he_lhs.he_idx = num_hyperedges_prefix_sum[bucket + 1]++;
          contracted_he_lhs.pin_idx = num_pins_prefix_sum[bucket + 1];
          num_pins_prefix_sum[bucket + 1] += contracted_he_lhs.hyperedge.size();
          // Aggregate the number of incident nets of each vertex
          for ( const HypernodeID& pin : contracted_he_lhs.hyperedge ) {
            ASSERT(pin < num_hypernodes);
            ++num_incident_nets[pin];
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
      num_pins_prefix_sum[i] += num_pins_prefix_sum[i - 1];
    }
    utils::Timer::instance().stop_timer("remove_parallel_hyperedges");

    // #################### STAGE 4 ####################
    // Initialize hypergraph
    utils::Timer::instance().start_timer("construct_contracted_hypergraph", "Construct Contracted Hypergraph");
    StaticHypergraph hypergraph;
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges_prefix_sum.back();
    hypergraph._num_pins = num_pins_prefix_sum.back();
    hypergraph._total_degree = num_pins_prefix_sum.back();

    // Compute start position of incident nets for each vertex
    // in incident net array
    utils::Timer::instance().start_timer("incident_net_prefix_sum", "Incident Net Prefix Sum");
    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>>
      incident_nets_prefix_sum(num_incident_nets);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, UI64(num_hypernodes)), incident_nets_prefix_sum);
    utils::Timer::instance().stop_timer("incident_net_prefix_sum");

    tbb::parallel_invoke([&] {
      // Setup hypernodes
      utils::Timer::instance().start_timer("setup_hypernodes", "Setup Hypernodes", true);
      hypergraph._hypernodes.resize(hypergraph._num_hypernodes);
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID hn) {
        const size_t incident_nets_pos = incident_nets_prefix_sum[hn];
        const size_t incident_nets_size = hn == 0 ? incident_nets_prefix_sum[hn + 1] :
          incident_nets_prefix_sum[hn + 1] - incident_nets_prefix_sum[hn];
        ASSERT(incident_nets_pos + incident_nets_size <=
              hypergraph._total_degree);
        Hypernode& hn_obj = hypergraph._hypernodes[hn];
        hn_obj.enable();
        hn_obj.setFirstEntry(incident_nets_pos);
        hn_obj.setSize(incident_nets_size);
        hn_obj.setOriginalNodeID(hn);
        hn_obj.setWeight(hn_weights[hn]);
        hn_obj.setCommunityID(community_ids[hn]);
      });
      utils::Timer::instance().stop_timer("setup_hypernodes");
    }, [&] {
      utils::Timer::instance().start_timer("setup_hyperedges", "Setup Hyperedges", true);
      // Setup hyperedges, incidence and incident nets array
      hypergraph._hyperedges.resize(hypergraph._num_hyperedges);
      hypergraph._incidence_array.resize(hypergraph._num_pins);
      hypergraph._incident_nets.resize(hypergraph._total_degree);

      parallel::scalable_vector<parallel::IntegralAtomicWrapper<HyperedgeID>> incident_nets_pos(
        num_hypernodes, parallel::IntegralAtomicWrapper<HyperedgeID>(0));
      tbb::parallel_for(0UL, hyperedge_buckets.size(), [&](const size_t bucket) {
        for ( const ContractedHyperedge& contracted_he : hyperedge_buckets[bucket] ) {
          if ( !contracted_he.is_parallel ) {
            const HyperedgeID he = num_hyperedges_prefix_sum[bucket] + contracted_he.he_idx;
            const size_t incidence_array_pos = num_pins_prefix_sum[bucket] + contracted_he.pin_idx;
            ASSERT(he < hypergraph._hyperedges.size());
            ASSERT(incidence_array_pos + contracted_he.hyperedge.size() <=
                  hypergraph._incidence_array.size());

            // Setup hyperedge
            Hyperedge& he_obj = hypergraph._hyperedges[he];
            he_obj.enable();
            he_obj.setFirstEntry(incidence_array_pos);
            he_obj.setSize(contracted_he.hyperedge.size());
            he_obj.setOriginalEdgeID(he);
            he_obj.setWeight(contracted_he.weight);
            he_obj.hash() = contracted_he.hash;

            // Copy content of hyperedge to incidence array
            memcpy(hypergraph._incidence_array.data() + incidence_array_pos,
                  contracted_he.hyperedge.data(), sizeof(HypernodeID) * he_obj.size());

            // Add hyperedge as incident net to all pins
            for ( const HypernodeID& pin : hypergraph.pins(he) ) {
              const size_t incident_nets_position =
                incident_nets_prefix_sum[pin] + incident_nets_pos[pin]++;
              ASSERT(incident_nets_position < hypergraph._incident_nets.size());
              hypergraph._incident_nets[incident_nets_position] = he;
            }
          }
        }
        // We free memory here in parallel, because this can become a major
        // bottleneck, if memory is freed sequential after function return
        parallel::scalable_vector<ContractedHyperedge> tmp_bucket;
        hyperedge_buckets[bucket] = std::move(tmp_bucket);
      });
      utils::Timer::instance().stop_timer("setup_hyperedges");
    });
    utils::Timer::instance().stop_timer("construct_contracted_hypergraph");

    // Initialize Communities and Update Total Weight
    utils::Timer::instance().start_timer("setup_communities", "Setup Communities");
    tbb::parallel_invoke([&] {
      if ( _community_support.isInitialized() ) {
        hypergraph.initializeCommunities(task_group_id);
        if ( _community_support.areCommunityHyperedgesInitialized() ) {
          hypergraph.initializeCommunityHyperedges(task_group_id);
        }
      }
    }, [&] {
      hypergraph.updateTotalWeight(task_group_id);
    });
    utils::Timer::instance().stop_timer("setup_communities");

    utils::Timer::instance().start_timer("free_internal_data", "Free Internal Data");
    // We free memory here in parallel, because this can become a major
    // bottleneck, if memory is freed sequential after function return
    parallel::parallel_free(hn_weights, community_ids, num_incident_nets);
    utils::Timer::instance().stop_timer("free_internal_data");

    return std::make_pair(std::move(hypergraph), std::move(mapping));
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
  void initializeCommunities(const TaskGroupID task_group_id,
                             const parallel::scalable_vector<StaticHypergraph>& hypergraphs = {}) {
    _community_support.initialize(*this, hypergraphs, task_group_id);
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

  // Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      tbb::parallel_invoke([&] {
        _community_support.freeInternalData();
      }, [&] {
        parallel::parallel_free(
          _hypernodes, _incident_nets,
          _hyperedges, _incidence_array);
      });
    }
    _num_hypernodes = 0;
    _num_hyperedges = 0;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
    parent->addChild("Incident Nets", sizeof(HyperedgeID) * _incident_nets.size());
    parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());

    utils::MemoryTreeNode* community_support_node = parent->addChild("Community Support");
    _community_support.memoryConsumption(community_support_node);
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

  // ####################### Contract / Uncontract #######################

  template<typename F>
  void contractHyperedges(const TaskGroupID task_group_id,
                          StreamingMap<size_t, ContractedHyperedge>& hash_to_hyperedge,
                          const F& mapping_to_coarse_hypergraph) const {
    doParallelForAllEdges(task_group_id, [&](const HyperedgeID& he) {
      parallel::scalable_vector<HypernodeID> hyperedge;
      hyperedge.reserve(edgeSize(he));
      for ( const HypernodeID pin : pins(he) ) {
        hyperedge.emplace_back(mapping_to_coarse_hypergraph(pin));
      }

      // Removing duplicates
      std::sort(hyperedge.begin(), hyperedge.end());
      hyperedge.erase(std::unique(hyperedge.begin(), hyperedge.end()), hyperedge.end());

      // Removing disabled hypernodes
      while ( !hyperedge.empty() && hyperedge.back() == kInvalidHypernode ) {
        hyperedge.pop_back();
      }

      // Remove single-pin hyperedges
      if ( hyperedge.size() > 1 ) {
        // Compute hash of hyperedge
        size_t he_hash = kEdgeHashSeed;
        for ( const HypernodeID& pin : hyperedge ) {
          he_hash += kahypar::math::hash(pin);
        }
        hash_to_hyperedge.stream(he_hash,
          ContractedHyperedge { he, he_hash, edgeWeight(he),
            false, _node, std::move(hyperedge), ID(0), ID(0) } );
      }
    });
  }

  // ####################### Initialization / Reset Functions #######################

  void finalizeCommunityNodeIds(const parallel::scalable_vector<StaticHypergraph>& hypergraphs,
                                const TaskGroupID task_group_id) {
    _community_support.finalizeCommunityNodeIds(*this, hypergraphs, task_group_id);
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
                                       _incident_nets.cbegin() + hn.firstInvalidEntry(), e) == 0,
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