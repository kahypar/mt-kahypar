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
#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/community_support.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/vector.h"
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

    Hypernode(const bool valid) :
      _begin(0),
      _size(0),
      _original_id(kInvalidHypernode),
      _weight(1),
      _community_id(0),
      _valid(valid) { }

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

  using IncidenceArray = Vector<HypernodeID>;
  using IncidentNets = Vector<HyperedgeID>;

  // ! Contains data structures that are needed during multilevel contractions.
  // ! Struct is allocated on top level hypergraph and passed to each contracted
  // ! hypergraph such that memory can be reused in consecutive contractions.
  struct TmpContractionBuffer {
    explicit TmpContractionBuffer(const HypernodeID num_hypernodes,
                                  const HyperedgeID num_hyperedges,
                                  const HyperedgeID num_pins,
                                  const bool is_numa_buffer = false) :
        is_initialized(false) {
      tbb::parallel_invoke([&] {
        tmp_hypernodes.resize(num_hypernodes);
      }, [&] {
        if ( !is_numa_buffer ) {
          tmp_incident_nets.resize(num_pins);
        }
      }, [&] {
        tmp_num_incident_nets.assign(num_hypernodes,
          parallel::IntegralAtomicWrapper<size_t>(0));
      }, [&] {
        hn_weights.assign(num_hypernodes,
          parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
      }, [&] {
        tmp_hyperedges.resize(num_hyperedges);
      }, [&] {
        tmp_incidence_array.resize(num_pins);
      }, [&] {
        valid_hyperedges.resize(num_hyperedges);
      });
      is_initialized = true;
    }

    void freeInternalData() {
      if ( is_initialized ) {
        parallel::parallel_free(tmp_hypernodes, tmp_num_incident_nets,
          hn_weights, tmp_hyperedges, valid_hyperedges);
        is_initialized = false;
      }
    }

    bool is_initialized;
    parallel::scalable_vector<Hypernode> tmp_hypernodes;
    IncidentNets tmp_incident_nets;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> tmp_num_incident_nets;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
    parallel::scalable_vector<Hyperedge> tmp_hyperedges;
    IncidenceArray tmp_incidence_array;
    parallel::scalable_vector<size_t> valid_hyperedges;
  };

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
    _num_removed_hyperedges(0),
    _num_pins(0),
    _total_degree(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
    _community_support(),
    _tmp_contraction_buffer(nullptr),
    _is_root_allocator(false) { }

  StaticHypergraph(const StaticHypergraph&) = delete;
  StaticHypergraph & operator= (const StaticHypergraph &) = delete;

  StaticHypergraph(StaticHypergraph&& other) :
    _node(other._node),
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_removed_hyperedges(other._num_removed_hyperedges),
    _num_pins(other._num_pins),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _community_support(std::move(other._community_support)),
    _tmp_contraction_buffer(std::move(other._tmp_contraction_buffer)),
    _is_root_allocator(other._is_root_allocator) {
    other._is_root_allocator = false;
  }

  StaticHypergraph & operator= (StaticHypergraph&& other) {
    _node = other._node;
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _hypernodes = std::move(other._hypernodes);
    _incident_nets = std::move(other._incident_nets);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
    _community_support = std::move(other._community_support);
    _tmp_contraction_buffer = std::move(other._tmp_contraction_buffer);
    _is_root_allocator = other._is_root_allocator;
    other._is_root_allocator = false;
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
  StaticHypergraph contract(parallel::scalable_vector<HypernodeID>& communities,
                            const TaskGroupID task_group_id) const {
    ASSERT(_tmp_contraction_buffer);
    ASSERT(communities.size() == _num_hypernodes);

    // AUXILLIARY BUFFERS - Reused during multilevel hierarchy to prevent expensive allocations
    parallel::scalable_vector<Hypernode>& tmp_hypernodes = _tmp_contraction_buffer->tmp_hypernodes;
    IncidentNets& tmp_incident_nets = _tmp_contraction_buffer->tmp_incident_nets;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>& tmp_num_incident_nets =
      _tmp_contraction_buffer->tmp_num_incident_nets;
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeWeight>>& hn_weights =
      _tmp_contraction_buffer->hn_weights;
    parallel::scalable_vector<Hyperedge>& tmp_hyperedges = _tmp_contraction_buffer->tmp_hyperedges;
    IncidenceArray& tmp_incidence_array = _tmp_contraction_buffer->tmp_incidence_array;
    parallel::scalable_vector<size_t>& valid_hyperedges = _tmp_contraction_buffer->valid_hyperedges;

    ASSERT(static_cast<size_t>(_num_hypernodes) <= tmp_hypernodes.size());
    ASSERT(static_cast<size_t>(_total_degree) <= tmp_incident_nets.size());
    ASSERT(static_cast<size_t>(_num_hypernodes) <= tmp_num_incident_nets.size());
    ASSERT(static_cast<size_t>(_num_hypernodes) <= hn_weights.size());
    ASSERT(static_cast<size_t>(_num_hyperedges) <= tmp_hyperedges.size());
    ASSERT(static_cast<size_t>(_num_pins) <= tmp_incidence_array.size());
    ASSERT(static_cast<size_t>(_num_hyperedges) <= valid_hyperedges.size());

    // #################### STAGE 1 ####################
    // Compute vertex ids of coarse hypergraph with a parallel prefix sum
    utils::Timer::instance().start_timer("preprocess_contractions", "Preprocess Contractions");
    parallel::scalable_vector<size_t> mapping(_num_hypernodes, 0);
    doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
      ASSERT(static_cast<size_t>(communities[originalNodeID(hn)]) < mapping.size());
      mapping[communities[hn]] = 1UL;
    });

    // Prefix sum determines vertex ids in coarse hypergraph
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, mapping.size()), mapping_prefix_sum);
    HypernodeID num_hypernodes = mapping_prefix_sum.total_sum();

    // Remap community ids
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
      if ( nodeIsEnabled(hn) ) {
        communities[hn] = mapping_prefix_sum[communities[hn]];
      } else {
        communities[hn] = kInvalidHypernode;
      }
    });

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
      ASSERT(hn < communities.size());
      return communities[hn];
    };

    tbb::parallel_invoke([&] {
      hn_weights.assign(num_hypernodes, parallel::IntegralAtomicWrapper<HypernodeWeight>(0));
    }, [&] {
      tmp_hypernodes.assign(num_hypernodes, Hypernode(true));
    }, [&] {
      tmp_num_incident_nets.assign(num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0));
    });

    doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
      const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
      ASSERT(coarse_hn < num_hypernodes, V(coarse_hn) << V(num_hypernodes));
      // Weight vector is atomic => thread-safe
      hn_weights[coarse_hn] += nodeWeight(hn);
      // In case community detection is enabled all vertices matched to one vertex
      // in the contracted hypergraph belong to same community. Otherwise, all communities
      // are default assigned to community 0
      tmp_hypernodes[coarse_hn].setCommunityID(communityID(hn));
      // Aggregate upper bound for number of incident nets of the contracted vertex
      tmp_num_incident_nets[coarse_hn] += nodeDegree(hn);
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
    ConcurrentBucketMap<HyperedgeHash> hyperedge_hash_map;
    hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);
    tbb::parallel_invoke([&] {
      // Contract Hyperedges
      utils::Timer::instance().start_timer("contract_hyperedges", "Contract Hyperedges", true);
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& id) {
        const HyperedgeID he = globalEdgeID(id);
        if ( edgeIsEnabled(he) ) {
          // Copy hyperedge and pins to temporary buffer
          const Hyperedge& e = _hyperedges[id];
          ASSERT(static_cast<size_t>(id) < tmp_hyperedges.size());
          ASSERT(e.firstInvalidEntry() <= tmp_incidence_array.size());
          tmp_hyperedges[id] = e;
          valid_hyperedges[id] = 1;

          // Map pins to vertex ids in coarse graph
          const size_t incidence_array_start = tmp_hyperedges[id].firstEntry();
          const size_t incidence_array_end = tmp_hyperedges[id].firstInvalidEntry();
          for ( size_t pos = incidence_array_start; pos < incidence_array_end; ++pos ) {
            const HypernodeID pin = _incidence_array[pos];
            ASSERT(pos < tmp_incidence_array.size());
            tmp_incidence_array[pos] = map_to_coarse_hypergraph(pin);
          }

          // Remove duplicates and disabled vertices
          auto first_entry_it = tmp_incidence_array.begin() + incidence_array_start;
          std::sort(first_entry_it, tmp_incidence_array.begin() + incidence_array_end);
          auto first_invalid_entry_it = std::unique(first_entry_it, tmp_incidence_array.begin() + incidence_array_end);
          while ( first_entry_it != first_invalid_entry_it && *(first_invalid_entry_it - 1) == kInvalidHypernode ) {
            --first_invalid_entry_it;
          }

          // Update size of hyperedge in temporary hyperedge buffer
          const size_t contracted_size = std::distance(
            tmp_incidence_array.begin() + incidence_array_start, first_invalid_entry_it);
          tmp_hyperedges[id].setSize(contracted_size);


          if ( contracted_size > 1 ) {
            // Compute hash of contracted hyperedge
            size_t he_hash = kEdgeHashSeed;
            for ( size_t pos = incidence_array_start; pos < incidence_array_start + contracted_size; ++pos ) {
              he_hash += kahypar::math::hash(tmp_incidence_array[pos]);
            }
            hyperedge_hash_map.insert(he_hash,
              HyperedgeHash { id, contracted_size, true });
          } else {
            // Hyperedge becomes a single-pin hyperedge
            valid_hyperedges[id] = 0;
            tmp_hyperedges[id].disable();
          }
        } else {
          valid_hyperedges[id] = 0;
        }
      });
      utils::Timer::instance().stop_timer("contract_hyperedges");
    }, [&] {
      // Contract Incident Nets
      utils::Timer::instance().start_timer("tmp_contract_incident_nets", "Tmp Contract Incident Nets", true);

      // Compute start position the incident nets of a coarse vertex in the
      // temporary incident nets array with a parallel prefix sum
      parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>> tmp_incident_nets_pos;
      parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>>
        tmp_incident_nets_prefix_sum(tmp_num_incident_nets);
      tbb::parallel_invoke([&] {
        tbb::parallel_scan(tbb::blocked_range<size_t>(
          0UL, UI64(num_hypernodes)), tmp_incident_nets_prefix_sum);
      }, [&] {
        tmp_incident_nets_pos.assign(num_hypernodes, parallel::IntegralAtomicWrapper<size_t>(0));
      });

      // Write the incident nets of each contracted vertex to the temporary incident net array
      doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
        const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
        const HyperedgeID node_degree = nodeDegree(hn);
        size_t incident_nets_pos = tmp_incident_nets_prefix_sum[coarse_hn] +
          tmp_incident_nets_pos[coarse_hn].fetch_add(node_degree);
        ASSERT(incident_nets_pos + node_degree <= tmp_incident_nets_prefix_sum[coarse_hn + 1]);
        memcpy(tmp_incident_nets.data() + incident_nets_pos,
               _incident_nets.data() + _hypernodes[originalNodeID(hn)].firstEntry(),
               sizeof(HyperedgeID) * node_degree);
      });

      // Setup temporary hypernodes
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID& coarse_hn) {
        // Remove duplicates
        const size_t incident_nets_start = tmp_incident_nets_prefix_sum[coarse_hn];
        const size_t incident_nets_end = tmp_incident_nets_prefix_sum[coarse_hn + 1];
        std::sort(tmp_incident_nets.begin() + incident_nets_start,
                  tmp_incident_nets.begin() + incident_nets_end);
        auto first_invalid_entry_it = std::unique(tmp_incident_nets.begin() + incident_nets_start,
                                                  tmp_incident_nets.begin() + incident_nets_end);

        // Setup pointers to temporary incident nets
        const size_t contracted_size = std::distance(tmp_incident_nets.begin() + incident_nets_start,
                                                     first_invalid_entry_it);
        tmp_hypernodes[coarse_hn].setWeight(hn_weights[coarse_hn]);
        tmp_hypernodes[coarse_hn].setFirstEntry(incident_nets_start);
        tmp_hypernodes[coarse_hn].setSize(contracted_size);
        tmp_hypernodes[coarse_hn].setOriginalNodeID(coarse_hn);
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
      const Hyperedge& lhs_he = tmp_hyperedges[lhs];
      const Hyperedge& rhs_he = tmp_hyperedges[rhs];
      if ( lhs_he.size() == rhs_he.size() ) {
        const size_t lhs_start = lhs_he.firstEntry();
        const size_t rhs_start = rhs_he.firstEntry();
        for ( size_t i = 0; i < lhs_he.size(); ++i ) {
          const size_t lhs_pos = lhs_start + i;
          const size_t rhs_pos = rhs_start + i;
          if ( tmp_incidence_array[lhs_pos] != tmp_incidence_array[rhs_pos] ) {
            return false;
          }
        }
        return true;
      } else {
        return false;
      }
    };

    using BucketElement = typename ConcurrentBucketMap<HyperedgeHash>::Element;
    tbb::parallel_for(0UL, hyperedge_hash_map.numBuckets(), [&](const size_t bucket) {
      auto& hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
      std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
        [&](const BucketElement& lhs, const BucketElement& rhs) {
          return lhs.key < rhs.key || (lhs.key == rhs.key && lhs.value.size < rhs.value.size);
        });

      // Parallel Hyperedge Detection
      for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
        const size_t hash_lhs = hyperedge_bucket[i].key;
        HyperedgeHash& contracted_he_lhs = hyperedge_bucket[i].value;
        if ( contracted_he_lhs.valid ) {
          const HyperedgeID lhs_he = contracted_he_lhs.he;
          HyperedgeWeight lhs_weight = tmp_hyperedges[lhs_he].weight();
          for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
            const size_t hash_rhs = hyperedge_bucket[j].key;
            HyperedgeHash& contracted_he_rhs = hyperedge_bucket[j].value;
            const HyperedgeID rhs_he = contracted_he_rhs.he;
            if ( contracted_he_rhs.valid && hash_lhs == hash_rhs &&
                 check_if_hyperedges_are_parallel(lhs_he, rhs_he) ) {
                // Hyperedges are parallel
                lhs_weight += tmp_hyperedges[rhs_he].weight();
                contracted_he_rhs.valid = false;
                valid_hyperedges[rhs_he] = false;
            } else if ( hash_lhs != hash_rhs ) {
              // In case, hash of both are not equal we go to the next hyperedge
              // because we compared it with all hyperedges that had an equal hash
              break;
            }
          }
          tmp_hyperedges[lhs_he].setWeight(lhs_weight);
        }
      }
      hyperedge_hash_map.clear(bucket);
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
    StaticHypergraph hypergraph;

    // Compute number of hyperedges in coarse graph (those flagged as valid)
    parallel::TBBPrefixSum<size_t> he_mapping(valid_hyperedges);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, UI64(_num_hyperedges)), he_mapping);
    }, [&] {
      hypergraph._hypernodes.resize(num_hypernodes);
    });

    const HyperedgeID num_hyperedges = he_mapping.total_sum();
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges;

    tbb::parallel_invoke([&] {
      utils::Timer::instance().start_timer("setup_hyperedges", "Setup Hyperedges", true);
      utils::Timer::instance().start_timer("compute_he_pointer", "Compute HE Pointer", true);
      // Compute start position of each hyperedge in incidence array
      parallel::scalable_vector<size_t> he_sizes;
      parallel::TBBPrefixSum<size_t> num_pins_prefix_sum(he_sizes);
      tbb::parallel_invoke([&] {
        he_sizes.assign(_num_hyperedges, 0);
        tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& id) {
          if ( he_mapping.value(id) ) {
            he_sizes[id] = tmp_hyperedges[id].size();
          }
        });

        tbb::parallel_scan(tbb::blocked_range<size_t>(
          0UL, UI64(_num_hyperedges)), num_pins_prefix_sum);
        const size_t num_pins = num_pins_prefix_sum.total_sum();
        hypergraph._num_pins = num_pins;
        hypergraph._incidence_array.resize(num_pins);
      }, [&] {
        hypergraph._hyperedges.resize(num_hyperedges);
      });
      utils::Timer::instance().stop_timer("compute_he_pointer");

      utils::Timer::instance().start_timer("setup_incidence_array", "Setup Incidence Array", true);
      // Write hyperedges from temporary buffers to incidence array
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& id) {
        if ( he_mapping.value(id) /* hyperedge is valid */ ) {
          const size_t he_pos = he_mapping[id];
          const size_t incidence_array_start = num_pins_prefix_sum[id];
          Hyperedge& he = hypergraph._hyperedges[he_pos];
          he = std::move(tmp_hyperedges[id]);
          const size_t tmp_incidence_array_start = he.firstEntry();
          std::memcpy(hypergraph._incidence_array.data() + incidence_array_start,
                      tmp_incidence_array.data() + tmp_incidence_array_start,
                      sizeof(HypernodeID) * he.size());
          he.setFirstEntry(incidence_array_start);
          he.setOriginalEdgeID(he_pos);
        }
      });
      utils::Timer::instance().stop_timer("setup_incidence_array");
      utils::Timer::instance().stop_timer("setup_hyperedges");
    }, [&] {
      utils::Timer::instance().start_timer("setup_hypernodes", "Setup Hypernodes", true);
      utils::Timer::instance().start_timer("compute_num_incident_nets", "Compute Num Incident Nets", true);
      // Remap hyperedge ids in temporary incident nets to hyperedge ids of the
      // coarse hypergraph and remove singple-pin/parallel hyperedges.
      parallel::scalable_vector<size_t> incident_nets_sizes(num_hypernodes, 0);
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID& id) {
        const size_t incident_nets_start =  tmp_hypernodes[id].firstEntry();
        size_t incident_nets_end = tmp_hypernodes[id].firstInvalidEntry();
        for ( size_t pos = incident_nets_start; pos < incident_nets_end; ++pos ) {
          const HyperedgeID he = tmp_incident_nets[pos];
          if ( he_mapping.value(he) ) {
            tmp_incident_nets[pos] = he_mapping[he];
          } else {
            std::swap(tmp_incident_nets[pos--], tmp_incident_nets[--incident_nets_end]);
          }
        }
        const size_t incident_nets_size = incident_nets_end - incident_nets_start;
        tmp_hypernodes[id].setSize(incident_nets_size);
        incident_nets_sizes[id] = incident_nets_size;
      });

      // Compute start position of the incident nets for each vertex inside
      // the coarsened incident net array
      parallel::TBBPrefixSum<size_t> num_incident_nets_prefix_sum(incident_nets_sizes);
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, UI64(num_hypernodes)), num_incident_nets_prefix_sum);
      const size_t total_degree = num_incident_nets_prefix_sum.total_sum();
      hypergraph._total_degree = total_degree;
      hypergraph._incident_nets.resize(total_degree);
      utils::Timer::instance().stop_timer("compute_num_incident_nets");

      utils::Timer::instance().start_timer("setup_incident_nets", "Setup Incidenct Nets", true);
      // Write incident nets from temporary buffer to incident nets array
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID& id) {
        const size_t incident_nets_start = num_incident_nets_prefix_sum[id];
        Hypernode& hn = hypergraph._hypernodes[id];
        hn = std::move(tmp_hypernodes[id]);
        const size_t tmp_incident_nets_start = hn.firstEntry();
        std::memcpy(hypergraph._incident_nets.data() + incident_nets_start,
                    tmp_incident_nets.data() + tmp_incident_nets_start,
                    sizeof(HyperedgeID) * hn.size());
        hn.setFirstEntry(incident_nets_start);
      });
      utils::Timer::instance().stop_timer("setup_incident_nets");
      utils::Timer::instance().stop_timer("setup_hypernodes");
    });
    utils::Timer::instance().stop_timer("contract_hypergraph");


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

    hypergraph._tmp_contraction_buffer = _tmp_contraction_buffer;
    return hypergraph;
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
    ++_num_removed_hyperedges;
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
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
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

    hypergraph.allocateTmpContractionBuffer();
    return hypergraph;
  }

  // ! Copy static hypergraph sequential
  StaticHypergraph copy() {
    StaticHypergraph hypergraph;

    hypergraph._node = _node;
    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
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

  // ! Allocate the temporary contraction buffer
  void allocateTmpContractionBuffer(const bool is_numa_buffer = false) {
    if ( !_is_root_allocator ) {
      _tmp_contraction_buffer = new TmpContractionBuffer(
        _num_hypernodes, _num_hyperedges, _num_pins, is_numa_buffer);
      _is_root_allocator = true;
    }
  }

  // ! Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      tbb::parallel_invoke([&] {
        _community_support.freeInternalData();
      }, [&] {
        if ( _tmp_contraction_buffer && _is_root_allocator ) {
          _tmp_contraction_buffer->freeInternalData();
          free(_tmp_contraction_buffer);
          _tmp_contraction_buffer = nullptr;
        }
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
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! Hypernodes
  Vector<Hypernode> _hypernodes;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Hyperedges
  Vector<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;

  // ! Community Information and Stats
  CommunitySupport<StaticHypergraph> _community_support;

  // ! Data that is reused throughout the multilevel hierarchy
  // ! to contract the hypergraph and to prevent expensive allocations
  TmpContractionBuffer* _tmp_contraction_buffer;
  // ! If true, than the contraction buffer was allocated within this hypergraph
  bool _is_root_allocator;
};

} // namespace ds
} // namespace mt_kahypar