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

#include <mutex>
#include <queue>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"
#include "tbb/concurrent_queue.h"

#include "kahypar/meta/mandatory.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/community_support.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/incident_net_vector.h"
#include "mt-kahypar/datastructures/contraction_tree.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"

namespace mt_kahypar {
namespace ds {

// Forward
class DynamicHypergraphFactory;
template <typename Hypergraph,
          typename HypergraphFactory>
class PartitionedHypergraph;

class DynamicHypergraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  // ! In order to update gain cache correctly for an uncontraction (u,v),
  // ! the partitioned hypergraph has to know wheter v replaces u in a hyperedge
  // ! or both a incident to that hyperedge after uncontraction. Therefore, the partitioned
  // ! hypergraph passes two lambda functions to the batch uncontraction function, one for
  // ! each case.
  using UncontractionFunction = std::function<void (const HypernodeID, const HypernodeID, const HyperedgeID)>;
  #define NOOP_BATCH_FUNC [] (const HypernodeID, const HypernodeID, const HyperedgeID) { }

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HypernodeID;

    Hypernode() :
      _weight(1),
      _community_id(0),
      _batch_idx(std::numeric_limits<HypernodeID>::max()),
      _valid(false) { }

    Hypernode(const bool valid) :
      _weight(1),
      _community_id(0),
      _batch_idx(std::numeric_limits<HypernodeID>::max()),
      _valid(valid) { }

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

    HyperedgeWeight weight() const {
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    PartitionID communityID() const {
      return _community_id;
    }

    void setCommunityID(const PartitionID community_id) {
      ASSERT(!isDisabled());
      _community_id = community_id;
    }

    HypernodeID batchIndex() const {
      return _batch_idx;
    }

    void setBatchIndex(const HypernodeID batch_idx) {
      _batch_idx = batch_idx;
    }

   private:
    // ! Hypernode weight
    HyperedgeWeight _weight;
    // ! Community id
    PartitionID _community_id;
    // ! Index of the uncontraction batch in which this hypernode is contained in
    HypernodeID _batch_idx;
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
      _weight(1),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    // Sentinel Constructor
    Hyperedge(const size_t begin) :
      _begin(begin),
      _size(0),
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

    void incrementSize() {
      ASSERT(!isDisabled());
      ++_size;
    }

    void decrementSize() {
      ASSERT(!isDisabled());
      --_size;
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

  enum class ContractionResult : uint8_t {
    CONTRACTED = 0,
    PENDING_CONTRACTIONS = 1,
    WEIGHT_LIMIT_REACHED = 2
  };

  using ContractionInterval = typename ContractionTree::Interval;
  using ChildIterator = typename ContractionTree::ChildIterator;

  struct PQBatchUncontractionElement {
    size_t _subtree_size;
    std::pair<ChildIterator, ChildIterator> _iterator;
  };

  struct PQElementComparator {
    bool operator()(const PQBatchUncontractionElement& lhs, const PQBatchUncontractionElement& rhs){
        return lhs._subtree_size < rhs._subtree_size;
    }
  };

  using IncidenceArray = Array<HypernodeID>;
  using IncidentNets = parallel::scalable_vector<IncidentNetVector<HyperedgeID>>;
  using OwnershipVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<bool>>;
  using ThreadLocalHyperedgeVector = tbb::enumerable_thread_specific<parallel::scalable_vector<HyperedgeID>>;
  using ThreadLocalBitset = tbb::enumerable_thread_specific<parallel::scalable_vector<bool>>;

 public:
  static constexpr bool is_static_hypergraph = false;
  static constexpr bool is_partitioned = false;
  static constexpr size_t SIZE_OF_HYPERNODE = sizeof(Hypernode);
  static constexpr size_t SIZE_OF_HYPEREDGE = sizeof(Hyperedge);

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename IncidenceArray::const_iterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename IncidentNetVector<HyperedgeID>::const_iterator;
  // ! Iterator to iterate over the set of communities contained in a hyperedge
  using CommunityIterator = typename CommunitySupport<StaticHypergraph>::CommunityIterator;


  explicit DynamicHypergraph() :
    _num_hypernodes(0),
    _num_removed_hypernodes(0),
    _num_hyperedges(0),
    _num_removed_hyperedges(0),
    _max_edge_size(0),
    _num_pins(0),
    _num_graph_edges(0),
    _total_degree(0),
    _total_weight(0),
    _version(0),
    _contraction_index(0),
    _hypernodes(),
    _contraction_tree(),
    _incident_nets(),
    _acquired_hns(),
    _hyperedges(),
    _incidence_array(),
    _acquired_hes(),
    _tmp_incident_nets(),
    _failed_hyperedge_contractions(),
    _removable_incident_nets(),
    _removable_single_pin_and_parallel_nets(),
    _num_graph_edges_up_to(),
    _community_support() { }

  DynamicHypergraph(const DynamicHypergraph&) = delete;
  DynamicHypergraph & operator= (const DynamicHypergraph &) = delete;

  DynamicHypergraph(DynamicHypergraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_removed_hyperedges(other._num_removed_hyperedges),
    _max_edge_size(other._max_edge_size),
    _num_pins(other._num_pins),
    _num_graph_edges(other._num_graph_edges),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _version(other._version),
    _contraction_index(0),
    _hypernodes(std::move(other._hypernodes)),
    _contraction_tree(std::move(other._contraction_tree)),
    _incident_nets(std::move(other._incident_nets)),
    _acquired_hns(std::move(other._acquired_hns)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _acquired_hes(std::move(other._acquired_hes)),
    _tmp_incident_nets(std::move(other._tmp_incident_nets)),
    _failed_hyperedge_contractions(std::move(other._failed_hyperedge_contractions)),
    _removable_incident_nets(std::move(other._removable_incident_nets)),
    _removable_single_pin_and_parallel_nets(std::move(other._removable_single_pin_and_parallel_nets)),
    _num_graph_edges_up_to(std::move(other._num_graph_edges_up_to)),
    _community_support(std::move(other._community_support)) { }

  DynamicHypergraph & operator= (DynamicHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _max_edge_size = other._max_edge_size;
    _num_pins = other._num_pins;
    _num_graph_edges = other._num_graph_edges;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _version = other._version;
    _contraction_index.store(other._contraction_index.load());
    _hypernodes = std::move(other._hypernodes);
    _contraction_tree = std::move(other._contraction_tree);
    _incident_nets = std::move(other._incident_nets);
    _acquired_hns = std::move(other._acquired_hns);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
    _acquired_hes = std::move(other._acquired_hes);
    _tmp_incident_nets = std::move(other._tmp_incident_nets);
    _failed_hyperedge_contractions = std::move(other._failed_hyperedge_contractions);
    _removable_incident_nets = std::move(other._removable_incident_nets);
    _removable_single_pin_and_parallel_nets = std::move(other._removable_single_pin_and_parallel_nets);
    _num_graph_edges_up_to = std::move(other._num_graph_edges_up_to);
    _community_support = std::move(other._community_support);
    return *this;
  }

  ~DynamicHypergraph() {
    freeInternalData();
  }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
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

  HyperedgeID numGraphEdges() const {
    return _num_graph_edges;
  }

  HyperedgeID numNonGraphEdges() const {
    return initialNumEdges() - _num_graph_edges;
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

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
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
          if ( nodeIsEnabled(hn) ) {
            weight += this->_hypernodes[hn].weight();
          }
        }
        return weight;
      }, std::plus<HypernodeWeight>());
  }

  // ! Recomputes the total weight of the hypergraph (sequential)
  void updateTotalWeight() {
    _total_weight = 0;
    for ( const HypernodeID& hn : nodes() ) {
      if ( nodeIsEnabled(hn) ) {
        _total_weight += nodeWeight(hn);
      }
    }
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    static_cast<const DynamicHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
      if ( nodeIsEnabled(hn) ) {
        f(hn);
      }
    });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const DynamicHypergraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
      if ( edgeIsEnabled(he) ) {
        f(he);
      }
    });
  }

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return IteratorRange<HypernodeIterator>(
      HypernodeIterator(_hypernodes.data(), ID(0), _num_hypernodes),
      HypernodeIterator(_hypernodes.data() + _num_hypernodes, _num_hypernodes, _num_hypernodes));
  }

  // ! Returns a range of the active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return IteratorRange<HyperedgeIterator>(
      HyperedgeIterator(_hyperedges.data(), ID(0), _num_hyperedges),
      HyperedgeIterator(_hyperedges.data() + _num_hyperedges, _num_hyperedges, _num_hyperedges));
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_nets[u].c_it();
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

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_nets[u].size();
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

  // ! Maximum size of a hyperedge
  HypernodeID maxEdgeSize() const {
    return _max_edge_size;
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

  bool isGraphEdge(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const bool is_initial_graph_edge = _num_graph_edges_up_to[e + 1] - _num_graph_edges_up_to[e];
    ASSERT(!is_initial_graph_edge || edgeSize(e) <= 2);
    return is_initial_graph_edge && edgeSize(e) == 2;
  }

  HyperedgeID graphEdgeID(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(edgeSize(e) == 2);
    return _num_graph_edges_up_to[e];
  }

  HyperedgeID nonGraphEdgeID(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(edgeSize(e) > 2);
    return e - _num_graph_edges_up_to[e];
  }

  HypernodeID graphEdgeHead(const HyperedgeID e, const HypernodeID tail) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(edgeSize(e) == 2);
    const size_t f = hyperedge(e).firstEntry();
    const size_t first_matches = static_cast<size_t>(_incidence_array[f] == tail);
    return _incidence_array[f + first_matches];
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
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
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

  // ####################### Contract / Uncontract #######################

  DynamicHypergraph contract(parallel::scalable_vector<HypernodeID>&,
                             const TaskGroupID) {
    ERROR("contract(c, id) is not supported in dynamic hypergraph");
    return DynamicHypergraph();
  }

  /**!
   * Registers a contraction in the hypergraph whereas vertex u is the representative
   * of the contraction and v its contraction partner. Several threads can call this function
   * in parallel. The function adds the contraction of u and v to a contraction tree that determines
   * a parallel execution order and synchronization points for all running contractions.
   * The contraction can be executed by calling function contract(v, max_node_weight).
   */
  bool registerContraction(const HypernodeID u, const HypernodeID v) {
    // Acquires ownership of vertex v that gives the calling thread exclusive rights
    // to modify the contraction tree entry of v
    acquireHypernode(v);

    // If there is no other contraction registered for vertex v
    // we try to determine its representative in the contraction tree
    if ( _contraction_tree.parent(v) == v ) {

      HypernodeID w = u;
      bool cycle_detected = false;
      while ( true ) {
        // Search for representative of u in the contraction tree.
        // It is either a root of the contraction tree or a vertex
        // with a reference count greater than zero, which indicates
        // that there are still ongoing contractions on this node that
        // have to be processed.
        while ( _contraction_tree.parent(w) != w &&
                _contraction_tree.pendingContractions(w) == 0 ) {
          w = _contraction_tree.parent(w);
          if ( w == v ) {
            cycle_detected = true;
            break;
          }
        }

        if ( !cycle_detected ) {
          // In case contraction of u and v does not induce any
          // cycle in the contraction tree we try to acquire vertex w
          if ( w < v ) {
            // Acquire ownership in correct order to prevent deadlocks
            releaseHypernode(v);
            acquireHypernode(w);
            acquireHypernode(v);
            if ( _contraction_tree.parent(v) != v ) {
              releaseHypernode(v);
              releaseHypernode(w);
              return false;
            }
          } else {
            acquireHypernode(w);
          }

          // Double-check condition of while loop above after acquiring
          // ownership of w
          if ( _contraction_tree.parent(w) != w &&
               _contraction_tree.pendingContractions(w) == 0 ) {
            // In case something changed, we release ownership of w and
            // search again for the representative of u.
            releaseHypernode(w);
          } else {
            // Otherwise we perform final cycle check to verify that
            // contraction of u and v will not introduce any new cycle.
            HypernodeID x = w;
            do {
              x = _contraction_tree.parent(x);
              if ( x == v ) {
                cycle_detected = true;
                break;
              }
            } while ( _contraction_tree.parent(x) != x );

            if ( cycle_detected ) {
              releaseHypernode(w);
              releaseHypernode(v);
              return false;
            }

            // All checks succeded, we can safely increment the
            // reference count of w and update the contraction tree
            break;
          }
        } else {
          releaseHypernode(v);
          return false;
        }
      }

      // Increment reference count of w indicating that there pending
      // contraction at vertex w and update contraction tree.
      _contraction_tree.registerContraction(w, v, _version);

      releaseHypernode(w);
      releaseHypernode(v);
      return true;
    } else {
      releaseHypernode(v);
      return false;
    }
  }

  /**!
   * Contracts a previously registered contraction. Representative u of vertex v is looked up
   * in the contraction tree and performed if there are no pending contractions in the subtree
   * of v and the contractions respects the maximum allowed node weight. If (u,v) is the last
   * pending contraction in the subtree of u then the function recursively contracts also
   * u (if any contraction is registered). Therefore, function can return several contractions
   * or also return an empty contraction vector.
   */
  size_t contract(const HypernodeID v,
                  const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
    ASSERT(_contraction_tree.parent(v) != v, "No contraction registered for hypernode" << v);

    HypernodeID x = _contraction_tree.parent(v);
    HypernodeID y = v;
    ContractionResult res = ContractionResult::CONTRACTED;
    size_t num_contractions = 0;
    // We perform all contractions registered in the contraction tree
    // as long as there are no pending contractions
    while ( x != y && res != ContractionResult::PENDING_CONTRACTIONS) {
      // Perform Contraction
      res = contract(x, y, max_node_weight);
      if ( res == ContractionResult::CONTRACTED ) {
        ++num_contractions;
      }
      y = x;
      x = _contraction_tree.parent(y);
    }
    return num_contractions;
  }

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   * The two uncontraction functions are required by the partitioned hypergraph to restore
   * pin counts and gain cache values.
   */
  void uncontract(const Batch& batch,
                  const UncontractionFunction& case_one_func = NOOP_BATCH_FUNC,
                  const UncontractionFunction& case_two_func = NOOP_BATCH_FUNC) {
    ASSERT(batch.size() > 0UL);
    ASSERT([&] {
      const HypernodeID expected_batch_index = hypernode(batch[0].v).batchIndex();
      for ( const Memento& memento : batch ) {
        if ( hypernode(memento.v).batchIndex() != expected_batch_index ) {
          LOG << "Batch contains uncontraction from different batches."
              << "Hypernode" << memento.v << "with version" << hypernode(memento.v).batchIndex()
              << "but expected is" << expected_batch_index;
          return false;
        }
        if ( _contraction_tree.version(memento.v) != _version ) {
          LOG << "Batch contains uncontraction from a different version."
              << "Hypernode" << memento.v << "with version" << _contraction_tree.version(memento.v)
              << "but expected is" << _version;
          return false;
        }
      }
      return true;
    }(), "Batch contains uncontractions from different batches or from a different hypergraph version");

    tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
      const Memento& memento = batch[i];
      ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode" << memento.u << "is disabled");
      ASSERT(hypernode(memento.v).isDisabled(), "Hypernode" << memento.v << "is not invalid");

      parallel::scalable_vector<HyperedgeID>& failed_hyperedge_uncontractions = _failed_hyperedge_contractions.local();
      parallel::scalable_vector<bool>& removable_incident_nets_of_u = _removable_incident_nets.local();
      for ( const HyperedgeID& he : incidentEdges(memento.v) ) {
        // Try to acquire ownership of hyperedge. In case of success, we perform the
        // uncontraction and otherwise, we remember the hyperedge and try later again.
        if ( tryAcquireHyperedge(he) ) {
          uncontractHyperedge(memento.u, memento.v, he,
            removable_incident_nets_of_u, case_one_func, case_two_func);
          releaseHyperedge(he);
        } else {
          failed_hyperedge_uncontractions.push_back(he);
        }
      }

      // Perform uncontractions on which we failed to acquire ownership on the first try
      for ( const HyperedgeID& he : failed_hyperedge_uncontractions ) {
        acquireHyperedge(he);
        uncontractHyperedge(memento.u, memento.v, he,
          removable_incident_nets_of_u, case_one_func, case_two_func);
        releaseHyperedge(he);
      }

      failed_hyperedge_uncontractions.clear();
      acquireHypernode(memento.u);
      // Restore hypernode v which includes enabling it and subtract its weight
      // from its representative
      hypernode(memento.v).enable();
      hypernode(memento.u).setWeight(hypernode(memento.u).weight() - hypernode(memento.v).weight());

      // Remove all incident nets of u marked as removable
      IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[memento.u];
      size_t incident_nets_size = incident_nets_of_u.size();
      ASSERT(incident_nets_of_u.active_iterators() == 0);
      for ( size_t i = 0; i < incident_nets_size; ++i ) {
        const HyperedgeID he = incident_nets_of_u[i];
        ASSERT(he < removable_incident_nets_of_u.size());
        if ( removable_incident_nets_of_u[he] ) {
          std::swap(incident_nets_of_u[i--], incident_nets_of_u[--incident_nets_size]);
          incident_nets_of_u.pop_back();
          removable_incident_nets_of_u[he] = false;
        }
      }
      releaseHypernode(memento.u);
    });
  }

  /**
   * Computes a batch uncontraction hierarchy. A batch is a vector of mementos
   * (uncontractions) that are uncontracted in parallel. The function returns a vector
   * of versioned batch vectors. A new version of the hypergraph is induced if we perform
   * single-pin and parallel net detection. Once we process all batches of a versioned
   * batch vector, we have to restore all previously removed single-pin and parallel nets
   * in order to uncontract the next batch vector. We create for each version of the
   * hypergraph a seperate batch uncontraction hierarchy (see createBatchUncontractionHierarchyOfVersion(...))
   */
  VersionedBatchVector createBatchUncontractionHierarchy(const size_t batch_size,
                                                         const bool test = false) {
    const size_t num_versions = _version + 1;
    utils::Timer::instance().start_timer("finalize_contraction_tree", "Finalize Contraction Tree");
    // Finalizes the contraction tree such that it is traversable in a top-down fashion
    // and contains subtree size for each  tree node
    _contraction_tree.finalize(num_versions);
    utils::Timer::instance().stop_timer("finalize_contraction_tree");

    utils::Timer::instance().start_timer("create_versioned_batches", "Create Versioned Batches");
    VersionedBatchVector versioned_batches(num_versions);
    parallel::scalable_vector<size_t> batch_sizes_prefix_sum(num_versions, 0);
    tbb::parallel_for(0UL, num_versions, [&](const size_t version) {
      versioned_batches[version] =
        createBatchUncontractionHierarchyForVersion(batch_size, version);
    });
    for ( size_t version = 0; version < num_versions; ++version ) {
      if ( version > 0 ) {
        batch_sizes_prefix_sum[version] =
          batch_sizes_prefix_sum[version - 1] + versioned_batches[version - 1].size();
      }
    }
    utils::Timer::instance().stop_timer("create_versioned_batches");

    if ( !test ) {
      utils::Timer::instance().start_timer("prepare_hg_for_uncontraction", "Prepare HG For Uncontraction");

      // Store the batch index of each vertex in its hypernode data structure
      tbb::parallel_for(0UL, num_versions, [&](const size_t version) {
        tbb::parallel_for(0UL, versioned_batches[version].size(), [&](const size_t local_batch_idx) {
          const size_t batch_idx = batch_sizes_prefix_sum[version] + local_batch_idx;
          for ( const Memento& memento : versioned_batches[version][local_batch_idx] ) {
            hypernode(memento.v).setBatchIndex(batch_idx);
          }
        });
      });

      // Sort the invalid part of each hyperedge according to the batch indices of its pins
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
        const size_t first_invalid_entry = hyperedge(he).firstInvalidEntry();
        const size_t last_invalid_entry = hyperedge(he + 1).firstEntry();
        std::sort(_incidence_array.begin() + first_invalid_entry,
                  _incidence_array.begin() + last_invalid_entry,
                  [&](const HypernodeID u, const HypernodeID v) {
                    ASSERT(hypernode(u).batchIndex() != std::numeric_limits<HypernodeID>::max(),
                      "Hypernode" << u << "is not contained in the uncontraction hierarchy");
                    ASSERT(hypernode(v).batchIndex() != std::numeric_limits<HypernodeID>::max(),
                      "Hypernode" << v << "is not contained in the uncontraction hierarchy");
                    return hypernode(u).batchIndex() > hypernode(v).batchIndex();
                  });
      });
      utils::Timer::instance().stop_timer("prepare_hg_for_uncontraction");
    }

    return versioned_batches;
  }

  // ! Only for testing
  VersionedBatchVector createBatchUncontractionHierarchy(ContractionTree&& tree,
                                                         const size_t batch_size,
                                                         const size_t num_versions = 1) {
    ASSERT(num_versions > 0);
    _version = num_versions - 1;
    _contraction_tree = std::move(tree);
    return createBatchUncontractionHierarchy(batch_size, true);
  }

  // ! Only for testing
  HypernodeID contractionTree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _contraction_tree.parent(u);
  }

  // ! Only for testing
  HypernodeID pendingContractions(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _contraction_tree.pendingContractions(u);
  }

  // ! Only for testing
  void decrementPendingContractions(const HypernodeID u) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    _contraction_tree.decrementPendingContractions(u);
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

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge. Note, in contrast to removeEdge, this function
  * removes hyperedge from all its pins in parallel.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeLargeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      removeIncidentEdgeFromHypernode(he, pin);
    });
    disableHyperedge(he);
  }

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "is enabled");
    enableHyperedge(he);
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      connectHypernodeWithIncidentEdge(he, pin);
    });
  }

  /**
   * Removes single-pin and parallel nets from the hypergraph. The weight
   * of a set of identical nets is aggregated in one representative hyperedge
   * and single-pin hyperedges are removed. Returns a vector of removed hyperedges.
   */
  parallel::scalable_vector<ParallelHyperedge> removeSinglePinAndParallelHyperedges() {
    _removable_single_pin_and_parallel_nets.reset();
    // Remove singple-pin hyperedges directly from the hypergraph and
    // insert all other hyperedges into a bucket data structure such that
    // hyperedges with the same hash/footprint are placed in the same bucket.
    utils::Timer::instance().start_timer("preprocess_hyperedges", "Preprocess Hyperedges");
    StreamingVector<ParallelHyperedge> tmp_removed_hyperedges;
    ConcurrentBucketMap<ContractedHyperedgeInformation> hyperedge_hash_map;
    hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);
    doParallelForAllEdges([&](const HyperedgeID& he) {
      const HypernodeID edge_size = edgeSize(he);
      if ( edge_size > 1 ) {
        const Hyperedge& e = hyperedge(he);
        const size_t footprint = e.hash();
        std::sort(_incidence_array.begin() + e.firstEntry(),
                  _incidence_array.begin() + e.firstInvalidEntry());
        hyperedge_hash_map.insert(footprint,
          ContractedHyperedgeInformation { he, footprint, edge_size, true });
      } else {
        hyperedge(he).disable();
        _removable_single_pin_and_parallel_nets.set(he, true);
        tmp_removed_hyperedges.stream(ParallelHyperedge { he, kInvalidHyperedge });
      }
    });
    utils::Timer::instance().stop_timer("preprocess_hyperedges");

    // Helper function that checks if two hyperedges are parallel
    // Note, pins inside the hyperedges are sorted.
    auto check_if_hyperedges_are_parallel = [&](const HyperedgeID lhs,
                                                const HyperedgeID rhs) {
      const Hyperedge& lhs_he = hyperedge(lhs);
      const Hyperedge& rhs_he = hyperedge(rhs);
      if ( lhs_he.size() == rhs_he.size() ) {
        const size_t lhs_start = lhs_he.firstEntry();
        const size_t rhs_start = rhs_he.firstEntry();
        for ( size_t i = 0; i < lhs_he.size(); ++i ) {
          const size_t lhs_pos = lhs_start + i;
          const size_t rhs_pos = rhs_start + i;
          if ( _incidence_array[lhs_pos] != _incidence_array[rhs_pos] ) {
            return false;
          }
        }
        return true;
      } else {
        return false;
      }
    };

    // In the step before we placed hyperedges within a bucket data structure.
    // Hyperedges with the same hash/footprint are stored inside the same bucket.
    // We iterate now in parallel over each bucket and sort each bucket
    // after its hash. A bucket is processed by one thread and parallel
    // hyperedges are detected by comparing the pins of hyperedges with
    // the same hash.
    utils::Timer::instance().start_timer("parallel_net_detection", "Parallel Net Detection");
    tbb::parallel_for(0UL, hyperedge_hash_map.numBuckets(), [&](const size_t bucket) {
      auto& hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
      std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
        [&](const ContractedHyperedgeInformation& lhs, const ContractedHyperedgeInformation& rhs) {
          return lhs.hash < rhs.hash || (lhs.hash == rhs.hash && lhs.size < rhs.size)||
            (lhs.hash == rhs.hash && lhs.size == rhs.size && lhs.he < rhs.he);
        });

      // Parallel Hyperedge Detection
      for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
        ContractedHyperedgeInformation& contracted_he_lhs = hyperedge_bucket[i];
        if ( contracted_he_lhs.valid ) {
          const HyperedgeID lhs_he = contracted_he_lhs.he;
          HyperedgeWeight lhs_weight = hyperedge(lhs_he).weight();
          for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
            ContractedHyperedgeInformation& contracted_he_rhs = hyperedge_bucket[j];
            const HyperedgeID rhs_he = contracted_he_rhs.he;
            if ( contracted_he_rhs.valid &&
                 contracted_he_lhs.hash == contracted_he_rhs.hash &&
                 check_if_hyperedges_are_parallel(lhs_he, rhs_he) ) {
                // Hyperedges are parallel
                lhs_weight += hyperedge(rhs_he).weight();
                hyperedge(rhs_he).disable();
                _removable_single_pin_and_parallel_nets.set(rhs_he, true);
                contracted_he_rhs.valid = false;
                tmp_removed_hyperedges.stream( ParallelHyperedge { rhs_he, lhs_he } );
            } else if ( contracted_he_lhs.hash != contracted_he_rhs.hash  ) {
              // In case, hash of both are not equal we go to the next hyperedge
              // because we compared it with all hyperedges that had an equal hash
              break;
            }
          }
          hyperedge(lhs_he).setWeight(lhs_weight);
        }
      }
      hyperedge_hash_map.free(bucket);
    });
    utils::Timer::instance().stop_timer("parallel_net_detection");

    // Remove single-pin and parallel nets from incident net vector of vertices
    utils::Timer::instance().start_timer("postprocess_incident_nets", "Postprocess Incident Nets");
    doParallelForAllNodes([&](const HypernodeID& u) {
      IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[u];
      size_t incident_nets_size = incident_nets_of_u.size();
      for ( size_t pos = 0; pos < incident_nets_size; ++pos ) {
        if ( _removable_single_pin_and_parallel_nets[incident_nets_of_u[pos]] ) {
          std::swap(incident_nets_of_u[pos--], incident_nets_of_u[--incident_nets_size]);
          incident_nets_of_u.pop_back();
        }
      }
    });
    utils::Timer::instance().stop_timer("postprocess_incident_nets");

    utils::Timer::instance().start_timer("store_removed_hyperedges", "Store Removed Hyperedges");
    parallel::scalable_vector<ParallelHyperedge> removed_hyperedges = tmp_removed_hyperedges.copy_parallel();
    tmp_removed_hyperedges.clear_parallel();
    utils::Timer::instance().stop_timer("store_removed_hyperedges");

    ++_version;
    return removed_hyperedges;
  }

  /**
   * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
   * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
   */
  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore) {
    // Restores all previously removed hyperedges
    utils::Timer::instance().start_timer("restore_removed_nets", "Restore Removed Hyperedges");
    ConcurrentBucketMap<Memento> incident_net_map;
    tbb::parallel_for(0UL, hes_to_restore.size(), [&](const size_t i) {
      const ParallelHyperedge& parallel_he = hes_to_restore[i];
      const HyperedgeID he = parallel_he.removed_hyperedge;
      ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "should be disabled");
      const bool is_parallel_net = parallel_he.representative != kInvalidHyperedge;
      hyperedge(he).enable();
      if ( is_parallel_net ) {
        const HyperedgeID rep = parallel_he.representative;
        ASSERT(edgeIsEnabled(rep), "Hyperedge" << rep << "should be enabled");
        Hyperedge& rep_he = hyperedge(rep);
        acquireHyperedge(rep);
        rep_he.setWeight(rep_he.weight() - hyperedge(he).weight());
        releaseHyperedge(rep);
      }

      for ( const HypernodeID& pin : pins(he) ) {
        incident_net_map.insert(pin, Memento { pin, he });
      }
    });
    utils::Timer::instance().stop_timer("restore_removed_nets");

    // Adds all restored hyperedges as incident net to its contained pins.
    // In the previous step we inserted all pins together with its hyperedges
    // into a bucket data structure. All hyperedges of the same pin are placed
    // within the same bucket and can be processed sequentially here without
    // locking.
    utils::Timer::instance().start_timer("add_to_incident_nets", "Add Restored HEs to Incident Nets");
    tbb::parallel_for(0UL, incident_net_map.numBuckets(), [&](const size_t bucket) {
      auto& incident_net_bucket = incident_net_map.getBucket(bucket);
      std::sort(incident_net_bucket.begin(), incident_net_bucket.end(),
        [&](const Memento& lhs, const Memento& rhs) {
          return lhs.u < rhs.u || ( lhs.u == rhs.u && lhs.v < rhs.v);
        });

      // No locking required since vertex u can only occur in one bucket
      for ( const Memento& memento : incident_net_bucket ) {
        const HypernodeID u = memento.u;
        const HyperedgeID he = memento.v;
        _incident_nets[u].push_back(he);
      }

      incident_net_map.free(bucket);
    });
    utils::Timer::instance().stop_timer("add_to_incident_nets");

    --_version;
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
  void initializeCommunities() {
    _community_support.initialize(*this);
  }

  /*!
  * Initializes community hyperedges.
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to a range of consecutive pins with
  *       same community in that hyperedge
  */
  void initializeCommunityHyperedges(const TaskGroupID) {
    _community_support.initializeCommunityHyperedges(*this);
  }

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   */
  void removeCommunityHyperedges(const TaskGroupID,
                                 const parallel::scalable_vector<HypernodeID>& contraction_index = {}) {
    _community_support.removeCommunityHyperedges(contraction_index);
  }

  // ! Reset internal community information
  void setCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
    if ( _community_support.isInitialized() ) {
      _community_support.freeInternalData();
    }

    ASSERT(community_ids.size() == UI64(_num_hypernodes));
    doParallelForAllNodes([&](const HypernodeID& hn) {
      hypernode(hn).setCommunityID(community_ids[hn]);
    });

    initializeCommunities();
  }

  // ####################### Copy #######################

  // ! Copy dynamic hypergraph in parallel
  DynamicHypergraph copy(const TaskGroupID task_group_id) {
    DynamicHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._num_graph_edges = _num_graph_edges;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;
    hypergraph._version = _version;
    hypergraph._contraction_index.store(_contraction_index.load());

    tbb::parallel_invoke([&] {
      hypergraph._hypernodes.resize(_hypernodes.size());
      memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
        sizeof(Hypernode) * _hypernodes.size());
    }, [&] {
      tbb::parallel_invoke([&] {
        hypergraph._incident_nets.resize(_incident_nets.size());
      }, [&] {
        hypergraph._acquired_hns.resize(_acquired_hns.size());
      });
      tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
        hypergraph._acquired_hns[hn] = _acquired_hns[hn];
        hypergraph._incident_nets[hn].resize(_incident_nets[hn].size());
        memcpy(hypergraph._incident_nets[hn].data(), _incident_nets[hn].data(),
          sizeof(HyperedgeID) * _incident_nets[hn].size());
      });
    }, [&] {
      hypergraph._contraction_tree = _contraction_tree.copy(task_group_id);
    }, [&] {
      hypergraph._hyperedges.resize(_hyperedges.size());
      memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
        sizeof(Hyperedge) * _hyperedges.size());
    }, [&] {
      hypergraph._incidence_array.resize(_incidence_array.size());
      memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
        sizeof(HypernodeID) * _incidence_array.size());
    }, [&] {
      hypergraph._acquired_hes.resize(_num_hyperedges);
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
        hypergraph._acquired_hes[he] = _acquired_hes[he];
      });
    }, [&] {
      hypergraph._removable_incident_nets = _removable_incident_nets;
    }, [&] {
      hypergraph._removable_single_pin_and_parallel_nets =
        kahypar::ds::FastResetFlagArray<>(_num_hyperedges);
    }, [&] {
      hypergraph._num_graph_edges_up_to.resize(_num_graph_edges_up_to.size());
      memcpy(hypergraph._num_graph_edges_up_to.data(), _num_graph_edges_up_to.data(),
             sizeof(HyperedgeID) * _num_graph_edges_up_to.size());
    }, [&] {
      hypergraph._community_support = _community_support.copy(task_group_id);
    });
    return hypergraph;
  }

  // ! Copy dynamic hypergraph sequential
  DynamicHypergraph copy() {
    DynamicHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._num_graph_edges = _num_graph_edges;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;
    hypergraph._version = _version;
    hypergraph._contraction_index.store(_contraction_index.load());

    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
      sizeof(Hypernode) * _hypernodes.size());
    hypergraph._incident_nets.resize(_incident_nets.size());
    hypergraph._acquired_hns.resize(_num_hypernodes);
    for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
      hypergraph._incident_nets[hn].resize(_incident_nets[hn].size());
      hypergraph._acquired_hns[hn] = _acquired_hns[hn];
      memcpy(hypergraph._incident_nets[hn].data(), _incident_nets[hn].data(),
        sizeof(HyperedgeID) * _incident_nets[hn].size());
    }
    hypergraph._contraction_tree = _contraction_tree.copy();
    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
      sizeof(Hyperedge) * _hyperedges.size());
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
      sizeof(HypernodeID) * _incidence_array.size());
    hypergraph._acquired_hes.resize(_num_hyperedges);
    for ( HyperedgeID he = 0; he < _num_hyperedges; ++he ) {
      hypergraph._acquired_hes[he] = _acquired_hes[he];
    }
    hypergraph._removable_incident_nets = _removable_incident_nets;
    hypergraph._removable_single_pin_and_parallel_nets =
      kahypar::ds::FastResetFlagArray<>(_num_hyperedges);
    hypergraph._num_graph_edges_up_to.resize(_num_graph_edges_up_to.size());
    memcpy(hypergraph._num_graph_edges_up_to.data(), _num_graph_edges_up_to.data(),
            sizeof(HyperedgeID) * _num_graph_edges_up_to.size());

    hypergraph._community_support = _community_support.copy();

    return hypergraph;
  }

  // ! Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      _community_support.freeInternalData();
    }
    _num_hypernodes = 0;
    _num_hyperedges = 0;
  }

  void freeTmpContractionBuffer() {
    ERROR("freeTmpContractionBuffer() is not supported in dynamic hypergraph");
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
    parent->addChild("Incident Nets", sizeof(HyperedgeID) * _incidence_array.size());
    parent->addChild("Hypernode Ownership Vector", sizeof(bool) * _acquired_hns.size());
    parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
    parent->addChild("Hyperedge Ownership Vector", sizeof(bool) * _acquired_hes.size());
    parent->addChild("Bitsets",
      ( _num_hyperedges * _removable_incident_nets.size() ) / 8UL + sizeof(uint16_t) * _num_hyperedges);
    parent->addChild("Graph Edge ID Mapping", sizeof(HyperedgeID) * _num_graph_edges_up_to.size());

    utils::MemoryTreeNode* contraction_tree_node = parent->addChild("Contraction Tree");
    _contraction_tree.memoryConsumption(contraction_tree_node);
    utils::MemoryTreeNode* community_support_node = parent->addChild("Community Support");
    _community_support.memoryConsumption(community_support_node);
  }

  // ! Only for testing
  bool verifyIncidenceArrayAndIncidentNets() {
    bool success = true;
    tbb::parallel_invoke([&] {
      doParallelForAllNodes([&](const HypernodeID& hn) {
        for ( const HyperedgeID& he : incidentEdges(hn) ) {
          bool found = false;
          for ( const HypernodeID& pin : pins(he) ) {
            if ( pin == hn ) {
              found = true;
              break;
            }
          }
          if ( !found ) {
            LOG << "Hypernode" << hn << "not found in incidence array of net" << he;
            success = false;
          }
        }
      });
    }, [&] {
      doParallelForAllEdges([&](const HyperedgeID& he) {
        for ( const HypernodeID& pin : pins(he) ) {
          bool found = false;
          for ( const HyperedgeID& e : incidentEdges(pin) ) {
            if ( e == he ) {
              found = true;
              break;
            }
          }
          if ( !found ) {
            LOG << "Hyperedge" << he << "not found in incident nets of vertex" << pin;
            success = false;
          }
        }
      });
    });
    return success;
  }

 private:
  friend class DynamicHypergraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;
  template <typename Hypergraph,
            typename HypergraphFactory>
  friend class PartitionedHypergraph;

  // ####################### Acquiring / Releasing Ownership #######################

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hns[u].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hns[u].compare_exchange_strong(expected, desired);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    ASSERT(_acquired_hns[u], "Hypernode" << u << "is not acquired!");
    _acquired_hns[u] = false;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hes[e].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hes[e].compare_exchange_strong(expected, desired);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    ASSERT(_acquired_hes[e], "Hyperedge" << e << "is not acquired!");
    _acquired_hes[e] = false;
  }

  // ####################### Hypernode Information #######################

  // ! Accessor for hypernode-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const DynamicHypergraph&>(*this).hypernode(u));
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const IncidentNetsIterator incident_nets_of(const HypernodeID u) const {
    return _incident_nets[u].cbegin();
  }

  // ####################### Hyperedge Information #######################

  // ! Accessor for hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
    ASSERT(e <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[e];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const DynamicHypergraph&>(*this).hyperedge(e));
  }

  // ####################### Contract / Uncontract #######################

  /**!
   * Contracts a previously registered contraction. The contraction of u and v is
   * performed if there are no pending contractions in the subtree of v and the
   * contractions respects the maximum allowed node weight. In case the contraction
   * was performed successfully, enum type CONTRACTED is returned. If contraction
   * was not performed either WEIGHT_LIMIT_REACHED (in case sum of both vertices is
   * greater than the maximum allowed node weight) or PENDING_CONTRACTIONS (in case
   * there are some unfinished contractions in the subtree of v) is returned.
   */
  ContractionResult contract(const HypernodeID u,
                             const HypernodeID v,
                             const HypernodeWeight max_node_weight) {

    // Acquire ownership in correct order to prevent deadlocks
    if ( u < v ) {
      acquireHypernode(u);
      acquireHypernode(v);
    } else {
      acquireHypernode(v);
      acquireHypernode(u);
    }

    // Contraction is valid if
    //  1.) Contraction partner v is enabled
    //  2.) There are no pending contractions on v
    //  4.) Resulting node weight is less or equal than a predefined upper bound
    const bool contraction_partner_valid = nodeIsEnabled(v) && _contraction_tree.pendingContractions(v) == 0;
    const bool less_or_equal_than_max_node_weight =
      hypernode(u).weight() + hypernode(v).weight() <= max_node_weight;
    if ( contraction_partner_valid && less_or_equal_than_max_node_weight ) {
      ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled!");
      hypernode(u).setWeight(nodeWeight(u) + nodeWeight(v));
      hypernode(v).disable();
      releaseHypernode(u);
      releaseHypernode(v);

      HypernodeID contraction_start = _contraction_index.load();
      parallel::scalable_vector<HyperedgeID>& tmp_incident_nets = _tmp_incident_nets.local();
      parallel::scalable_vector<HyperedgeID>& failed_hyperedge_contractions = _failed_hyperedge_contractions.local();
      for ( const HyperedgeID& he : incidentEdges(v) ) {
        // Try to acquire ownership of hyperedge. In case of success, we perform the
        // contraction and otherwise, we remember the hyperedge and try later again.
        if ( tryAcquireHyperedge(he) ) {
          contractHyperedge(u, v, he, tmp_incident_nets);
          releaseHyperedge(he);
        } else {
          failed_hyperedge_contractions.push_back(he);
        }
      }

      // Perform contraction on which we failed to acquire ownership on the first try
      for ( const HyperedgeID& he : failed_hyperedge_contractions ) {
        acquireHyperedge(he);
        contractHyperedge(u, v, he, tmp_incident_nets);
        releaseHyperedge(he);
      }

      // tmp_incident_nets contains all hyperedges to which vertex u is
      // adjacent after the contraction, we use a special insert function to make
      // sure that iterators are not invalidated while inserting into the vector.
      if ( tmp_incident_nets.size() > 0 ) {
        _incident_nets[u].bulk_insert(tmp_incident_nets);
      }
      tmp_incident_nets.clear();
      failed_hyperedge_contractions.clear();

      HypernodeID contraction_end = ++_contraction_index;
      acquireHypernode(u);
      _contraction_tree.unregisterContraction(u, v, contraction_start, contraction_end);
      releaseHypernode(u);
      return ContractionResult::CONTRACTED;
    } else {
      ContractionResult res = ContractionResult::PENDING_CONTRACTIONS;
      if ( !less_or_equal_than_max_node_weight && nodeIsEnabled(v) ) {
        _contraction_tree.unregisterContraction(u, v,
          kInvalidHypernode, kInvalidHypernode, true /* failed */);
        res = ContractionResult::WEIGHT_LIMIT_REACHED;
      }
      releaseHypernode(u);
      releaseHypernode(v);
      return res;
    }
  }

  // ! Performs the contraction of (u,v) inside hyperedge he
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void contractHyperedge(const HypernodeID u, const HypernodeID v, const HyperedgeID he,
                                                         parallel::scalable_vector<HyperedgeID>& tmp_incident_nets) {
    Hyperedge& e = hyperedge(he);
    const HypernodeID pins_begin = e.firstEntry();
    const HypernodeID pins_end = e.firstInvalidEntry();
    HypernodeID slot_of_u = pins_end - 1;
    HypernodeID last_pin_slot = pins_end - 1;

    for (HypernodeID idx = pins_begin; idx != last_pin_slot; ++idx) {
      const HypernodeID pin = _incidence_array[idx];
      if (pin == v) {
        std::swap(_incidence_array[idx], _incidence_array[last_pin_slot]);
        --idx;
      } else if (pin == u) {
        slot_of_u = idx;
      }
    }

    ASSERT(_incidence_array[last_pin_slot] == v, "v is not last entry in adjacency array!");

    if (slot_of_u != last_pin_slot) {
      // Case 1:
      // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
      // can just cut off the last entry in the edge array of e that now contains v.
      DBG << V(he) << ": Case 1";
      e.hash() -= kahypar::math::hash(v);
      e.decrementSize();
    } else {
      DBG << V(he) << ": Case 2";
      // Case 2:
      // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
      // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
      e.hash() -= kahypar::math::hash(v);
      e.hash() += kahypar::math::hash(u);
      _incidence_array[last_pin_slot] = u;
      tmp_incident_nets.push_back(he);
    }
  }

  // ! Uncontracts u and v in hyperedge he.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void uncontractHyperedge(const HypernodeID u,
                                                           const HypernodeID v,
                                                           const HyperedgeID he,
                                                           parallel::scalable_vector<bool>& removable_incident_nets_of_u,
                                                           const UncontractionFunction& case_one_func,
                                                           const UncontractionFunction& case_two_func) {
    const HypernodeID batch_index = hypernode(v).batchIndex();

    // We search for hypernode v in the invalid part of hyperedge he. Note, that
    // pins in the invalid part of the hyperedge are sorted in decreasing order
    // of their batch index which they are contained in as contraction partner.
    // Therefore, we only visit the pins that have the same batch index as v
    // in the invalid part.
    const size_t first_invalid_entry = hyperedge(he).firstInvalidEntry();
    const size_t last_invalid_entry = hyperedge(he + 1).firstEntry();
    size_t slot_of_v = last_invalid_entry;
    for ( size_t pos = first_invalid_entry; pos < last_invalid_entry; ++pos ) {
      const HypernodeID pin = _incidence_array[pos];
      ASSERT(hypernode(pin).batchIndex() <= batch_index, V(he));
      if ( pin == v ) {
        slot_of_v = pos;
        break;
      } else if ( hypernode(pin).batchIndex() != batch_index ) {
        break;
      }
    }

    if ( slot_of_v != last_invalid_entry ) {
      // In that case v was found in the invalid part and u and v were part of hyperedge
      // he also before its contraction. Therefore, we swap v to the first invalid
      // entry and increment the size of the hyperedge.
      std::swap(_incidence_array[first_invalid_entry], _incidence_array[slot_of_v]);
      hyperedge(he).incrementSize();
      case_one_func(u, v, he);
    } else {
      // In that case v was not found in the invalid part of hyperedge he, which means that
      // u was not incident to hyperedge he before the contraction. Therefore, we have to find
      // u and replace it with v.
      const size_t first_valid_entry = hyperedge(he).firstEntry();
      size_t slot_of_u = first_invalid_entry;
      for ( size_t pos = first_invalid_entry - 1; pos != first_valid_entry - 1; --pos ) {
        if ( u == _incidence_array[pos] ) {
          slot_of_u = pos;
          break;
        }
      }

      ASSERT(slot_of_u != first_invalid_entry, "Hypernode" << u << "is not incident to hyperedge" << he);

      _incidence_array[slot_of_u] = v;
      // u is not incident to hyperedge he after this uncontraction
      // => we mark hyperedge he as removable and remove it in a postprocessing step
      removable_incident_nets_of_u[he] = true;
      case_two_func(u, v, he);
    }
  }

  /**
   * Computes a batch uncontraction hierarchy for a specific version of the hypergraph.
   * A batch is a vector of mementos (uncontractions) that are uncontracted in parallel.
   * Each time we perform single-pin and parallel net detection we create a new version of
   * the hypergraph.
   * A batch of uncontractions that is uncontracted in parallel must satisfy two conditions:
   *  1.) All representatives must be active vertices of the hypergraph
   *  2.) For a specific representative its contraction partners must be uncontracted in reverse
   *      contraction order. Meaning that a contraction (u, v) that happens before a contraction (u, w)
   *      must be uncontracted in a batch that is part of the same batch or a batch uncontracted after the
   *      batch which (u, w) is part of. This ensures that a parallel batch uncontraction does not
   *      increase the objective function.
   * We use the contraction tree to create a batch uncontraction order. Note, uncontractions from
   * different subtrees can be interleaved abitrary. To ensure condition 1.) we peform a BFS starting
   * from all roots of the contraction tree. Each BFS level induces a new batch. Since we contract
   * vertices in parallel its not possible to create a relative order of the contractions which is
   * neccassary for condition 2.). However, during a contraction we store a start and end "timepoint"
   * of a contraction. If two contractions time intervals do not intersect, we can determine
   * which contraction happens strictly before the other. If they intersect, it is not possible to
   * give a relative order. To ensure condition 2.) we sort the childs of a vertex in the contraction tree
   * after its time intervals. Once we add a uncontraction (u,v) to a batch, we also add all uncontractions
   * (u,w) to the batch which intersect the time interval of (u,v). To merge uncontractions of different
   * subtrees in a batch, we insert all eligble uncontractions into a max priority queue with the subtree
   * size of the contraction partner as key. We insert uncontractions into the current batch as long
   * as the maximum batch size is not reached or the PQ is empty. Once the batch reaches its maximum
   * batch size, we create a new empty batch. If the PQ is empty, we replace it with the PQ of the next
   * BFS level. With this approach heavy vertices are uncontracted earlier (subtree size in the PQ as key = weight of
   * a vertex for an unweighted hypergraph) such that average node weight of the hypergraph decreases faster and
   * local searches are more effective in early stages of the uncontraction hierarchy where hyperedge sizes are
   * usually smaller than on the original hypergraph.
   */
  BatchVector createBatchUncontractionHierarchyForVersion(const size_t batch_size,
                                                          const size_t version) {

    using PQ = std::priority_queue<PQBatchUncontractionElement,
                                   parallel::scalable_vector<PQBatchUncontractionElement>,
                                   PQElementComparator>;

    // Checks if two contraction intervals intersect
    auto does_interval_intersect = [&](const ContractionInterval& i1, const ContractionInterval& i2) {
      if (i1.start == kInvalidHypernode || i2.start == kInvalidHypernode) {
        return false;
      }
      return (i1.start <= i2.end && i1.end >= i2.end) ||
             (i2.start <= i1.end && i2.end >= i1.end);
    };

    auto push_into_pq = [&](PQ& prio_q, const HypernodeID& u) {
      auto it = _contraction_tree.childs(u);
      auto current = it.begin();
      auto end = it.end();
      while ( current != end && _contraction_tree.version(*current) != version ) {
        ++current;
      }
      if ( current != end ) {
        prio_q.push(PQBatchUncontractionElement {
          _contraction_tree.subtreeSize(*current), std::make_pair(current, end) } );
      }
    };

    BatchVector batches(1);
    // Contains eligble uncontractions of the current BFS level
    PQ pq;
    // Contains eligble uncontractions of the next BFS level
    PQ next_pq;

    // Insert all roots of the current version into the priority queue
    const parallel::scalable_vector<HypernodeID>& roots = _contraction_tree.roots_of_version(version);
    for ( const HypernodeID& root : roots ) {
      push_into_pq(pq, root);
    }

    while ( !pq.empty() ) {
      // Iterator over the childs of a active vertex
      auto it = pq.top()._iterator;
      ASSERT(it.first != it.second);
      const HypernodeID v = *it.first;
      pq.pop();

      // If the current batch reaches the maximum batch size we
      // create a new batch
      ASSERT(_contraction_tree.version(v) == version);
      if ( batches.back().size() >= batch_size ) {
        batches.emplace_back();
      }

      // Insert uncontraction (u,v) into the current batch
      const HypernodeID u = _contraction_tree.parent(v);
      batches.back().push_back(Memento { u, v });
      // Push contraction partner into pq for the next BFS level
      push_into_pq(next_pq, v);

      // Insert all childs of u that intersect the contraction time interval of
      // (u,v) into the current batch
      ++it.first;
      ContractionInterval current_ival = _contraction_tree.interval(v);
      while ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
        const HypernodeID w = *it.first;
        const ContractionInterval w_ival = _contraction_tree.interval(w);
        if ( does_interval_intersect(current_ival, w_ival) ) {
          ASSERT(_contraction_tree.parent(w) == u);
          batches.back().push_back(Memento { u, w });
          current_ival.start = std::min(current_ival.start, w_ival.start);
          current_ival.end = std::max(current_ival.end, w_ival.end);
          push_into_pq(next_pq, w);
        } else {
          break;
        }
        ++it.first;
      }

      // If there are still childs left of u, we push the iterator again into the
      // priority queue of the current BFS level.
      if ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
        pq.push(PQBatchUncontractionElement {
          _contraction_tree.subtreeSize(*it.first), it });
      }

      if ( pq.empty() ) {
        // Processing of current BFS level finishes
        // => create new batch and swap current BFS pq with next level BFS pq
        if ( !batches.back().empty() ) {
          batches.emplace_back();
        }
        std::swap(pq, next_pq);
      }
    }

    while ( !batches.empty() && batches.back().empty() ) {
      batches.pop_back();
    }
    std::reverse(batches.begin(), batches.end());

    return batches;
  }

  // ####################### Remove / Restore Hyperedges #######################

  // ! Removes hyperedge e from the incident nets of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID e,
                                                                       const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");

    IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[u];
    size_t pos = 0;
    for ( ; pos < incident_nets_of_u.size(); ++pos ) {
      if ( incident_nets_of_u[pos] == e ) {
        break;
      }
    }
    ASSERT(pos < incident_nets_of_u.size());
    swap(incident_nets_of_u[pos], incident_nets_of_u.back());
    incident_nets_of_u.pop_back();
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHypernodeWithIncidentEdge(const HyperedgeID e,
                                                                        const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[u];
    HEAVY_REFINEMENT_ASSERT(std::count(incident_nets_of_u.cbegin(), incident_nets_of_u.cend(), e) == 0,
                        "HN" << u << "is already connected to HE" << e);
    incident_nets_of_u.push_back(e);
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of removed hypernodes
  HypernodeID _num_removed_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
  // ! Maximum size of a hyperedge
  HypernodeID _max_edge_size;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Number of graph edges (hyperedges of size two)
  HyperedgeID _num_graph_edges;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;
  // ! Version of the hypergraph, each time we remove a single-pin and parallel nets,
  // ! we create a new version
  size_t _version;
  // ! Contraction Index, increment whenever a contraction terminates
  std::atomic<HypernodeID> _contraction_index;

  // ! Hypernodes
  Array<Hypernode> _hypernodes;
  // ! Contraction Tree
  ContractionTree _contraction_tree;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Atomic bool vector used to acquire unique ownership of hypernodes
  OwnershipVector _acquired_hns;


  // ! Hyperedges
  Array<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;
  // ! Atomic bool vector used to acquire unique ownership of hyperedges
  OwnershipVector _acquired_hes;
  // ! Collects hyperedes that will be adjacent to a vertex after a contraction
  ThreadLocalHyperedgeVector _tmp_incident_nets;
  // ! Collects hyperedge contractions that failed due to failed acquired ownership
  ThreadLocalHyperedgeVector _failed_hyperedge_contractions;
  // ! Marks incident nets that have to be removed from vertex u after an uncontraction (u,v)
  ThreadLocalBitset _removable_incident_nets;
  // ! Single-pin and parallel nets are marked within that vector during the algorithm
  kahypar::ds::FastResetFlagArray<> _removable_single_pin_and_parallel_nets;
  // ! Number of graph edges with smaller ID than the access ID
  Array<HyperedgeID> _num_graph_edges_up_to;

  // ! Community Information and Stats
  CommunitySupport<DynamicHypergraph> _community_support;

};

} // namespace ds
} // namespace mt_kahypar