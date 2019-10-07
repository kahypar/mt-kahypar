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

#include <algorithm>
#include <thread>
#include <type_traits>
#include <atomic>
#include <functional>
#include <chrono>
#include <set>

#include "tbb/task_scheduler_observer.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_scan.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"
#include "tbb/queuing_mutex.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"


namespace mt_kahypar {
namespace ds {

// forward
template <typename HypernodeType_,
          typename HyperedgeType_,
          typename HypernodeWeightType_,
          typename HyperedgeWeightType_,
          typename PartitionIDType_,
          typename HardwareTopology,
          typename TBBNumaArena>
class Hypergraph;

// ! Helper function to allow range-based for loops
template <typename Iterator>
Iterator begin(const std::pair<Iterator, Iterator>& x) {
  return x.first;
}

// ! Helper function to allow range-based for loops
template <typename Iterator>
Iterator end(const std::pair<Iterator, Iterator>& x) {
  return x.second;
}

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class StreamingHypergraph {

  static constexpr bool debug = false;
  static constexpr size_t NUMA_NODE_INDENTIFIER = 48;
  // seed for edge hashes used for parallel net detection
  static constexpr size_t kEdgeHashSeed = 42;

  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using PartitionID = PartitionIDType_;

  static constexpr PartitionID kInvalidPartition = -1;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

  using Self = StreamingHypergraph<HypernodeID, HyperedgeID, HypernodeWeight,
                                   HyperedgeWeight, PartitionID, HardwareTopology,
                                   TBBNumaArena>;

  using IncidentNets = parallel::scalable_vector<parallel::scalable_vector<HyperedgeID>>;

  using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

  static_assert( sizeof(HypernodeID) == 8, "Hypernode ID must be 8 byte" );
  static_assert( std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned" );
  static_assert( sizeof(HyperedgeID) == 8, "Hyperedge ID must be 8 byte" );
  static_assert( std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned" );

  class Hypernode {

    public:
      using IDType = HyperedgeID;

      Hypernode() :
        _id(kInvalidHypernode),
        _original_id(kInvalidHypernode),
        _community_node_id(kInvalidHypernode),
        _weight(1),
        _community_id(kInvalidPartition),
        _part_id(kInvalidPartition),
        _invalid_incident_nets(0),
        _single_pin_community_nets(0),
        _invalid_community_nets(0),
        _valid(false) { }

      Hypernode(const HypernodeID id,
                const HypernodeID original_id,
                const HypernodeWeight weight) :
        _id(id),
        _original_id(original_id),
        _community_node_id(kInvalidHypernode),
        _weight(weight),
        _community_id(kInvalidPartition),
        _part_id(kInvalidPartition),
        _invalid_incident_nets(0),
        _single_pin_community_nets(0),
        _invalid_community_nets(0),
        _valid(true) { }

      HypernodeID nodeId() const {
        return _id;
      }

      HypernodeID originalNodeId() const {
        return _original_id;
      }

      HypernodeID communityNodeId() const {
        return _community_node_id;
      }

      void setCommunityNodeId(const HypernodeID community_node_id) {
        _community_node_id = community_node_id;
      }

      // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
      // ! when iterating over the set of all nodes/edges.
      void disable() {
        ASSERT(!isDisabled());
        _valid = false;
      }

      bool isDisabled() const {
        return _valid == false;
      }

      void enable() {
        ASSERT(isDisabled());
        _valid = true;
      }

      IDType size() const {
        ASSERT(!isDisabled());
        return _incident_nets.size();
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
        //ASSERT(!isDisabled());
        return _community_id;
      }

      void setCommunityID(const PartitionID community_id) {
        ASSERT(!isDisabled());
        _community_id = community_id;
      }

      PartitionID partID() const {
        ASSERT(!isDisabled());
        return _part_id;
      }

      void setPartID(const PartitionID part_id) {
        ASSERT(!isDisabled());
        _part_id = part_id;
      }

      size_t invalidIncidentNets() const {
        return _invalid_incident_nets;
      }

      void setInvalidIncidentNets(const size_t invalid_incident_nets) {
        _invalid_incident_nets = invalid_incident_nets;
      }

      void incrementInvalidIncidentNets() {
        ++_invalid_incident_nets;
      }

      void decrementInvalidIncidentNets() {
        --_invalid_incident_nets;
      }

      size_t singlePinCommunityNets() const {
        ASSERT(!isDisabled());
        return _single_pin_community_nets;
      }

      void setSinglePinCommunityNets(const size_t single_pin_community_nets) {
        ASSERT(!isDisabled());
        _single_pin_community_nets = single_pin_community_nets;
      }

      void incrementSinglePinCommunityNets() {
        ASSERT(!isDisabled());
        ++_single_pin_community_nets;
      }

      void decrementSinglePinCommunityNets() {
        ASSERT(!isDisabled());
        --_single_pin_community_nets;
      }

      size_t invalidCommunityNets() const {
        ASSERT(!isDisabled());
        return _invalid_community_nets;
      }

      void setInvalidCommunityNets(const size_t invalid_community_nets) {
        ASSERT(!isDisabled());
        _invalid_community_nets = invalid_community_nets;
      }

      void incrementInvalidCommunityNets() {
        ASSERT(!isDisabled());
        ++_invalid_community_nets;
      }

      void decrementInvalidCommunityNets() {
        ASSERT(!isDisabled());
        --_invalid_community_nets;
      }

      bool operator== (const Hypernode& rhs) const {
        return _incident_nets.size() == rhs._incident_nets.size() &&
              _weight == rhs._weight &&
              _valid == rhs._valid &&
              std::is_permutation(_incident_nets.begin(),
                                  _incident_nets.end(),
                                  rhs._incident_nets.begin());
      }

      bool operator!= (const Hypernode& rhs) const {
        return !operator== (this, rhs);
      }

    private:
      // ! Hypernode id
      HypernodeID _id;
      // ! Original hypernode id
      HypernodeID _original_id;
      // ! Hypernode id within a community
      HypernodeID _community_node_id;
      // ! Hypernode weight
      HyperedgeWeight _weight;
      // ! Community id
      PartitionID _community_id;
      // ! Part Id
      PartitionID _part_id;
      // ! Pointer to incident nets indicating that all hyperedges up to
      // ! that position are invalid.
      size_t _invalid_incident_nets;

      // ! During community coarsening we maintain the following order
      // ! on the incident nets of a hypernode
      // !  1.) All single-pin community hyperedges
      // !  2.) All valid community hyperedges
      // !  3.) All invalidated community hyperedges (due to parallel to an other edge)
      // ! Reason for that is:
      // !  1.) During rating we only want to iterate over all valid community hyperedges
      // !  2.) During contraction we want to iterate over all community hyperedges
      // !  3.) During parallel hyperedge detection we only want to iterate over all
      // !      single-pin and valid hyperedges

      // ! All incident nets from 0 to _single_pin_community_nets are single-pin community hyperedges
      size_t _single_pin_community_nets;
      // ! All incident nets from _invalid_community_nets to |I(v)| are invalid community hyperedges
      size_t _invalid_community_nets;

      // ! Flag indicating whether or not the element is active.
      bool _valid;
  };

  class Hyperedge {

    public:
      using IDType = HyperedgeID;

      Hyperedge() :
        _begin(0),
        _size(0),
        _original_id(kInvalidHypernode),
        _weight(1),
        _hash(kEdgeHashSeed),
        _valid(false),
        _init_community_hyperedges(false) { }

      Hyperedge(const size_t begin,
                const size_t size,
                const HyperedgeID original_id,
                const HyperedgeWeight weight) :
        _begin(begin),
        _size(size),
        _original_id(original_id),
        _weight(weight),
        _hash(kEdgeHashSeed),
        _valid(true),
        _init_community_hyperedges(false) { }

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

      void initializeCommunityHyperedges() {
        ASSERT(!isInitCommunityHyperedges());
        _init_community_hyperedges = true;
      }

      void deinitializeCommunityHyperedges() {
        ASSERT(isInitCommunityHyperedges());
        _init_community_hyperedges = false;
      }

      bool isInitCommunityHyperedges() const {
        return _init_community_hyperedges;
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

      HypernodeID originalEdgeId() const {
        return _original_id;
      }

      void incrementSize() {
        ++_size;
      }

      void decrementSize() {
        ASSERT(_size > 0);
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

      size_t & hash() {
        return _hash;
      }

      size_t hash() const {
        return _hash;
      }

      bool operator== (const Hyperedge& rhs) const {
        return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
      }

      bool operator!= (const Hyperedge& rhs) const {
        return !operator== (this, rhs);
      }

    private:
      // ! Index of the first element in _incidence_array
      size_t _begin;
      // ! Number of _incidence_array elements
      size_t _size;
      // ! Original hyperedge id
      HyperedgeID _original_id;
      // ! hyperedge weight
      HyperedgeWeight _weight;
      // ! Hash of pins
      size_t _hash;
      // ! Flag indicating whether or not the element is active.
      bool _valid;
      // ! Flag indicating wheter community hyperedges are initialized
      // ! or not. In case of true some functions of streaming hypergraph
      // ! are not callable any more.
      bool _init_community_hyperedges;
  };

  class CommunityHyperedge {

    public:

      CommunityHyperedge() :
        _begin(0),
        _size(0),
        _weight(0),
        _hash(kEdgeHashSeed),
        _valid(false) { }

      CommunityHyperedge(const size_t begin,
                         const size_t size,
                         const HypernodeWeight weight) :
        _begin(begin),
        _size(size),
        _weight(weight),
        _hash(kEdgeHashSeed),
        _valid(true) { }

      void disable() {
        ASSERT(!isDisabled());
        _valid = false;
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
        _begin = begin;
      }

      // ! Returns the index of the first element in _incidence_array
      size_t firstInvalidEntry() const {
        return _begin + _size;
      }

      size_t size() const {
        return _size;
      }

      void setSize(size_t size) {
        _size = size;
      }

      void incrementSize() {
        ++_size;
      }

      void decrementSize() {
        ASSERT(_size > 0);
        --_size;
      }

      HyperedgeWeight weight() const {
        return _weight;
      }

      void setWeight(HyperedgeWeight weight) {
        _weight = weight;
      }

      size_t & hash() {
        return _hash;
      }

      size_t hash() const {
        return _hash;
      }

      bool operator== (const CommunityHyperedge& rhs) const {
        return _begin == rhs._begin && _size == rhs._size;
      }

      bool operator!= (const CommunityHyperedge& rhs) const {
        return !operator== (this, rhs);
      }

    private:
      // ! Index of the first element in _incidence_array
      size_t _begin;
      // ! Number of _incidence_array elements
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
                         typename ElementType::IDType>{   // reference

   public:
    using IDType = typename ElementType::IDType;

    HypergraphElementIterator() = default;

    HypergraphElementIterator(const HypergraphElementIterator& other) = default;
    HypergraphElementIterator& operator= (const HypergraphElementIterator& other) = default;

    HypergraphElementIterator(HypergraphElementIterator&& other) = default;
    HypergraphElementIterator& operator= (HypergraphElementIterator&& other) = default;

    ~HypergraphElementIterator() = default;

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
    HypergraphElementIterator& operator++ () {
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

    // ! Convenience function for range-based for-loops
    friend HypergraphElementIterator end<>(const std::pair<HypergraphElementIterator,
                                                           HypergraphElementIterator>& iter_pair);
    // ! Convenience function for range-based for-loops
    friend HypergraphElementIterator begin<>(const std::pair<HypergraphElementIterator,
                                                             HypergraphElementIterator>& iter_pair);

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
    const IDType _max_id = 0;
    // HypergraphElement the iterator currently points to
    const ElementType* _element = nullptr;
  };

  // ! Iterator to iterator over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterator over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the set of incident nets of a hypernode
  // ! or the set of pins of a hyperedge
  using IncidenceIterator = typename parallel::scalable_vector<HypernodeID>::const_iterator;
  using CommunityIterator = typename parallel::scalable_vector<PartitionID>::const_iterator;
  // ! Community Hyperedges
  using CommunityHyperedges = parallel::scalable_vector<parallel::scalable_vector<CommunityHyperedge>>;

 public:
  enum class UncontractionCase : uint8_t {
    CASE_1 = 0,
    CASE_2 = 1
  };

  /*!
  * A memento stores all information necessary to undo the contraction operation
  * of a vertex pair \f$(u,v)\f$.
  *
  * A contraction operations can increase the set \f$I(u)\f$ of nets incident
  * to \f$u\f$. This in turn leads to changes in _incidence_array. Therefore
  * the memento stores the initial starting index of \f$u\f$'s incident nets
  * as well as the old size, i.e. \f$|I(u)|\f$.
  *
  */
  struct Memento {

    Memento() :
      u(kInvalidHypernode),
      v(kInvalidHypernode),
      community_id(kInvalidPartition),
      one_pin_hes_begin(0),
      one_pin_hes_size(0),
      parallel_hes_begin(0),
      parallel_hes_size(0) { }

    Memento( HypernodeID representative, HypernodeID contraction_partner ) :
      u(representative),
      v(contraction_partner),
      community_id(kInvalidPartition),
      one_pin_hes_begin(0),
      one_pin_hes_size(0),
      parallel_hes_begin(0),
      parallel_hes_size(0) { }

    Memento( HypernodeID representative, HypernodeID contraction_partner, PartitionID community ) :
      u(representative),
      v(contraction_partner),
      community_id(community),
      one_pin_hes_begin(0),
      one_pin_hes_size(0),
      parallel_hes_begin(0),
      parallel_hes_size(0) { }

    // ! The representative hypernode that remains in the hypergraph
    HypernodeID u;
    // ! The contraction partner of u that is removed from the hypergraph after the contraction.
    HypernodeID v;
    // ! Community id of hypernodes
    PartitionID community_id;
    // ! start of removed single pin hyperedges
    int one_pin_hes_begin;
    // ! # removed single pin hyperedges
    int one_pin_hes_size;
    // ! start of removed parallel hyperedges
    int parallel_hes_begin;
    // ! # removed parallel hyperedges
    int parallel_hes_size;
  };

  explicit StreamingHypergraph(const int node, const bool remove_single_pin_community_nets = true) :
    _node(node),
    _arena(TBBNumaArena::instance().numa_task_arena(node)),
    _remove_single_pin_community_nets(remove_single_pin_community_nets),
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
    _community_hyperedge_ids(),
    _community_hyperedges(),
    _incident_nets_of_v(),
    _vertex_pin_count(),
    _pin_stream(),
    _hyperedge_stream(),
    _next_node_id(0),
    _hypernode_stream(),
    _incident_net_stream() {
    // Make sure constructor is called on corresponding numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Only allowed to allocate numa hypergraph on node" << _node << ", but it is"
      << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));
  }

  StreamingHypergraph(const StreamingHypergraph&) = delete;
  StreamingHypergraph& operator= (const StreamingHypergraph&) = delete;

  StreamingHypergraph(StreamingHypergraph&& other) :
    _node(other._node),
    _arena(other._arena),
    _remove_single_pin_community_nets(other._remove_single_pin_community_nets),
    _num_hypernodes(other._num_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_pins(other._num_pins),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _community_hyperedge_ids(std::move(other._community_hyperedge_ids)),
    _community_hyperedges(std::move(other._community_hyperedges)),
    _incident_nets_of_v(std::move(other._incident_nets_of_v)),
    _vertex_pin_count(std::move(other._vertex_pin_count)),
    _pin_stream(std::move(other._pin_stream)),
    _hyperedge_stream(std::move(other._hyperedge_stream)),
    _next_node_id(other._next_node_id.load()),
    _hypernode_stream(std::move(other._hypernode_stream)),
    _incident_net_stream(std::move(other._incident_net_stream)) { }

  StreamingHypergraph& operator= (StreamingHypergraph&&) = default;

  ~StreamingHypergraph() = default;

  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  HypernodeID originalNodeId(const HypernodeID u) const {
    return hypernode(u).originalNodeId();
  }

  HypernodeID originalEdgeId(const HyperedgeID e) const {
    return hyperedge(e).originalEdgeId();
  }

  HypernodeID communityNodeId(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityNodeId();
  }

  size_t numCommunitiesOfHyperedge(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    return _community_hyperedges[local_id].size();
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return std::make_pair(incident_nets(u).cbegin() + hypernode(u).invalidIncidentNets(),
                          incident_nets(u).cend());
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  // ! Note, in contrast to first iterator, this iterator skips all invalidated
  // ! community hyperedges in incident nets of hypernode u.
  std::pair<IncidenceIterator, IncidenceIterator> validIncidentEdges(const HypernodeID u, const PartitionID) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    ASSERT(hypernode(u).invalidCommunityNets() <= incident_nets(u).size());
    return std::make_pair(incident_nets(u).cbegin(),
                          incident_nets(u).cbegin() + hypernode(u).invalidCommunityNets());
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  // ! Note, in contrast to first iterator, this iterator skips all single-pin and invalidated
  // ! community hyperedges in incident nets of hypernode u.
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u, const PartitionID) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    ASSERT(hypernode(u).singlePinCommunityNets() <= hypernode(u).invalidCommunityNets());
    ASSERT(hypernode(u).invalidCommunityNets() <= incident_nets(u).size());
    return std::make_pair(incident_nets(u).cbegin() + hypernode(u).singlePinCommunityNets(),
                          incident_nets(u).cbegin() + hypernode(u).invalidCommunityNets());
  }

  // ! Returns a const reference to the incident net vector of hypernode u
  const parallel::scalable_vector<HyperedgeID>& incidentNets(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return incident_nets(u);
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(!hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are initialized");
    return std::make_pair(_incidence_array.cbegin() + hyperedge(e).firstEntry(),
                          _incidence_array.cbegin() + hyperedge(e).firstInvalidEntry());
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e in a community.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    const CommunityHyperedge& community_he = community_hyperedge(e, community_id);
    return std::make_pair(_incidence_array.cbegin() + community_he.firstEntry(),
                          _incidence_array.cbegin() + community_he.firstInvalidEntry());
  }

  std::pair<CommunityIterator, CommunityIterator> communities(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(local_id < _community_hyperedge_ids.size());
    return std::make_pair(_community_hyperedge_ids[local_id].cbegin(),
                          _community_hyperedge_ids[local_id].cend());
  }

  /*!
   * Returns a for-each iterator-pair to loop over the set of all hypernodes.
   * Since the iterator just skips over disabled hypernodes, iteration always
   * iterates over all _num_hypernodes hypernodes.
   */
  std::pair<HypernodeIterator, HypernodeIterator> nodes() const {
    HypernodeID start = get_global_node_id(0);
    HypernodeID end = get_global_node_id(_num_hypernodes);
    return std::make_pair(HypernodeIterator(_hypernodes.data(), start, end),
                          HypernodeIterator((_hypernodes.data() + _num_hypernodes), end, end));
  }

  /*!
   * Returns a for-each iterator-pair to loop over the set of all hyperedges.
   * Since the iterator just skips over disabled hyperedges, iteration always
   * iterates over all _num_hyperedges hyperedges.
   */
  std::pair<HyperedgeIterator, HyperedgeIterator> edges() const {
    HyperedgeID start = get_global_edge_id(0UL);
    HyperedgeID end = get_global_edge_id(_num_hyperedges);
    return std::make_pair(HyperedgeIterator(_hyperedges.data(), start, end),
                          HyperedgeIterator((_hyperedges.data() + _num_hyperedges), end, end));
  }


  size_t vertexPinCount(const HypernodeID hn) const {
    ASSERT(hn < _vertex_pin_count.size());
    return _vertex_pin_count[hn];
  }

  void contract(const HypernodeID u, const HypernodeID v,
                const HyperedgeID e, Self& hypergraph_of_u) {
    ASSERT(get_numa_node_of_hyperedge(e) == _node);
    using std::swap;
    const HypernodeID pins_begin = hyperedge(e).firstEntry();
    const HypernodeID pins_end = hyperedge(e).firstInvalidEntry();
    HypernodeID slot_of_u = pins_end - 1;
    HypernodeID last_pin_slot = pins_end - 1;

    // Swap contraction partner v to end of the pin list of he and find slot of
    // representative u.
    for (HypernodeID pin_iter = pins_begin; pin_iter != last_pin_slot; ++pin_iter) {
      const HypernodeID pin = _incidence_array[pin_iter];
      if (pin == v) {
        swap(_incidence_array[pin_iter], _incidence_array[last_pin_slot]);
        --pin_iter;
      } else if (pin == u) {
        slot_of_u = pin_iter;
      }
    }

    ASSERT(_incidence_array[last_pin_slot] == v, "v is not last entry in incidence array!");

    if ( slot_of_u != last_pin_slot ) {
      // Case 1:
      // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
      // can just cut off the last entry in the edge array of e that now contains v.
      DBG << V(e) << ": Case 1";
      hyperedge(e).hash() -= kahypar::math::hash(v);
      hyperedge(e).decrementSize();
      // TODO(heuer): Update pin count in part
    } else {
      // Case 2:
      // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
      // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
      DBG << V(e) << ": Case 2";
      hyperedge(e).hash() -= kahypar::math::hash(v);
      hyperedge(e).hash() += kahypar::math::hash(u);
      connectHyperedgeToRepresentative(e, u, hypergraph_of_u);
    }
  }

  void contract(const HypernodeID u, const HypernodeID v,
                const HyperedgeID e, const PartitionID community_id,
                Self& hypergraph_of_u) {
    ASSERT(get_numa_node_of_hyperedge(e) == _node);

    using std::swap;
    CommunityHyperedge& community_he = community_hyperedge(e, community_id);
    const HypernodeID pins_begin = community_he.firstEntry();
    const HypernodeID pins_end = community_he.firstInvalidEntry();
    HypernodeID slot_of_u = pins_end - 1;
    HypernodeID last_pin_slot = pins_end - 1;

    // Swap contraction partner v to end of the pin list of he and find slot of
    // representative u.
    for (HypernodeID pin_iter = pins_begin; pin_iter != last_pin_slot; ++pin_iter) {
      const HypernodeID pin = _incidence_array[pin_iter];
      if (pin == v) {
        swap(_incidence_array[pin_iter], _incidence_array[last_pin_slot]);
        --pin_iter;
      } else if (pin == u) {
        slot_of_u = pin_iter;
      }
    }

    ASSERT(_incidence_array[last_pin_slot] == v, "v is not last entry in incidence array!");

    if ( slot_of_u != last_pin_slot ) {
      // Case 1:
      // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
      // can just cut off the last entry in the edge array of e that now contains v.
      DBG << V(e) << ": Case 1";
      community_he.hash() -= kahypar::math::hash(v);
      community_he.decrementSize();
      // TODO(heuer): Update pin count in part
    } else {
      // Case 2:
      // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
      // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
      DBG << V(e) << ": Case 2";
      community_he.hash() -= kahypar::math::hash(v);
      community_he.hash() += kahypar::math::hash(u);
      _incidence_array[community_he.firstInvalidEntry() - 1] = u;
      connectHyperedgeToRepresentative(e, u, community_he, hypergraph_of_u);
    }
  }

  bool uncontract(const HypernodeID u, const HypernodeID v,
                  const HyperedgeID e, const size_t incident_nets_pos,
                  std::vector<Self>& hypergraphs) {
    using std::swap;

    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if ( containsIncidentNet(e) ) {
      // ... then we have to do some kind of restore operation.
      if ( get_uncontraction_case(e, hyperedge(e).size(), v) == UncontractionCase::CASE_1 ) {
        // hyperedge(he + 1) always exists because of sentinel
        // Undo case 1 operation (i.e. Pin v was just cut off by decreasing size of HE e)
        DBG << V(e) << " -> case 1";
        hyperedge(e).incrementSize();

        // TODO(heuer): Increment pin count in part
        return false;
      } else {
        // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
        // Set incidence entry containing u for this HE e back to v, because this slot was used
        // to store the new edge to representative u during contraction as u was not a pin of e.
        DBG << V(e) << " -> case 2";
        Self& hypergraph_of_u = hypergraph_of_vertex(u, hypergraphs);
        ASSERT(e == hypergraph_of_u.incident_nets(u)[incident_nets_pos]);

        size_t incident_nets_end = hypergraph_of_u.incident_nets(u).size();
        swap(hypergraph_of_u.incident_nets(u)[incident_nets_pos],
             hypergraph_of_u.incident_nets(u)[incident_nets_end - 1]);
        hypergraph_of_u.incident_nets(u).pop_back();

        DBG << "resetting reused Pinslot of HE" << e << "from" << u << "to" << v;
        resetReusedPinSlotToOriginalValue(e, u, v);

        // TODO(heuer): Increment pin count in part
        return true;
      }
    }
    return false;
  }

  #ifndef NDEBUG
  /**
   * This function uncontracts a disabled parallel hyperedge. Hyperedge e and
   * representative are parallel and e is disabled (whereas representative is enabled).
   */
  bool uncontract(const HypernodeID u, const HypernodeID v,
                  const HyperedgeID e, const HyperedgeID representative,
                  const size_t incident_nets_pos,
                  std::vector<Self>& hypergraphs) {
    using std::swap;
    ASSERT(hyperedge(e).isDisabled(), "Hyperedge" << e << "is enabled");
    const Self& hypergraph_of_rep = hypergraph_of_hyperedge(representative, hypergraphs);
    ASSERT(hypergraph_of_rep.edgeIsEnabled(representative), "Hyperedge" << representative << "is disabled");

    if ( containsIncidentNet(e) ) {
      ASSERT(hypergraph_of_rep.edgeSize(representative) <=
             hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry(),
             V(hypergraph_of_rep.edgeSize(representative))
             << V((hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry())));

      size_t edge_size = hypergraph_of_rep.edgeSize(representative);
      if ( get_uncontraction_case(e, edge_size, v) == UncontractionCase::CASE_2 ) {
        // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
        // Set incidence entry containing u for this HE e back to v, because this slot was used
        // to store the new edge to representative u during contraction as u was not a pin of e.
        DBG << V(e) << " -> case 2";
        Self& hypergraph_of_u = hypergraph_of_vertex(u, hypergraphs);
        ASSERT(e == hypergraph_of_u.incident_nets(u)[incident_nets_pos]);

        size_t invalid_incident_nets = hypergraph_of_u.hypernode(u).invalidIncidentNets();
        ASSERT(incident_nets_pos < invalid_incident_nets);
        swap(hypergraph_of_u.incident_nets(u)[incident_nets_pos],
             hypergraph_of_u.incident_nets(u)[invalid_incident_nets - 1]);
        hypergraph_of_u.hypernode(u).decrementInvalidIncidentNets();

        size_t incident_nets_end = hypergraph_of_u.incident_nets(u).size();
        swap(hypergraph_of_u.incident_nets(u)[invalid_incident_nets - 1],
             hypergraph_of_u.incident_nets(u)[incident_nets_end - 1]);
        hypergraph_of_u.incident_nets(u).pop_back();

        DBG << "resetting reused Pinslot of HE" << e << "from" << u << "to" << v;
        resetReusedPinSlotToOriginalValue(e, edge_size, u, v);

        // TODO(heuer): Increment pin count in part
        return true;
      }
    }

    return false;
  }
  #endif

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).weight();
  }

  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    hypernode(u).setWeight(weight);
  }

  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).weight();
  }

  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).weight();
  }

  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    hyperedge(e).setWeight(weight);
  }

  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).setWeight(weight);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return incident_nets(u).size() - hypernode(u).invalidIncidentNets();
  }

  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if ( hyperedge(e).isInitCommunityHyperedges() ) {
      HypernodeID size = 0;
      HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
      ASSERT(local_id < _community_hyperedges.size());
      for ( const CommunityHyperedge& community_he : _community_hyperedges[local_id] ) {
        size += community_he.size();
      }
      return size;
    } else {
      return hyperedge(e).size();
    }
  }

  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).size();
  }

  size_t edgeHash(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if ( hyperedge(e).isInitCommunityHyperedges() ) {
      size_t hash = 0;
      HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
      ASSERT(local_id < _community_hyperedges.size());
      for ( const CommunityHyperedge& community_he : _community_hyperedges[local_id] ) {
        hash += community_he.hash();
      }
      return hash;
    } else {
      return hyperedge(e).hash();
    }
  }

  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).hash();
  }

  PartitionID communityID(const HypernodeID u) const {
    //ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityID();
  }

  void setPartInfo(const HypernodeID u, const PartitionID id) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    hypernode(u).setPartID(id);
  }

  PartitionID partID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).partID();
  }

  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  bool edgeIsEnabled(const HyperedgeID e) const {
    return !hyperedge(e).isDisabled();
  }

  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  void enableHyperedge(const HyperedgeID e) {
    hyperedge(e).enable();
  }

  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
  }

  void disableHyperedge(const HyperedgeID e, const PartitionID community_id) {
    community_hyperedge(e, community_id).disable();
  }

  HypernodeID streamHypernode(HypernodeID original_id, HypernodeWeight weight) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node);
    HypernodeID node_id = get_global_node_id();
    _hypernode_stream.stream(node_id, original_id, weight);
    return node_id;
  }

  void streamHyperedge(const parallel::scalable_vector<HypernodeID>& hyperedge,
                       const HyperedgeID original_id,
                       const HyperedgeWeight& weight) {
    int cpu_id = sched_getcpu();
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(cpu_id) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was CPU" << cpu_id
       << "is on node" << HardwareTopology::instance().numa_node_of_cpu(cpu_id));
    _hyperedge_stream.stream(_pin_stream.size(cpu_id), hyperedge.size(), original_id, weight);
    for ( const HypernodeID& pin : hyperedge  ) {
      _pin_stream.stream(pin);
    }
  }

  void streamIncidentNet(const HypernodeID hn, const HyperedgeID he) {
    ASSERT(get_numa_node_of_vertex(hn) == _node);
    _incident_net_stream.stream(hn, he);
  }

  void streamCommunityID(const HypernodeID hn, const PartitionID community_id) {
    hypernode(hn).setCommunityID(community_id);
  }

  void initializeHyperedges(const HypernodeID num_hypernodes) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    // Copy streamed data into global vectors
     HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    _incidence_array = _pin_stream.copy(_arena);
    _hyperedges = _hyperedge_stream.copy(_arena);
    _num_pins = _incidence_array.size();
    _num_hyperedges = _hyperedges.size();
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("copy_incidencce_array_and_he", "Copy Incidence Array and HEs",
      "initialize_hyperedges", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());


    kahypar::ds::FastResetFlagArray<> tmp_incidence_nets_of_v(_num_hyperedges);
    _incident_nets_of_v = std::move(tmp_incidence_nets_of_v);

    // Update start position of each hyperedge to correct one in global incidence array
    // Note, start positions are stored relative to the local buffer there are streamed into.
    // However, memcpy does hyperedges invalidates those local positions.
    start = std::chrono::high_resolution_clock::now();
    tbb::task_group group;
    _arena.execute([&] {
      for ( size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
        group.run([&, cpu_id] {
          size_t start = _hyperedge_stream.prefix_sum(cpu_id);
          size_t end = start + _hyperedge_stream.size(cpu_id);
          size_t delta = _pin_stream.prefix_sum(cpu_id);
          ASSERT(end <= _hyperedges.size());
          for ( size_t pos = start; pos < end; ++pos ) {
            _hyperedges[pos].setFirstEntry( _hyperedges[pos].firstEntry() + delta );
          }
        });
      }
    });
    group.wait();
    // Emplace Back Sentinel
    _hyperedges.emplace_back(_incidence_array.size(), 0, 0UL, 0);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("update_start_position", "Update Start Positions of HEs",
      "initialize_hyperedges", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());


    ASSERT([&] {
      for ( size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id ) {
        size_t start = _hyperedge_stream.prefix_sum(cpu_id);
        for ( size_t idx = 0; idx < _hyperedge_stream.size(cpu_id); ++idx ) {
          const Hyperedge& stream_he = _hyperedge_stream.value(cpu_id, idx);
          const Hyperedge& he = _hyperedges[idx + start];
          if ( stream_he.size() != he.size() ) {
            LOG << "Unequal size";
            return false;
          }
          size_t global_pos = he.firstEntry();
          for ( size_t pos = stream_he.firstEntry(); pos < stream_he.firstEntry() + stream_he.size(); ++pos ) {
            if ( _incidence_array[global_pos] != _pin_stream.value(cpu_id, pos) ) {
              LOG << "Pins in stream and global incidence array are not equal"
                  << V(_incidence_array[global_pos]) << V(_pin_stream.value(cpu_id, pos));
              return false;
            }
            ++global_pos;
          }
        }
      }
      return true;
    }(), "Failed to copy buffer to hypergraph");

    // Compute how many times a hypernode occurs on this node
    // as pin. Will be later important to compute node assignment.
    // TODO(heuer): Think how to parallelize this
    start = std::chrono::high_resolution_clock::now();
    _vertex_pin_count.resize(num_hypernodes);
    for ( size_t i = 0; i < _incidence_array.size(); ++i ) {
      const HypernodeID& pin = _incidence_array[i];
      ASSERT(pin < _vertex_pin_count.size());
      _vertex_pin_count[pin]++;
    }
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_vertex_pin_count", "Compute Vertex Pin Counts",
      "initialize_hyperedges", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());


    _pin_stream.clear();
    _hyperedge_stream.clear();
  }

  void initializeHypernodes(std::vector<Self>& hypergraphs,
                            const std::vector<HypernodeID>& node_mapping) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _hypernodes = _hypernode_stream.copy(_arena);
    _num_hypernodes = _hypernodes.size();

    // Sort hypernodes in increasing order of their node id and ...
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_sort(_hypernodes.begin(), _hypernodes.end(),
          [&](const Hypernode& lhs, const Hypernode& rhs) {
          return lhs.nodeId() < rhs.nodeId();
        });
      });
    });


    // ... parallel to that remap node ids in incidence array to new node ids
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, _incidence_array.size()),
          [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t pos = range.begin(); pos < range.end(); ++pos ) {
            HypernodeID pin = _incidence_array[pos];
            ASSERT(pin < node_mapping.size());
            _incidence_array[pos] = node_mapping[pin];
          }
        });
      });
    });
    group.wait();
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("sort_and_remap_node_ids", "Sort and Remap Nodes",
      "initialize_numa_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());

    // Compute Total Hypergraph Weight
    start = std::chrono::high_resolution_clock::now();
    updateTotalWeight();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_total_weight", "Compute Total Weight",
      "initialize_numa_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());

    ASSERT([&]{
      for ( size_t i = 0; i < _hypernodes.size(); ++i ) {
        const Hypernode& hn = _hypernodes[i];
        const int node = get_numa_node_of_vertex(hn.nodeId());
        const HypernodeID local_id = get_local_node_id_of_vertex(hn.nodeId());
        if ( node != _node ) {
          LOG << "Hypernode" << hn.nodeId() << "should be on numa node" << node
              << ", but is on numa node" << _node;
          return false;
        } else if ( i != local_id ) {
          LOG << "Hypernode" << local_id << "should have local id" << i;
          return false;
        }

        if ( i > 0 ) {
          const Hypernode& previous_hn = _hypernodes[i-1];
          const HypernodeID previous_local_id = get_local_node_id_of_vertex(previous_hn.nodeId());
          if ( local_id <= previous_local_id ) {
            LOG << "Local node ids not sorted" << V(previous_local_id) << V(local_id);
            return false;
          } else if ( previous_local_id + 1 != local_id ) {
            LOG << "Local node ids not consecutive" << V(previous_local_id) << V(local_id);
            return false;
          }
        }
      }
      return true;
    }(), "Initialization of hypernodes failed");

    // Stream incident nets
    start = std::chrono::high_resolution_clock::now();
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, _hyperedges.size()),
          [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t pos = range.begin(); pos < range.end(); ++pos ) {
            Hyperedge& he = _hyperedges[pos];
            const HyperedgeID he_id = get_global_edge_id(pos);
            for ( size_t incidence_array_pos = he.firstEntry();
                  incidence_array_pos < he.firstEntry() + he.size();
                  ++incidence_array_pos ) {
              const HypernodeID pin = _incidence_array[incidence_array_pos];
              hypergraph_of_vertex(pin, hypergraphs).streamIncidentNet(pin, he_id);
              // Initialize edge hash
              he.hash() += kahypar::math::hash(pin);
            }
          }
        });
      });
    });
    group.wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("stream_incident_nets", "Stream Incident Nets",
      "initialize_numa_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());

    _hypernode_stream.clear();
  }

  void updateTotalWeight() {
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 0,
          [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
            HypernodeWeight weight = init;
            for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
              weight += this->_hypernodes[hn].weight();
            }
            return weight;
          },
          std::plus<HypernodeWeight>());
      });
    });
    group.wait();
  }

  void initializeIncidentNets() {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was on node"
        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _incident_nets.resize(_hypernodes.size());
    _incident_net_stream.copy(_arena, _incident_nets, [&](const HypernodeID& u) {
      ASSERT(get_numa_node_of_vertex(u) == _node);
      return get_local_node_id_of_vertex(u);
    });
    _incident_net_stream.clear();
  }

  void initializeCommunityHyperedges(const std::vector<Self>& hypergraphs) {

    auto add_community_hyperedge = [&](const HyperedgeID he,
                                       const PartitionID community_id,
                                       const size_t start,
                                       const size_t end,
                                       const HyperedgeWeight weight) {
      if ( community_id != kInvalidPartition ) {
        ASSERT(he < this->_community_hyperedges.size());
        ASSERT(start < end);
        this->_community_hyperedge_ids[he].emplace_back(community_id);
        this->_community_hyperedges[he].emplace_back(start, end - start, weight);

        // Compute community hyperedge hash
        for ( size_t pos = start; pos < end; ++pos ) {
          this->_community_hyperedges[he].back().hash() +=
            kahypar::math::hash(this->_incidence_array[pos]);
        }
      }
    };

    // Add community hyperedges
    _community_hyperedge_ids.resize(_num_hyperedges);
    _community_hyperedges.resize(_num_hyperedges);
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, this->_num_hyperedges),
        [&](const tbb::blocked_range<size_t>& range) {
          for ( HyperedgeID he = range.begin(); he < range.end(); ++he ) {
            Hyperedge& e = this->_hyperedges[he];
            if ( !e.isDisabled() ) {
              // Sort pins of hyperedge in increasing order of their community ids
              size_t incidence_array_start = e.firstEntry();
              size_t incidence_array_end = e.firstInvalidEntry();
              std::sort(this->_incidence_array.begin() + incidence_array_start,
                        this->_incidence_array.begin() + incidence_array_end,
                        [&](const HypernodeID& lhs, const HypernodeID& rhs) {
                          return hypergraph_of_vertex(lhs, hypergraphs).communityID(lhs) <
                                hypergraph_of_vertex(rhs, hypergraphs).communityID(rhs);
                        });

              // Add community hyperedges for each consecutive range of pins with
              // the same community id
              size_t last_community_start = incidence_array_start;
              PartitionID last_community_id = kInvalidPartition;
              for ( size_t incidence_array_pos = incidence_array_start;
                    incidence_array_pos < incidence_array_end;
                    ++incidence_array_pos ) {
                const HypernodeID pin = this->_incidence_array[incidence_array_pos];
                const PartitionID community_id = hypergraph_of_vertex(pin, hypergraphs).communityID(pin);
                if ( community_id != last_community_id ) {
                  add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_pos, e.weight());
                  last_community_start = incidence_array_pos;
                  last_community_id = community_id;
                }
              }
              add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_end, e.weight());
              e.initializeCommunityHyperedges();

              ASSERT([&] {
                if ( e.firstEntry() != _community_hyperedges[he][0].firstEntry() ) {
                  return false;
                }
                for ( size_t i = 1; i < _community_hyperedges[he].size(); ++i ) {
                  if ( _community_hyperedges[he][i - 1].firstInvalidEntry() !=
                      _community_hyperedges[he][i].firstEntry() ) {
                    return false;
                  }
                }
                if ( e.firstInvalidEntry() != _community_hyperedges[he].back().firstInvalidEntry() ) {
                  return false;
                }
                return true;
              }(), "Initialization of community hyperedges failed!");
            }
          }
        });
      });
    });
    group.wait();
  }

  void initializeCommunityHypernodes(const std::vector<Self>& hypergraphs) {
    // Partition incident nets of each hypernode such that in first part are all hyperedges
    // which contains at least to pins and in second part are all single-pin community hyperedges
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, this->_num_hypernodes),
        [&](const tbb::blocked_range<HypernodeID>& range) {
          for ( HypernodeID v = range.begin(); v < range.end(); ++v ) {
            HypernodeID hn = get_global_node_id(v);
            if ( !this->hypernode(hn).isDisabled() ) {
              PartitionID community_id = this->communityID(hn);
              int single_pin_community_nets = 0;
              if ( _remove_single_pin_community_nets ) {
                for ( int i = 0; i < (int) _incident_nets[v].size(); ++i ) {
                  HyperedgeID he = this->_incident_nets[v][i];
                  if ( hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) == 1 ) {
                    std::swap(this->_incident_nets[v][i],
                              this->_incident_nets[v][single_pin_community_nets]);
                    ++single_pin_community_nets;
                  }
                }
              }
              this->hypernode(hn).setSinglePinCommunityNets(single_pin_community_nets);
              this->hypernode(hn).setInvalidCommunityNets(_incident_nets[v].size());

              ASSERT([&] {
                if ( _remove_single_pin_community_nets ) {
                  size_t single_pin_community_hyperedges = this->hypernode(hn).singlePinCommunityNets();
                  for ( size_t i = 0; i < single_pin_community_hyperedges; ++i ) {
                    const HyperedgeID he = incident_nets(hn)[i];
                    if ( hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) > 1 ) {
                      LOG << "Hyperedge" << he << "is a non single-pin commnunity hyperedge";
                      return false;
                    }
                  }
                }
                return true;
              }(), "There non single-pin community hyperedges in single-pin part of incident nets");

              ASSERT([&] {
                if ( _remove_single_pin_community_nets ) {
                  size_t single_pin_community_hyperedges = this->hypernode(hn).singlePinCommunityNets();
                  for ( size_t i = single_pin_community_hyperedges; i < incident_nets(hn).size(); ++i ) {
                    const HyperedgeID he = incident_nets(hn)[i];
                    if ( hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) <= 1 ) {
                      LOG << "Hyperedge" << he << "is a single-pin commnunity hyperedge";
                      return false;
                    }
                  }
                }
                return true;
              }(), "There single-pin community hyperedges in non-single-pin part of incident nets");
            }
          }
        });
      });
    });
    group.wait();
  }

  void resetCommunityHyperedges(const std::vector<Memento>& mementos,
                                const HypernodeID num_hypernodes,
                                const std::vector<Self>& hypergraphs) {
    // All disabled hypernodes have to follow a specific order in invalid part of the incidence array
    // such that they can be successfully uncontracted. They have be sorted in decreasing order of their
    // contraction. In order to realize this we compute the contraction index of a hypernode inside the
    // contraction history and use it later for sorting them.
    parallel::scalable_vector<HypernodeID> contraction_index(num_hypernodes, std::numeric_limits<HypernodeID>::max());
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, mementos.size()),
        [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t i = range.begin(); i < range.end(); ++i ) {
            const HypernodeID v = mementos[i].v;
            ASSERT(hypergraph_of_vertex(v, hypergraphs).originalNodeId(v) < num_hypernodes);
            ASSERT(contraction_index[hypergraph_of_vertex(v, hypergraphs).originalNodeId(v)] == std::numeric_limits<HypernodeID>::max(),
              "Hypernode" << v << "occurs more than once as contraction partner in hierarchy");
            contraction_index[hypergraph_of_vertex(v, hypergraphs).originalNodeId(v)] = i;
          }
        });
      });
    });
    group.wait();

    // The incidence array of a hyperedge is constructed as follows: The first part consists
    // of all enabled pins and the remainder of all invalid pins. The invalid pins in the
    // remainder are sorted in decreasing order of their contraction index.
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, _num_hyperedges),
        [&](const tbb::blocked_range<HyperedgeID>& range) {
          for ( HyperedgeID he = range.begin(); he < range.end(); ++he ) {
            Hyperedge& e = this->_hyperedges[he];
            if ( e.isInitCommunityHyperedges() ) {
              int64_t first_entry = e.firstEntry();
              int64_t last_entry = e.firstInvalidEntry();
              // Sort pins of hyperedge in decreasing order of their contraction index
              std::sort(this->_incidence_array.begin() + first_entry,
                        this->_incidence_array.begin() + last_entry,
                        [&](const HypernodeID& lhs, const HypernodeID& rhs) {
                          HypernodeID original_lhs = hypergraph_of_vertex(lhs, hypergraphs).originalNodeId(lhs);
                          HypernodeID original_rhs = hypergraph_of_vertex(rhs, hypergraphs).originalNodeId(rhs);
                          return contraction_index[original_lhs] > contraction_index[original_rhs];
                        });

              // Count number of enabled hypernodes
              --last_entry;
              ASSERT(first_entry <= last_entry);
              for ( ; last_entry >= first_entry; --last_entry ) {
                const HypernodeID pin = this->_incidence_array[last_entry];
                if ( !hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin) ) {
                  e.decrementSize();
                }
              }

              e.deinitializeCommunityHyperedges();

              ASSERT([&] {
                for ( size_t i = e.firstEntry(); i < e.firstInvalidEntry(); ++i ) {
                  const HypernodeID& pin = this->_incidence_array[i];
                  if ( !hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin) ) {
                    LOG << "Hypernode" << pin << "is disabled";
                    return false;
                  }
                }
                return true;
              }(), "There are disabled hypernodes in valid part of hyperedge");

              ASSERT([&] {
                for ( size_t i = e.firstInvalidEntry(); i < _hyperedges[he + 1].firstEntry(); ++i ) {
                  const HypernodeID& pin = this->_incidence_array[i];
                  if ( hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin) ) {
                    LOG << "Hypernode" << pin << "is enabled";
                    return false;
                  }
                }
                return true;
              }(), "There are enabled hypernodes in invalid part of hyperedge");
            }
          }
        });
      });
    });
    group.wait();

    CommunityHyperedges community_hyperedges;
    _community_hyperedges = std::move(community_hyperedges);
  }

  void resetPinsToOriginalNodeIds(const std::vector<Self>& hypergraphs) {
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<size_t>(0UL, this->_num_pins),
        [&](const tbb::blocked_range<size_t>& range) {
          for ( size_t i = range.begin(); i < range.end(); ++i ) {
            HypernodeID pin = this->_incidence_array[i];
            this->_incidence_array[i] = hypergraph_of_vertex(pin, hypergraphs).hypernode(pin).originalNodeId();
          }
        });
      });
    });
    group.wait();
  }

  /**
   * During parallel hyperedge removal it can happen that some community hypergraph pruner
   * do not detect that a hyperedge becomes parallel with an other. We fix this in a postprocessing step
   * after coarsening separate. However, the corresponding pins of these hyperedges might contain
   * the disabled hyperedge, which is fixed here.
   */
  void removeDisabledHyperedgesFromIncidentNets(const std::vector<Self>& hypergraphs) {
    tbb::task_group group;
    _arena.execute([&] {
      group.run([&] {
        tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, this->_num_hypernodes),
        [&](const tbb::blocked_range<HypernodeID>& range) {
          for ( HypernodeID id = range.begin(); id < range.end(); ++id ) {
            const HypernodeID hn = get_global_node_id(id);
            removeDisabledHyperedgesFromIncidentNets(hn, hypergraphs);
          }
        });
      });
    });
    group.wait();
  }

  // ! Only for assertion
  bool verify_incident_nets_of_hypergraph(const std::vector<Self>& hypergraphs) const {
    for ( size_t pos = 0; pos < _incident_nets.size(); ++pos ) {
      const HypernodeID& hn = _hypernodes[pos].nodeId();
      for ( const HyperedgeID& he : _incident_nets[pos] ) {
        const Self& hypergraph_of_he = hypergraph_of_hyperedge(he, hypergraphs);
        const HyperedgeID local_edge_id = get_local_edge_id_of_hyperedge(he);
        ASSERT(local_edge_id < hypergraph_of_he._hyperedges.size());
        const Hyperedge& e = hypergraph_of_he._hyperedges[local_edge_id];
        const auto first = hypergraph_of_he._incidence_array.begin() + e.firstEntry();
        const auto last = first + e.size();
        if ( std::find(first, last, hn) == last ) {
          LOG << "Hypernode" << hn << "not part of hyperedge" << he << "on numa node" << get_numa_node_of_hyperedge(he);
          return false;
        }
      }
    }
    return true;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_global_edge_id(const int node, const size_t edge_pos) {
    return ( ( (HyperedgeID) node ) << NUMA_NODE_INDENTIFIER ) | edge_pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID u) {
    return (int) (u >> NUMA_NODE_INDENTIFIER);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_hyperedge(const HyperedgeID e) {
    return (int) (e >> NUMA_NODE_INDENTIFIER);
  }

  void printHyperedgeInfo(const HyperedgeID e) const {
    LOG << "Hyperedge:" << e;
    LOG << "Original Size:" << (hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry());
    if ( edgeIsEnabled(e) ) {
      LOG << "Current Size:" << hyperedge(e).size();
    }
    for ( size_t pos = hyperedge(e).firstEntry(); pos < hyperedge(e + 1).firstEntry(); ++pos ) {
      std::cout << _incidence_array[pos] << " ";
    }
    std::cout << std::endl;
  }

 private:

  template <typename NodeType_,
            typename EdgeType_,
            typename NodeWeightType_,
            typename EdgeWeightType_,
            typename PartitionID_,
            typename HardwareTopology_,
            typename TBBNumaArena_>
  friend class Hypergraph;

  /*!
   * Connect hyperedge e to representative hypernode u.
   */
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHyperedgeToRepresentative(const HyperedgeID e,
                                                                        const HypernodeID u,
                                                                        Self& hypergraph_of_u) {
    ASSERT(hypergraph_of_u.nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(_node == get_numa_node_of_hyperedge(e));
    // Hyperedge e does not contain u. Therefore we use the entry of v (i.e. the last entry
    // -- this is ensured by the contract method) in e's edge array to store the information
    // that u is now connected to e and add the edge (u,e) to indicate this conection also from
    // the hypernode's point of view.
    _incidence_array[hyperedge(e).firstInvalidEntry() - 1] = u;
    hypergraph_of_u.incident_nets(u).push_back(e);
  }

  /*!
   * Connect hyperedge e to representative hypernode u.
   */
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHyperedgeToRepresentative(const HyperedgeID e,
                                                                        const HypernodeID u,
                                                                        CommunityHyperedge& community_he,
                                                                        Self& hypergraph_of_u) {
    ASSERT(hypergraph_of_u.nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(_node == get_numa_node_of_hyperedge(e));

    auto& incident_nets_of_u = hypergraph_of_u.incident_nets(u);
    incident_nets_of_u.push_back(e);
    size_t he_pos = incident_nets_of_u.size() - 1;

    // If community hyperedge is not disabled we swap it to valid part of incident nets
    if ( !community_he.isDisabled() ) {
      size_t invalid_community_nets = hypergraph_of_u.hypernode(u).invalidCommunityNets();
      std::swap(incident_nets_of_u[invalid_community_nets], incident_nets_of_u[he_pos]);
      he_pos = invalid_community_nets;
      hypergraph_of_u.hypernode(u).incrementInvalidCommunityNets();

      // If community hyperedge is a single-pin community hyperedge we swap it to single-pin
      // part of incident nets of u
      if ( _remove_single_pin_community_nets && community_he.size() == 1 ) {
        size_t single_pin_community_nets = hypergraph_of_u.hypernode(u).singlePinCommunityNets();
        std::swap(incident_nets_of_u[single_pin_community_nets], incident_nets_of_u[he_pos]);
        hypergraph_of_u.hypernode(u).incrementSinglePinCommunityNets();
      }
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void markAllIncidentNetsOf(const HypernodeID v,
                                                             std::vector<Self>& hypergraphs) {
    // Reset incident nets on each node
    for ( size_t node = 0; node < hypergraphs.size(); ++node ) {
      hypergraphs[node]._incident_nets_of_v.reset();
    }

    for ( const HyperedgeID& he : incident_nets(v) ) {
      Self& hypergraph_of_he = hypergraph_of_hyperedge(he, hypergraphs);
      HyperedgeID local_he_id = get_local_edge_id_of_hyperedge(he);
      ASSERT(local_he_id < hypergraph_of_he._num_hyperedges);
      hypergraph_of_he._incident_nets_of_v.set(local_he_id, true);
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool containsIncidentNet(const HyperedgeID e) const {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of this numa node");
    ASSERT(local_id <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _incident_nets_of_v[local_id];
  }

  /*!
   * Resets the pin slot containing u back to contain v.
   * If hyperedge he only contained v before the contraction, then the array entry of
   * v in the incidence structure of he is used to store u after the contraction.
   * This method undoes this operation.
   */
  void resetReusedPinSlotToOriginalValue(const HyperedgeID he,
                                         const HypernodeID u,
                                         const HypernodeID v) {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge" << he << "is disabled");

    int64_t pins_start = hyperedge(he).firstEntry();
    int64_t pins_end = hyperedge(he).firstInvalidEntry();
    int64_t slot_of_u = pins_end;

    for ( int64_t pos = pins_end - 1; pos >= pins_start; --pos ) {
      if ( _incidence_array[pos] == u ) {
        slot_of_u = pos;
        break;
      }
    }
    ASSERT(slot_of_u != pins_end, "Hypernode" << u << "not found in hyperedge" << he);
    ASSERT(_incidence_array[slot_of_u] == u, "Hypernode" << u << "not found in hyperedge" << he);
    _incidence_array[slot_of_u] = v;
  }

  void resetReusedPinSlotToOriginalValue(const HyperedgeID he,
                                         const size_t size,
                                         const HypernodeID u,
                                         const HypernodeID v) {
    ASSERT(hyperedge(he).isDisabled(), "Hyperedge" << he << "is enabled");

    int64_t pins_start = hyperedge(he).firstEntry();
    int64_t pins_end = pins_start + size;
    int64_t slot_of_u = pins_end;

    for ( int64_t pos = pins_end - 1; pos >= pins_start; --pos ) {
      if ( _incidence_array[pos] == u ) {
        slot_of_u = pos;
        break;
      }
    }

    ASSERT(slot_of_u != pins_end, "Hypernode" << u << "not found in hyperedge" << he);
    ASSERT(_incidence_array[slot_of_u] == u, "Hypernode" << u << "not found in hyperedge" << he);
    _incidence_array[slot_of_u] = v;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID he,
                                                                       const HypernodeID hn) {
    using std::swap;
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode" << hn << "is disabled");

    auto begin = incident_nets(hn).begin();
    ASSERT(incident_nets(hn).size() > 0);
    auto last_entry = incident_nets(hn).end() - 1;
    while (*begin != he) {
      ++begin;
    }
    ASSERT(begin < incident_nets(hn).end());
    swap(*begin, *last_entry);
    incident_nets(hn).pop_back();
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID he,
                                                                       const HypernodeID hn,
                                                                       const PartitionID community_id,
                                                                       const std::vector<Self>& hypergraphs,
                                                                       const bool invalidate_only = false) {
    unused(community_id);
    unused(hypergraphs);
    using std::swap;
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode" << hn << "is disabled");

    auto& incident_nets_of_hn = incident_nets(hn);
    int64_t invalid_community_nets = hypernode(hn).invalidCommunityNets();
    int64_t pos = invalid_community_nets - 1;
    for ( ; pos >= 0; --pos ) {
      if ( incident_nets_of_hn[pos] == he ) {
        break;
      }
    }
    ASSERT(incident_nets_of_hn[pos] == he, "Hyperedge" << he << "not found");

    // If hyperedge is in single-pin part of incident nets we first swap it
    // to valid part incident nets
    int64_t single_pin_community_nets = hypernode(hn).singlePinCommunityNets();
    ASSERT(pos != -1);
    if ( pos < single_pin_community_nets ) {
      ASSERT(single_pin_community_nets > 0);
      std::swap(incident_nets_of_hn[pos], incident_nets_of_hn[single_pin_community_nets - 1]);
      pos = single_pin_community_nets - 1;
      hypernode(hn).decrementSinglePinCommunityNets();
    }

    // Afterwards, we swap incident nets to invalid part ...
    std::swap(incident_nets_of_hn[pos], incident_nets_of_hn[invalid_community_nets - 1]);
    pos = invalid_community_nets - 1;
    hypernode(hn).decrementInvalidCommunityNets();

    if ( !invalidate_only ) {
      // ... and if hyperedge should be removed from incident nets, we swap it to the end
      // and pop back.
      std::swap(incident_nets_of_hn[pos], incident_nets_of_hn.back());
      incident_nets_of_hn.pop_back();
    } else {
      ASSERT(hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled());
    }

    ASSERT([&] {
      size_t invalid_community_nets = hypernode(hn).invalidCommunityNets();
      for ( size_t i = 0; i < invalid_community_nets; ++i) {
        const HyperedgeID he = incident_nets_of_hn[i];
        if ( hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled() ) {
          LOG << "HE" << he << "should be in invalid part of incident nets of HN" << hn;
          return false;
        }
      }
      return true;
    }(), "There is an invalidated community hyperedge in valid part of incident nets");

    ASSERT([&] {
      size_t invalid_community_nets = hypernode(hn).invalidCommunityNets();
      for ( size_t i = invalid_community_nets; i < incident_nets_of_hn.size(); ++i) {
        const HyperedgeID he = incident_nets_of_hn[i];
        if ( !hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled() ) {
          LOG << "HE" << he << "should be in valid part of incident nets of HN" << hn;
          return false;
        }
      }
      return true;
    }(), "There is an valid community hyperedge in invalid part of incident nets");
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID he,
                                                                     const HypernodeID hn) {
    size_t invalid_incident_nets = hypernode(hn).invalidIncidentNets();
    ASSERT(std::count(incident_nets(hn).begin() + invalid_incident_nets,
                      incident_nets(hn).end(), he)
            == 0,
            "HN" << hn << "is already connected to HE" << he);

    auto& incident_nets_of_hn = incident_nets(hn);
    size_t slot_of_he = invalid_incident_nets;
    for ( size_t pos = 0; pos < invalid_incident_nets; ++pos ) {
      if ( incident_nets_of_hn[pos] == he ) {
        slot_of_he = pos;
        break;
      }
    }

    if ( slot_of_he < invalid_incident_nets ) {
      ASSERT(incident_nets_of_hn[slot_of_he] == he);
      std::swap(incident_nets_of_hn[slot_of_he],
                incident_nets_of_hn[invalid_incident_nets - 1]);
      hypernode(hn).decrementInvalidIncidentNets();
    } else {
      incident_nets_of_hn.push_back(he);
    }

    // TODO(heuer): increment pin count in part
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeDisabledHyperedgesFromIncidentNets(const HypernodeID hn,
                                                                                const std::vector<Self>& hypergraphs) {
    size_t invalid_incident_nets = 0;
    size_t incident_nets_end = incident_nets(hn).size();
    for ( size_t pos = 0; pos < incident_nets_end; ++pos  ) {
      const HyperedgeID he = incident_nets(hn)[pos];
      if ( !hypergraph_of_hyperedge(he, hypergraphs).edgeIsEnabled(he) ) {
        std::swap(incident_nets(hn)[pos], incident_nets(hn)[invalid_incident_nets++]);
      }
    }
    hypernode(hn).setInvalidIncidentNets(invalid_incident_nets);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE UncontractionCase get_uncontraction_case(const HyperedgeID he,
                                                                           const size_t size,
                                                                           const HypernodeID v) const {
    size_t incidence_array_start = hyperedge(he).firstEntry();
    if ( incidence_array_start + size < hyperedge(he + 1).firstEntry() &&
         _incidence_array[incidence_array_start + size] == v ) {
      return UncontractionCase::CASE_1;
    } else {
      return UncontractionCase::CASE_2;
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID get_global_node_id() {
    HypernodeID local_node_id = _next_node_id++;
    return ( ( (HypernodeID) _node ) << NUMA_NODE_INDENTIFIER ) | local_node_id;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID get_global_node_id(const HypernodeID local_id) const {
    return ( ( (HypernodeID) _node ) << NUMA_NODE_INDENTIFIER ) | local_id;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID get_global_edge_id(const size_t edge_pos) const {
    return ( ( (HyperedgeID) _node ) << NUMA_NODE_INDENTIFIER ) | edge_pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_node_id_of_vertex(const HypernodeID u) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & u;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_edge_id_of_hyperedge(const HyperedgeID e) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & e;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Self& hypergraph_of_vertex(const HypernodeID u,
                                                                          const std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_vertex(u);
    ASSERT(node < (int) hypergraph.size());
    return hypergraph[node];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Self& hypergraph_of_vertex(const HypernodeID u,
                                                                    std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_vertex(u);
    ASSERT(node < (int) hypergraph.size());
    return hypergraph[node];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Self& hypergraph_of_hyperedge(const HyperedgeID e,
                                                                             const std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_hyperedge(e);
    ASSERT(node < (int) hypergraph.size());
    return hypergraph[node];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Self& hypergraph_of_hyperedge(const HyperedgeID e,
                                                                       std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_hyperedge(e);
    ASSERT(node < (int) hypergraph.size());
    return hypergraph[node];
  }

  // ! Accessor for hypernode-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode & hypernode(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[local_id];
  }

  // ! Accessor for hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge & hyperedge(const HyperedgeID e) const {
    // <= instead of < because of sentinel
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[local_id];
  }

  // ! Accessor for community hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const CommunityHyperedge & community_hyperedge(const HyperedgeID e, const PartitionID community_id) const {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hyperedges, "Hyperedge" << e << "does not exist");

    size_t community_hyperedge_position = find_position_of_community_hyperedge(e, community_id);
    ASSERT(community_hyperedge_position < _community_hyperedges[local_id].size(),
           "Community hyperedge" << e << "with community id" << community_id << "not found");
    return _community_hyperedges[local_id][community_hyperedge_position];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t find_position_of_community_hyperedge(const HyperedgeID e, const PartitionID community_id) const {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hyperedges, "Hyperedge" << e << "does not exist");

    size_t pos = 0;
    for ( ; pos < _community_hyperedge_ids[local_id].size(); ++pos ) {
      if ( _community_hyperedge_ids[local_id][pos] == community_id ) {
        break;
      }
    }

    ASSERT(pos < _community_hyperedges[local_id].size(),
           "Community hyperedge" << e << "with community id" << community_id << "not found");
    return pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_nets[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode & hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const StreamingHypergraph&>(*this).hypernode(u));
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CommunityHyperedge & community_hyperedge(const HyperedgeID e, const PartitionID community_id) {
    return const_cast<CommunityHyperedge&>(static_cast<const StreamingHypergraph&>(*this).community_hyperedge(e, community_id));
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge & hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const StreamingHypergraph&>(*this).hyperedge(e));
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) {
    return const_cast<parallel::scalable_vector<HyperedgeID>&>(static_cast<const StreamingHypergraph&>(*this).incident_nets(u));
  }

  const int _node;
  // ! task arena for numa node
  tbb::task_arena& _arena;
  // ! If true, than single-pin community hyperedges are removed
  bool _remove_single_pin_community_nets;

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;

  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! Hypernodes
  parallel::scalable_vector<Hypernode> _hypernodes;
  // ! Incident nets of hypernodes
  IncidentNets _incident_nets;
  // ! Hyperedges
  parallel::scalable_vector<Hyperedge> _hyperedges;
  // ! Pins of hyperedges
  parallel::scalable_vector<HypernodeID> _incidence_array;
  // ! Community Ids of community hyperedges
  std::vector<parallel::scalable_vector<PartitionID>> _community_hyperedge_ids;
  // ! Community Hyperedges
  CommunityHyperedges _community_hyperedges;

  // ! Will be used during uncontraction to mark
  // ! all incident nets of contraction partner v
  kahypar::ds::FastResetFlagArray<> _incident_nets_of_v;

  // ! Contains for each hypernode, how many time it
  // ! occurs as pin in incidence array
  parallel::scalable_vector<size_t> _vertex_pin_count;

  // ! Stream for pins of hyperedges
  StreamingVector<HypernodeID> _pin_stream;
  // ! Stream for hyperedges
  StreamingVector<Hyperedge> _hyperedge_stream;

  // ! Atomic counter for unique local node ids
  std::atomic<HypernodeID> _next_node_id;
  // ! Stream for hypernodes
  StreamingVector<Hypernode> _hypernode_stream;

  // ! Stream for incident nets
  StreamingMap<HypernodeID, HyperedgeID> _incident_net_stream;
};

} // namespace ds
} // namespace mt_kahypar