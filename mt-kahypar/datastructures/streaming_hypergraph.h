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
#include <atomic>
#include <chrono>
#include <functional>
#include <set>
#include <thread>
#include <type_traits>

#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"
#include "tbb/queuing_mutex.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {
// forward
template <typename HardwareTopology,
          typename TBBNumaArena>
class Hypergraph;

template <typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class StreamingHypergraph {
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr size_t NUMA_NODE_INDENTIFIER = 48;
  // seed for edge hashes used for parallel net detection
  static constexpr size_t kEdgeHashSeed = 42;

  using HypernodeAtomic = parallel::IntegralAtomicWrapper<HypernodeID>;
  using HyperedgeAtomic = parallel::IntegralAtomicWrapper<HyperedgeID>;
  using PartitionAtomic = parallel::IntegralAtomicWrapper<PartitionID>;

  static constexpr PartitionID kInvalidPartition = -1;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

  using Self = StreamingHypergraph<HardwareTopology, TBBNumaArena>;

  using IncidentNets = parallel::scalable_vector<parallel::scalable_vector<HyperedgeID> >;
  using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;

  using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

  static_assert(sizeof(HypernodeID) == 8, "Hypernode ID must be 8 byte");
  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(sizeof(HyperedgeID) == 8, "Hyperedge ID must be 8 byte");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  struct AtomicHypernodeData {
    AtomicHypernodeData() :
      part_id(kInvalidPartition),
      num_incident_cut_hes(0) { }

    PartitionAtomic part_id;
    HyperedgeAtomic num_incident_cut_hes;
  };

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex (except incident nets and block information).
   */
  class Hypernode {
   public:
    using IDType = HyperedgeID;

    Hypernode() :
      _id(kInvalidHypernode),
      _original_id(kInvalidHypernode),
      _community_node_id(kInvalidHypernode),
      _weight(1),
      _community_id(0),
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
      _community_id(0),
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
      return _invalid_community_nets;
    }

    void setInvalidCommunityNets(const size_t invalid_community_nets) {
      _invalid_community_nets = invalid_community_nets;
    }

    void incrementInvalidCommunityNets() {
      ++_invalid_community_nets;
    }

    void decrementInvalidCommunityNets() {
      --_invalid_community_nets;
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

  /**
   * Represents a community hyperedge of the hypergraph and contains all information
   * associated with it.
   */
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

    size_t& hash() {
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

  // ! Iterator to iterator over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterator over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the set of incident nets of a hypernode
  // ! or the set of pins of a hyperedge
  using IncidenceIterator = typename parallel::scalable_vector<HypernodeID>::const_iterator;
  // ! Iterator over the set of community ids of a hyperedge
  using CommunityIterator = typename parallel::scalable_vector<PartitionID>::const_iterator;
  // ! Community Hyperedges
  using CommunityHyperedges = parallel::scalable_vector<parallel::scalable_vector<CommunityHyperedge> >;

 public:
  enum class UncontractionCase : uint8_t {
    CASE_1 = 0,
    CASE_1_FORWARD = 1,
    CASE_1_BACKWARD = 2,
    CASE_2 = 3
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

    Memento(HypernodeID representative, HypernodeID contraction_partner) :
      u(representative),
      v(contraction_partner),
      community_id(kInvalidPartition),
      one_pin_hes_begin(0),
      one_pin_hes_size(0),
      parallel_hes_begin(0),
      parallel_hes_size(0) { }

    Memento(HypernodeID representative, HypernodeID contraction_partner, PartitionID community) :
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

  explicit StreamingHypergraph(const int node,
                               const PartitionID k,
                               tbb::task_arena& arena,
                               const bool remove_single_pin_community_nets = true) :
    _node(node),
    _k(k),
    _arena(arena),
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
    _atomic_hn_data(),
    _pins_in_part(),
    _connectivity_sets(0,0),
    _vertex_pin_count(),
    _pin_stream(),
    _hyperedge_stream(),
    _next_node_id(0),
    _hypernode_stream(),
    _incident_net_stream() {
    // Make sure constructor is called on corresponding numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Only allowed to allocate numa hypergraph on node" << _node << ", but it is"
                                                              << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));
  }

  StreamingHypergraph(const StreamingHypergraph&) = delete;
  StreamingHypergraph & operator= (const StreamingHypergraph &) = delete;

  StreamingHypergraph(StreamingHypergraph&& other) :
    _node(other._node),
    _k(other._k),
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
    _atomic_hn_data(std::move(other._atomic_hn_data)),
    _pins_in_part(std::move(other._pins_in_part)),
    _connectivity_sets(std::move(other._connectivity_sets)),
    _vertex_pin_count(std::move(other._vertex_pin_count)),
    _pin_stream(std::move(other._pin_stream)),
    _hyperedge_stream(std::move(other._hyperedge_stream)),
    _next_node_id(other._next_node_id.load()),
    _hypernode_stream(std::move(other._hypernode_stream)),
    _incident_net_stream(std::move(other._incident_net_stream)) { }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Recomputes the total weight of the hypergraph (in parallel)
  void updateTotalWeight() {
    tbb::task_group group;
    _arena.execute([&] {
          group.run([&] {
            _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 0,
                                                 [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
              HypernodeWeight weight = init;
              for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
                weight += this->_hypernodes[hn].weight();
              }
              return weight;
            }, std::plus<HypernodeWeight>());
          });
        });
    group.wait();
  }

    // ! Recomputes the total weight of the hypergraph (in sequential)
  void updateTotalWeightSequential() {
    _total_weight = 0;
    for ( const HypernodeID& hn : nodes() ) {
      _total_weight += nodeWeight(hn);
    }
  }

  // ####################### Iterators #######################

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    HypernodeID start = get_global_node_id(0);
    HypernodeID end = get_global_node_id(_num_hypernodes);
    return IteratorRange<HypernodeIterator>(HypernodeIterator(_hypernodes.data(), start, end),
                          HypernodeIterator((_hypernodes.data() + _num_hypernodes), end, end));
  }

  // ! Returns a range of the active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    HyperedgeID start = get_global_edge_id(0UL);
    HyperedgeID end = get_global_edge_id(_num_hyperedges);
    return IteratorRange<HyperedgeIterator>(HyperedgeIterator(_hyperedges.data(), start, end),
                          HyperedgeIterator((_hyperedges.data() + _num_hyperedges), end, end));
  }

  /**
   * For a more detailed overview of different incidentEdges(..) function please have a look
   * at the documentation in hypergraph.
   */

  // ! During parallel community coarsening we do not remove parallel (identical) hyperedges that span more than one community.
  // ! Instead, we just invalidate those hyperedges, which means swapping them to the end of incident nets and storing a pointer to all invalidated hyperedges.

  // ! Returns a range to loop over the set of all VALID and INVALID incident hyperedges of hypernode u.
  IteratorRange<IncidenceIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return IteratorRange<IncidenceIterator>(incident_nets(u).cbegin() + hypernode(u).invalidIncidentNets(), incident_nets(u).cend());
  }

  // ! Returns a range to loop over the set of all VALID incident hyperedges of hypernode u.
  IteratorRange<IncidenceIterator> validIncidentEdges(const HypernodeID u, const PartitionID) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    ASSERT(hypernode(u).invalidCommunityNets() <= incident_nets(u).size());
    return IteratorRange<IncidenceIterator>(incident_nets(u).cbegin(), incident_nets(u).cbegin() + hypernode(u).invalidCommunityNets());
  }


  // TODO function name should reflect its purpose
  // ! Returns a range to loop over the set of all VALID incident hyperedges of hypernode u that are not single-pin community hyperedges.
  IteratorRange<IncidenceIterator> incidentEdges(const HypernodeID u, const PartitionID) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    ASSERT(hypernode(u).singlePinCommunityNets() <= hypernode(u).invalidCommunityNets());
    ASSERT(hypernode(u).invalidCommunityNets() <= incident_nets(u).size());
    return IteratorRange<IncidenceIterator>(incident_nets(u).cbegin() + hypernode(u).singlePinCommunityNets(), incident_nets(u).cbegin() + hypernode(u).invalidCommunityNets());
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  // ! Note, this function fails if community hyperedges are initialized.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(!hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are initialized");
    return IteratorRange<IncidenceIterator>(_incidence_array.cbegin() + hyperedge(e).firstEntry(), _incidence_array.cbegin() + hyperedge(e).firstInvalidEntry());
  }

  // ! Returns a range to loop over the pins of hyperedge e that belong to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    const CommunityHyperedge& community_he = community_hyperedge(e, community_id);
    return IteratorRange<IncidenceIterator>(_incidence_array.cbegin() + community_he.firstEntry(), _incidence_array.cbegin() + community_he.firstInvalidEntry());
  }

  // ! Returns a range to loop over the set of communities contained in hyperedge e.
  IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(local_id < _community_hyperedge_ids.size());
    return IteratorRange<CommunityIterator>(_community_hyperedge_ids[local_id].cbegin(), _community_hyperedge_ids[local_id].cend());
  }

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivitySets::Iterator> connectivitySet(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    return _connectivity_sets.connectivitySet(local_id);
  }

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! defined in the input file
  HypernodeID originalNodeId(const HypernodeID u) const {
    return hypernode(u).originalNodeId();
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return incident_nets(u).size() - hypernode(u).invalidIncidentNets();
  }

  // ! Number of invalid incident nets
  HyperedgeID numInvalidIncidentNets(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).invalidIncidentNets();
  }

  // ! Number of hyperedges a vertex occurs as a pin in this hypergraph
  size_t vertexPinCount(const HypernodeID hn) const {
    ASSERT(hn < _vertex_pin_count.size());
    return _vertex_pin_count[hn];
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  // ####################### Hyperedge Information #######################

  // ! Returns for a hyperedge of the hypergraph its original hyperedge id
  // ! defined in the input file
  HypernodeID originalEdgeId(const HyperedgeID e) const {
    return hyperedge(e).originalEdgeId();
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).weight();
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    hyperedge(e).setWeight(weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if (hyperedge(e).isInitCommunityHyperedges()) {
      HypernodeID size = 0;
      HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
      ASSERT(local_id < _community_hyperedges.size());
      for (const CommunityHyperedge& community_he : _community_hyperedges[local_id]) {
        size += community_he.size();
      }
      return size;
    } else {
      return hyperedge(e).size();
    }
  }

  // ! Hash value defined over the pins of a hyperedge
  size_t edgeHash(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if (hyperedge(e).isInitCommunityHyperedges()) {
      size_t hash = 0;
      HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
      ASSERT(local_id < _community_hyperedges.size());
      for (const CommunityHyperedge& community_he : _community_hyperedges[local_id]) {
        hash += community_he.hash();
      }
      return hash;
    } else {
      return hyperedge(e).hash();
    }
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return !hyperedge(e).isDisabled();
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hyperedge(e).enable();
  }

  // ! Disables a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
  }

  // ####################### Community Hyperedge Information #######################

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).weight();
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).setWeight(weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).size();
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    return community_hyperedge(e, community_id).hash();
  }

  // ! Disables a community hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e, const PartitionID community_id) {
    community_hyperedge(e, community_id).disable();
  }

  // ####################### Community Information #######################

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityID();
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID hn, const PartitionID community_id) {
    hypernode(hn).setCommunityID(community_id);
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  HypernodeID communityNodeId(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityNodeId();
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(hyperedge(e).isInitCommunityHyperedges(), "Community hyperedges of HE" << e << "are not initialized");
    HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    return _community_hyperedges[local_id].size();
  }

  // ####################### Partition Information #######################

  // ! Sets the block id of an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u, PartitionID id) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    PartitionID invalid_id = kInvalidPartition;
    return _atomic_hn_data[get_local_node_id_of_vertex(u)].part_id.compare_and_exchange_strong(invalid_id, id);
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u, PartitionID from, PartitionID to) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    // TODO why do we need CAS for part IDs if each vertex can only be moved by one processor?
    return _atomic_hn_data[get_local_node_id_of_vertex(u)].part_id.compare_and_exchange_strong(from, to);
  }

  // ! Block which vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _atomic_hn_data[get_local_node_id_of_vertex(u)].part_id;
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _atomic_hn_data[get_local_node_id_of_vertex(u)].num_incident_cut_hes > 0;
  }

  // ! Number of incident cut hyperedges of vertex u
  HyperedgeID numIncidentCutHyperedges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _atomic_hn_data[get_local_node_id_of_vertex(u)].num_incident_cut_hes;
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  void initializeNumCutHyperedges(const std::vector<Self>& hypergraphs) {
    _arena.execute([&] {
          tbb::parallel_for(0UL, this->_num_hypernodes, [&](const HypernodeID& id) {
            const HypernodeID hn = get_global_node_id(id);
            if ( nodeIsEnabled(hn) ) {
              ASSERT(_atomic_hn_data[get_local_node_id_of_vertex(hn)].num_incident_cut_hes == 0);
              for (const HyperedgeID& he : incidentEdges(hn)) {
                if (hypergraph_of_hyperedge(he, hypergraphs).connectivity(he) > 1) {
                  incrementIncidentNumCutHyperedges(hn);
                }
              }
            }
          });
        });
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    return _connectivity_sets.connectivity(local_id);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID id) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(id < _k && id != kInvalidPartition, "Partition ID" << id << "is out of bounds");
    const HyperedgeID local_id = get_local_edge_id_of_hyperedge(e);
    return _pins_in_part[static_cast<size_t>(local_id) * _k + id];
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    // Reset partition ids
    for ( AtomicHypernodeData& hn_data : _atomic_hn_data ) {
      hn_data.part_id = kInvalidPartition;
      hn_data.num_incident_cut_hes = 0;
    }

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      const HyperedgeID local_id = get_local_edge_id_of_hyperedge(he);
      for ( const PartitionID& part_id : connectivitySet(he) ) {
        _pins_in_part[static_cast<size_t>(local_id) * _k + part_id] = 0;
      }
      _connectivity_sets.clear(local_id);
    }
  }

  // ####################### Contract / Uncontract #######################

  /*!
   * Contracts the vertex pair (u,v) inside hyperedge e. The representative u remains
   * in the hypergraph. The contraction partner v is removed from the hypergraph.
   *
   * For each hyperedge e incident to v, a contraction lead to one of two operations:
   * 1.) If e contained both u and v, then v is removed from e.
   * 2.) If e only contained v, than the slot of v in the incidence structure of e
   *     is reused to store u.
   *
   * NOTE, this function is not thread-safe and should be only called in a single-threaded
   * setting.
   *
   * \param u Representative hypernode that will remain in the hypergraph
   * \param v Contraction partner that will be removed from the hypergraph
   * \param e hyperedge which contains v
   * \param hypergraph_of_u corresponding StreamingHypergraph that contains u
   */
  void contract(const HypernodeID u, const HypernodeID v,
                const HyperedgeID e, Self& hypergraph_of_u) {
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node);
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

    if (slot_of_u != last_pin_slot) {
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

  /*!
   * Contracts the vertex pair (u,v) belonging to the same community inside hyperedge e.
   * The representative u remains in the hypergraph. The contraction partner
   * v is removed from the hypergraph.
   *
   * For each hyperedge e incident to v, a contraction lead to one of two operations:
   * 1.) If e contained both u and v, then v is removed from e.
   * 2.) If e only contained v, than the slot of v in the incidence structure of e
   *     is reused to store u.
   *
   * NOTE, in order that this function works correct, community hyperedges have to be
   * initialized beforehand. This function is thread-safe as long as only one thread
   * performs contractions in one community.
   *
   * \param u Representative hypernode that will remain in the hypergraph
   * \param v Contraction partner that will be removed from the hypergraph
   * \param e hyperedge which contains v
   * \param community_id Community to which u and v belongs to
   * \param hypergraph_of_u corresponding StreamingHypergraph that contains u
   */
  void contract(const HypernodeID u, const HypernodeID v,
                const HyperedgeID e, const PartitionID community_id,
                Self& hypergraph_of_u) {
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node);

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

    if (slot_of_u != last_pin_slot) {
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

  /*!
  * Undoes a contraction operation (u,v) inside hyperedge e.
  *
  * NOTE, this function is only thread-safe, if we are in batch uncontraction mode
  * and the batch fullfils the requirement decribed in the batch uncontractions
  * function of hypergraph.h
  */
  bool uncontract(const HypernodeID u, const HypernodeID v,
                  const HyperedgeID e, const size_t incident_nets_pos,
                  std::vector<Self>& hypergraphs,
                  const kahypar::ds::FastResetFlagArray<>* batch_hypernodes = nullptr) {
    using std::swap;

    bool is_batch_uncontraction = batch_hypernodes != nullptr;
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if (containsIncidentNet(e)) {
      // ... then we have to do some kind of restore operation.
      UncontractionCase uncontraction_case = is_batch_uncontraction ?
                                             get_uncontraction_case(e, hyperedge(e).size(), v, hypergraphs, *batch_hypernodes) :
                                             get_uncontraction_case(e, hyperedge(e).size(), v);
      if (uncontraction_case != UncontractionCase::CASE_2) {
        // hyperedge(he + 1) always exists because of sentinel
        // Undo case 1 operation (i.e. Pin v was just cut off by decreasing size of HE e)
        DBG << V(e) << " -> case 1";

        if (uncontraction_case == UncontractionCase::CASE_1) {
          hyperedge(e).incrementSize();
        } else if (uncontraction_case == UncontractionCase::CASE_1_FORWARD) {
          ASSERT(is_batch_uncontraction);
          incrementHyperedgeSizeDuringBatchUncontraction(e, hypergraphs, *batch_hypernodes);
        }

        incrementPinCountInPart(e, hypergraph_of_vertex(v, hypergraphs).partID(v));

        if (connectivity(e) > 1) {
          hypergraph_of_vertex(v, hypergraphs).incrementIncidentNumCutHyperedges(v);
        }

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

        if (connectivity(e) > 1) {
          hypergraph_of_vertex(u, hypergraphs).decrementIncidentNumCutHyperedges(u);
          hypergraph_of_vertex(v, hypergraphs).incrementIncidentNumCutHyperedges(v);
        }

        return true;
      }
    }
    return false;
  }

  #ifndef NDEBUG
  /*!
  * Undoes a contraction operation (u,v) inside hyperedge e.
  *
  * This function uncontracts a disabled parallel hyperedge. Hyperedge e and
  * representative are parallel and e is disabled (whereas representative is enabled).
  * This function is only executed in debug to verify, that e is parallel to its representative
  * before restoring e. In Release, we just copy the active pins (via memcpy) to hyperedge e.
  * For more detailed information how we detect and process parallel hyperedges during coarsening
  * please have a look at hypergraph.h, community_coarsener_base.h and hypergraph_pruner.h
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  bool uncontract(const HypernodeID u, const HypernodeID v,
                  const HyperedgeID e, const HyperedgeID representative,
                  const size_t incident_nets_pos,
                  std::vector<Self>& hypergraphs,
                  const kahypar::ds::FastResetFlagArray<>* batch_hypernodes = nullptr) {
    using std::swap;
    ASSERT(hyperedge(e).isDisabled(), "Hyperedge" << e << "is enabled");
    const Self& hypergraph_of_rep = hypergraph_of_hyperedge(representative, hypergraphs);
    ASSERT(hypergraph_of_rep.edgeIsEnabled(representative), "Hyperedge" << representative << "is disabled");

    bool is_batch_uncontraction = batch_hypernodes != nullptr;
    if (containsIncidentNet(e)) {
      ASSERT(hypergraph_of_rep.edgeSize(representative) <=
             hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry(),
             V(hypergraph_of_rep.edgeSize(representative))
             << V((hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry())));

      size_t edge_size = hypergraph_of_rep.edgeSize(representative);
      UncontractionCase uncontraction_case = is_batch_uncontraction ?
                                             get_uncontraction_case(e, edge_size, v, hypergraphs, *batch_hypernodes) :
                                             get_uncontraction_case(e, edge_size, v);
      if (uncontraction_case == UncontractionCase::CASE_2) {
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

        return true;
      }
    }

    return false;
  }
  #endif

  // ####################### Streaming Functions #######################

  // ! Adds a hypernode to the hypergraph
  // ! Note, after all hypernodes are streamed into the streaming hypergraphs
  // ! initializeHypernodes() have to be called.
  HypernodeID streamHypernode(const HypernodeID original_id,
                              const HypernodeWeight weight) {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node);
    HypernodeID node_id = get_global_node_id();
    _hypernode_stream.stream(node_id, original_id, weight);
    return node_id;
  }

  // ! Adds a hyperedge to the hypergraph
  // ! Note, after all hyperedges are streamed into the streaming hypergraphs
  // ! initializeHyperedges() have to be called.
  void streamHyperedge(const parallel::scalable_vector<HypernodeID>& hyperedge,
                       const HyperedgeID original_id,
                       const HyperedgeWeight weight) {
    int cpu_id = sched_getcpu();
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(cpu_id) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was CPU" << cpu_id
                                                        << "is on node" << HardwareTopology::instance().numa_node_of_cpu(cpu_id));
    _hyperedge_stream.stream(_pin_stream.size(cpu_id), hyperedge.size(), original_id, weight);
    for (const HypernodeID& pin : hyperedge) {
      _pin_stream.stream(pin);
    }
  }

  // ! Adds a hyperedge to the hypergraph
  // ! Note, adding a hyperedge, the function initializeHyperedges(num_he, num_pins) has
  // ! to be called.
  void addHyperedge(const parallel::scalable_vector<HypernodeID>& hyperedge,
                    const HyperedgeID original_id,
                    const HyperedgeWeight weight,
                    const HyperedgeID he_idx,
                    const HypernodeID pin_idx) {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was CPU" << sched_getcpu()
                                                        << "is on node" << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));
    ASSERT(he_idx < _hyperedges.size(), V(he_idx) << V(_hyperedges.size()));
    ASSERT(pin_idx + hyperedge.size() <= _incidence_array.size());
    _hyperedges[he_idx] = Hyperedge(pin_idx, hyperedge.size(), original_id, weight);
    memcpy(_incidence_array.data() + pin_idx, hyperedge.data(), sizeof(HypernodeID) * hyperedge.size());
  }

  // ! Indicates that hypernode hn is incident to net he.
  // ! Note, after all incident nets are streamed into the streaming hypergraphs
  // ! initializeIncidentNets() have to be called.
  void streamIncidentNet(const HypernodeID hn,
                         const HyperedgeID he) {
    ASSERT(_node == -1 || get_numa_node_of_vertex(hn) == _node);
    _incident_net_stream.stream(hn, he);
  }

  // ####################### Initialization / Reset Functions #######################

  /*!
  * This function initializes all hypernodes that are streamed into this hypergraph.
  * This includes:
  *   1.) Replace vertex ids in incidence array with new numa vertex ids
  *   2.) Stream all incident nets part of this hypergraph to their corresponding
  *       streaming hypergraphs.
  */
  void initializeHypernodes(std::vector<Self>& hypergraphs,
                            const std::vector<HypernodeID>& node_mapping) {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was on node"
                                                        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _hypernodes = _hypernode_stream.copy(_arena);
    _num_hypernodes = _hypernodes.size();
    _atomic_hn_data.resize(_num_hypernodes);

    // Sort hypernodes in increasing order of their node id and ...
    utils::Timer::instance().start_timer("sort_and_remap_node_ids", "Sort and Remap Nodes", true);
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
            tbb::parallel_for(0UL, _incidence_array.size(), [&](const size_t& pos) {
              HypernodeID pin = _incidence_array[pos];
              ASSERT(pin < node_mapping.size());
              _incidence_array[pos] = node_mapping[pin];
            });
          });
        });
    group.wait();
    utils::Timer::instance().stop_timer("sort_and_remap_node_ids");

    // Compute Total Hypergraph Weight
    utils::Timer::instance().start_timer("compute_total_weight", "Compute Total Weight", true);
    updateTotalWeight();
    utils::Timer::instance().stop_timer("compute_total_weight");

    ASSERT([&] {
          for (size_t i = 0; i < _hypernodes.size(); ++i) {
            const Hypernode& hn = _hypernodes[i];
            const int node = get_numa_node_of_vertex(hn.nodeId());
            const HypernodeID local_id = get_local_node_id_of_vertex(hn.nodeId());
            if (_node != -1 && node != _node) {
              LOG << "Hypernode" << hn.nodeId() << "should be on numa node" << node
                  << ", but is on numa node" << _node;
              return false;
            } else if (i != local_id) {
              LOG << "Hypernode" << local_id << "should have local id" << i;
              return false;
            }

            if (i > 0) {
              const Hypernode& previous_hn = _hypernodes[i - 1];
              const HypernodeID previous_local_id = get_local_node_id_of_vertex(previous_hn.nodeId());
              if (local_id <= previous_local_id) {
                LOG << "Local node ids not sorted" << V(previous_local_id) << V(local_id);
                return false;
              } else if (previous_local_id + 1 != local_id) {
                LOG << "Local node ids not consecutive" << V(previous_local_id) << V(local_id);
                return false;
              }
            }
          }
          return true;
        } (), "Initialization of hypernodes failed");

    // Stream incident nets
    utils::Timer::instance().start_timer("stream_incident_nets", "Stream Incident Nets", true);
    _arena.execute([&] {
          group.run([&] {
            tbb::parallel_for(0UL, _hyperedges.size(), [&](const size_t& pos) {
              Hyperedge& he = _hyperedges[pos];
              const HyperedgeID he_id = get_global_edge_id(pos);
              for (size_t incidence_array_pos = he.firstEntry();
                   incidence_array_pos < he.firstEntry() + he.size();
                   ++incidence_array_pos) {
                const HypernodeID pin = _incidence_array[incidence_array_pos];
                hypergraph_of_vertex(pin, hypergraphs).streamIncidentNet(pin, he_id);
                // Initialize edge hash
                he.hash() += kahypar::math::hash(pin);
              }
            });
          });
        });
    group.wait();
    utils::Timer::instance().stop_timer("stream_incident_nets");

    _hypernode_stream.clear();
  }

  /*!
  * This function initializes all hypernodes that are streamed into this hypergraph.
  * This includes:
  *   1.) Replace vertex ids in incidence array with new numa vertex ids
  *   2.) Stream all incident nets part of this hypergraph to their corresponding
  *       streaming hypergraphs.
  */
  void initializeHypernodesSequential(const std::vector<HypernodeID>& node_mapping) {
    ASSERT(_node == -1);

    _hypernodes = _hypernode_stream.copy();
    _num_hypernodes = _hypernodes.size();
    _atomic_hn_data.resize(_num_hypernodes);

    // Sort hypernodes in increasing order of their node id and ...
    utils::Timer::instance().start_timer("sort_and_remap_node_ids", "Sort and Remap Nodes", true);
    std::sort(_hypernodes.begin(), _hypernodes.end(),
              [&](const Hypernode& lhs, const Hypernode& rhs) {
                return lhs.nodeId() < rhs.nodeId();
              });

    // ... remap node ids in incidence array to new node ids
    for ( size_t pos = 0; pos < _incidence_array.size(); ++pos ) {
      HypernodeID pin = _incidence_array[pos];
      ASSERT(pin < node_mapping.size());
      _incidence_array[pos] = node_mapping[pin];
    }
    utils::Timer::instance().stop_timer("sort_and_remap_node_ids");

    // Compute Total Hypergraph Weight
    utils::Timer::instance().start_timer("compute_total_weight", "Compute Total Weight", true);
    updateTotalWeightSequential();
    utils::Timer::instance().stop_timer("compute_total_weight");

    ASSERT([&] {
          for (size_t i = 0; i < _hypernodes.size(); ++i) {
            const Hypernode& hn = _hypernodes[i];
            const int node = get_numa_node_of_vertex(hn.nodeId());
            const HypernodeID local_id = get_local_node_id_of_vertex(hn.nodeId());
            if (_node != -1 && node != _node) {
              LOG << "Hypernode" << hn.nodeId() << "should be on numa node" << node
                  << ", but is on numa node" << _node;
              return false;
            } else if (i != local_id) {
              LOG << "Hypernode" << local_id << "should have local id" << i;
              return false;
            }

            if (i > 0) {
              const Hypernode& previous_hn = _hypernodes[i - 1];
              const HypernodeID previous_local_id = get_local_node_id_of_vertex(previous_hn.nodeId());
              if (local_id <= previous_local_id) {
                LOG << "Local node ids not sorted" << V(previous_local_id) << V(local_id);
                return false;
              } else if (previous_local_id + 1 != local_id) {
                LOG << "Local node ids not consecutive" << V(previous_local_id) << V(local_id);
                return false;
              }
            }
          }
          return true;
        } (), "Initialization of hypernodes failed");

    // Stream incident nets
    utils::Timer::instance().start_timer("stream_incident_nets", "Stream Incident Nets", true);
    _incident_nets.resize(_hypernodes.size());
    for ( size_t pos = 0; pos < _hyperedges.size(); ++pos ) {
      Hyperedge& he = _hyperedges[pos];
      const HyperedgeID he_id = get_global_edge_id(pos);
      for (size_t incidence_array_pos = he.firstEntry();
            incidence_array_pos < he.firstEntry() + he.size();
            ++incidence_array_pos) {
        const HypernodeID pin = _incidence_array[incidence_array_pos];
        const HypernodeID local_pin_id = get_local_node_id_of_vertex(pin);
        ASSERT(local_pin_id < _incident_nets.size());
        _incident_nets[local_pin_id].push_back(he_id);
        // Initialize edge hash
        he.hash() += kahypar::math::hash(pin);
      }
    }
    utils::Timer::instance().stop_timer("stream_incident_nets");

    _hypernode_stream.clear();
  }

  /*!
  * This function initializes all hyperedges that are streamed into this hypergraph in parallel.
  * This includes:
  *   1.) Setting up incidence array of this hypergraph
  *   2.) Compute vertex pin counts of all pins contained in this hypergraph
  */
  void initializeHyperedges(const HypernodeID num_hypernodes,
                            const bool compute_vertex_pin_counts = true) {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was on node"
                                                        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    // Copy streamed data into global vectors
    utils::Timer::instance().start_timer("copy_incidence_array_and_he", "Copy Incidence Array and HEs", true);
    _incidence_array = _pin_stream.copy(_arena);
    _hyperedges = _hyperedge_stream.copy(_arena);
    _num_pins = _incidence_array.size();
    _num_hyperedges = _hyperedges.size();
    utils::Timer::instance().stop_timer("copy_incidence_array_and_he");

    ASSERT(_k > 0);
    tbb::parallel_invoke([&] {
      _pins_in_part.assign(_num_hyperedges * _k, HypernodeAtomic(0));
    }, [&] {
      _connectivity_sets = ConnectivitySets(_num_hyperedges, _k);
    }, [&] {
      ThreadLocalFastResetFlagArray tmp_incidence_nets_of_v(_num_hyperedges);
      _incident_nets_of_v = std::move(tmp_incidence_nets_of_v);
    });

    // Update start position of each hyperedge to correct one in global incidence array
    // Note, start positions are stored relative to the local buffer there are streamed into.
    // However, memcpy does hyperedges invalidates those local positions.
    utils::Timer::instance().start_timer("update_start_position", "Update Start Positions of HEs", true);
    tbb::task_group group;
    _arena.execute([&] {
          for (size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id) {
            group.run([&, cpu_id] {
              size_t start = _hyperedge_stream.prefix_sum(cpu_id);
              size_t end = start + _hyperedge_stream.size(cpu_id);
              size_t delta = _pin_stream.prefix_sum(cpu_id);
              ASSERT(end <= _hyperedges.size());
              for (size_t pos = start; pos < end; ++pos) {
                _hyperedges[pos].setFirstEntry(_hyperedges[pos].firstEntry() + delta);
              }
            });
          }
        });
    group.wait();
    // Emplace Back Sentinel
    _hyperedges.emplace_back(_incidence_array.size(), 0, 0UL, 0);
    utils::Timer::instance().stop_timer("update_start_position");

    ASSERT([&] {
          for (size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id) {
            size_t start = _hyperedge_stream.prefix_sum(cpu_id);
            for (size_t idx = 0; idx < _hyperedge_stream.size(cpu_id); ++idx) {
              const Hyperedge& stream_he = _hyperedge_stream.value(cpu_id, idx);
              const Hyperedge& he = _hyperedges[idx + start];
              if (stream_he.size() != he.size()) {
                LOG << "Unequal size";
                return false;
              }
              size_t global_pos = he.firstEntry();
              for (size_t pos = stream_he.firstEntry(); pos < stream_he.firstEntry() + stream_he.size(); ++pos) {
                if (_incidence_array[global_pos] != _pin_stream.value(cpu_id, pos)) {
                  LOG << "Pins in stream and global incidence array are not equal"
                      << V(_incidence_array[global_pos]) << V(_pin_stream.value(cpu_id, pos));
                  return false;
                }
                ++global_pos;
              }
            }
          }
          return true;
        } (), "Failed to copy buffer to hypergraph");

    // Compute how many times a hypernode occurs on this node
    // as pin. Will be later important to compute node assignment.
    // TODO(heuer): Think how to parallelize this
    if ( compute_vertex_pin_counts ) {
      utils::Timer::instance().start_timer("compute_vertex_pin_count", "Compute Vertex Pin Counts", true);
      _vertex_pin_count.resize(num_hypernodes);
      for (size_t i = 0; i < _incidence_array.size(); ++i) {
        const HypernodeID& pin = _incidence_array[i];
        ASSERT(pin < _vertex_pin_count.size());
        _vertex_pin_count[pin]++;
      }
      utils::Timer::instance().stop_timer("compute_vertex_pin_count");
    }

    _pin_stream.clear();
    _hyperedge_stream.clear();
  }

  /*!
   * This function allocates the data structures needed for setting up hyperedge related
   * information. Afterwards, one can call addHyperedge(...) to add hyperedges in a more
   * efficient manner than with streamHyperedge(...)
   */
  void initializeHyperedges(const HyperedgeID num_hyperedges,
                            const HypernodeID num_pins) {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was on node"
                                                        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _num_hyperedges = num_hyperedges;
    _num_pins = num_pins;

    ASSERT(_k > 0);
    tbb::parallel_invoke([&] {
      _hyperedges.resize(num_hyperedges);
    }, [&] {
      _incidence_array.resize(num_pins);
    }, [&] {
      _pins_in_part.assign(_num_hyperedges * _k, HypernodeAtomic(0));
    }, [&] {
      _connectivity_sets = ConnectivitySets(_num_hyperedges, _k);
    }, [&] {
      ThreadLocalFastResetFlagArray tmp_incidence_nets_of_v(_num_hyperedges);
      _incident_nets_of_v = std::move(tmp_incidence_nets_of_v);
    });
  }

  /*!
  * This function initializes all hyperedges that are streamed into this hypergraph non parallel.
  * This includes:
  *   1.) Setting up incidence array of this hypergraph
  *   2.) Compute vertex pin counts of all pins contained in this hypergraph
  */
  void initializeHyperedgesSequential(const HypernodeID num_hypernodes,
                                      const bool compute_vertex_pin_counts = true) {
    // Copy streamed data into global vectors
    utils::Timer::instance().start_timer("copy_incidence_array_and_he", "Copy Incidence Array and HEs", true);
    _incidence_array = _pin_stream.copy();
    _hyperedges = _hyperedge_stream.copy();
    _num_pins = _incidence_array.size();
    _num_hyperedges = _hyperedges.size();
    utils::Timer::instance().stop_timer("copy_incidence_array_and_he");

    ASSERT(_k > 0);
    _pins_in_part.assign(_num_hyperedges * _k, HypernodeAtomic(0));

    _connectivity_sets = ConnectivitySets(_num_hyperedges, _k);

    ThreadLocalFastResetFlagArray tmp_incidence_nets_of_v(_num_hyperedges);
    _incident_nets_of_v = std::move(tmp_incidence_nets_of_v);

    // Update start position of each hyperedge to correct one in global incidence array
    // Note, start positions are stored relative to the local buffer there are streamed into.
    // However, memcpy does hyperedges invalidates those local positions.
    utils::Timer::instance().start_timer("update_start_position", "Update Start Positions of HEs", true);
    for (size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id) {
      size_t start = _hyperedge_stream.prefix_sum(cpu_id);
      size_t end = start + _hyperedge_stream.size(cpu_id);
      size_t delta = _pin_stream.prefix_sum(cpu_id);
      ASSERT(end <= _hyperedges.size());
      for (size_t pos = start; pos < end; ++pos) {
        _hyperedges[pos].setFirstEntry(_hyperedges[pos].firstEntry() + delta);
      }
    }
    // Emplace Back Sentinel
    _hyperedges.emplace_back(_incidence_array.size(), 0, 0UL, 0);
    utils::Timer::instance().stop_timer("update_start_position");

    ASSERT([&] {
          for (size_t cpu_id = 0; cpu_id < std::thread::hardware_concurrency(); ++cpu_id) {
            size_t start = _hyperedge_stream.prefix_sum(cpu_id);
            for (size_t idx = 0; idx < _hyperedge_stream.size(cpu_id); ++idx) {
              const Hyperedge& stream_he = _hyperedge_stream.value(cpu_id, idx);
              const Hyperedge& he = _hyperedges[idx + start];
              if (stream_he.size() != he.size()) {
                LOG << "Unequal size";
                return false;
              }
              size_t global_pos = he.firstEntry();
              for (size_t pos = stream_he.firstEntry(); pos < stream_he.firstEntry() + stream_he.size(); ++pos) {
                if (_incidence_array[global_pos] != _pin_stream.value(cpu_id, pos)) {
                  LOG << "Pins in stream and global incidence array are not equal"
                      << V(_incidence_array[global_pos]) << V(_pin_stream.value(cpu_id, pos));
                  return false;
                }
                ++global_pos;
              }
            }
          }
          return true;
        } (), "Failed to copy buffer to hypergraph");

    // Compute how many times a hypernode occurs on this node
    // as pin. Will be later important to compute node assignment.
    // TODO(heuer): Think how to parallelize this
    if ( compute_vertex_pin_counts ) {
      utils::Timer::instance().start_timer("compute_vertex_pin_count", "Compute Vertex Pin Counts", true);
      _vertex_pin_count.resize(num_hypernodes);
      for (size_t i = 0; i < _incidence_array.size(); ++i) {
        const HypernodeID& pin = _incidence_array[i];
        ASSERT(pin < _vertex_pin_count.size());
        _vertex_pin_count[pin]++;
      }
      utils::Timer::instance().stop_timer("compute_vertex_pin_count");
    }

    _pin_stream.clear();
    _hyperedge_stream.clear();
  }

  // ! This function initializes all incident nets that are streamed into this hypergraph.
  void initializeIncidentNets() {
    // Make sure calling thread is part of correct numa node
    ASSERT(_node == -1 || HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node,
           "Expected that assigned cpu is on numa node" << _node << ", but was on node"
                                                        << HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()));

    _incident_nets.resize(_hypernodes.size());
    _incident_net_stream.copy(_arena, _incident_nets, [&](const HypernodeID& u) {
          ASSERT(_node == -1 || get_numa_node_of_vertex(u) == _node);
          return get_local_node_id_of_vertex(u);
        });
    _incident_net_stream.clear();
  }

  /*!
  * Initializes community hyperedges (must be called before community coarsening).
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to the sorted incidence array
  */
  void initializeCommunityHyperedges(const std::vector<Self>& hypergraphs) {
    auto add_community_hyperedge =
      [&](const HyperedgeID he,
          const PartitionID community_id,
          const size_t start,
          const size_t end,
          const HyperedgeWeight weight) {
        if (community_id != kInvalidPartition) {
          ASSERT(he < this->_community_hyperedges.size());
          ASSERT(start < end);
          this->_community_hyperedge_ids[he].emplace_back(community_id);
          this->_community_hyperedges[he].emplace_back(start, end - start, weight);

          // Compute community hyperedge hash
          for (size_t pos = start; pos < end; ++pos) {
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
            tbb::parallel_for(0UL, this->_num_hyperedges, [&](const HyperedgeID& he) {
              Hyperedge& e = this->_hyperedges[he];
              if (!e.isDisabled()) {
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
                for (size_t incidence_array_pos = incidence_array_start;
                     incidence_array_pos < incidence_array_end;
                     ++incidence_array_pos) {
                  const HypernodeID pin = this->_incidence_array[incidence_array_pos];
                  const PartitionID community_id = hypergraph_of_vertex(pin, hypergraphs).communityID(pin);
                  if (community_id != last_community_id) {
                    add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_pos, e.weight());
                    last_community_start = incidence_array_pos;
                    last_community_id = community_id;
                  }
                }
                add_community_hyperedge(he, last_community_id, last_community_start, incidence_array_end, e.weight());
                e.initializeCommunityHyperedges();

                HEAVY_COARSENING_ASSERT([&] {
                  if (e.firstEntry() != _community_hyperedges[he][0].firstEntry()) {
                    return false;
                  }
                  for (size_t i = 1; i < _community_hyperedges[he].size(); ++i) {
                    if (_community_hyperedges[he][i - 1].firstInvalidEntry() !=
                        _community_hyperedges[he][i].firstEntry()) {
                      return false;
                    }
                  }
                  if (e.firstInvalidEntry() != _community_hyperedges[he].back().firstInvalidEntry()) {
                    return false;
                  }
                  return true;
                } (), "Initialization of community hyperedges failed!");
              }
            });
          });
        });
    group.wait();
  }

  /*!
  * Setting up the incident net vector of each hypernode for community coarsening
  * (must be called before community coarsening and after initializeCommunityHyperedges()).
  * This includes:
  *   1.) Partition the set of incident nets of each vertex in two parts
  *     a.) All single-pin community hyperedges
  *     b.) All other community hyperedges
  * One can iterate over the different part of incident nets via the corresponding
  * incidentEdges(...) function.
  */
  void initializeCommunityHypernodes(const std::vector<Self>& hypergraphs) {
    // Partition incident nets of each hypernode such that in first part are all hyperedges
    // which contains at least to pins and in second part are all single-pin community hyperedges
    tbb::task_group group;
    _arena.execute([&] {
          group.run([&] {
            tbb::parallel_for(0UL, this->_num_hypernodes, [&](const HypernodeID& v) {
              HypernodeID hn = get_global_node_id(v);
              if (!this->hypernode(hn).isDisabled()) {
                PartitionID community_id = this->communityID(hn);
                int single_pin_community_nets = 0;
                if (_remove_single_pin_community_nets) {
                  for (int i = 0; i < (int)_incident_nets[v].size(); ++i) {
                    HyperedgeID he = this->_incident_nets[v][i];
                    if (hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) == 1) {
                      std::swap(this->_incident_nets[v][i],
                                this->_incident_nets[v][single_pin_community_nets]);
                      ++single_pin_community_nets;
                    }
                  }
                }
                this->hypernode(hn).setSinglePinCommunityNets(single_pin_community_nets);
                this->hypernode(hn).setInvalidCommunityNets(_incident_nets[v].size());

                HEAVY_COARSENING_ASSERT([&] {
                  if (_remove_single_pin_community_nets) {
                    size_t single_pin_community_hyperedges = this->hypernode(hn).singlePinCommunityNets();
                    for (size_t i = 0; i < single_pin_community_hyperedges; ++i) {
                      const HyperedgeID he = incident_nets(hn)[i];
                      if (hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) > 1) {
                        LOG << "Hyperedge" << he << "is a non single-pin commnunity hyperedge";
                        return false;
                      }
                    }
                  }
                  return true;
                } (), "There non single-pin community hyperedges in single-pin part of incident nets");

                HEAVY_COARSENING_ASSERT([&] {
                  if (_remove_single_pin_community_nets) {
                    size_t single_pin_community_hyperedges = this->hypernode(hn).singlePinCommunityNets();
                    for (size_t i = single_pin_community_hyperedges; i < incident_nets(hn).size(); ++i) {
                      const HyperedgeID he = incident_nets(hn)[i];
                      if (hypergraph_of_hyperedge(he, hypergraphs).edgeSize(he, community_id) <= 1) {
                        LOG << "Hyperedge" << he << "is a single-pin commnunity hyperedge";
                        return false;
                      }
                    }
                  }
                  return true;
                } (), "There single-pin community hyperedges in non-single-pin part of incident nets");
              }
            });
          });
        });
    group.wait();
  }

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   *
   * The pins of the original hyperedge are sorted in decreasing order of their
   * contraction index. The contraction index of a vertex v is defined as the index
   * of the contraction (u,v) in the contraction history, where v occurs as contraction
   * partner. This is done to fullfil the invariants required by the uncontraction method.
   *
   * Note this function have to be called after parallel community coarsening such
   * that uncontractions can be performed correctly.
   *
   * @params history contraction history
   */
  void removeCommunityHyperedges(const std::vector<HypernodeID>& contraction_index,
                                 const std::vector<Self>& hypergraphs) {
    // The incidence array of a hyperedge is constructed as follows: The first part consists
    // of all enabled pins and the remainder of all invalid pins. The invalid pins in the
    // remainder are sorted in decreasing order of their contraction index.
    tbb::task_group group;
    _arena.execute([&] {
          group.run([&] {
            tbb::parallel_for(0UL, _num_hyperedges, [&](const HyperedgeID& he) {
              Hyperedge& e = this->_hyperedges[he];
              if (e.isInitCommunityHyperedges()) {
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
                for ( ; last_entry >= first_entry; --last_entry) {
                  const HypernodeID pin = this->_incidence_array[last_entry];
                  if (!hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin)) {
                    e.decrementSize();
                  }
                }

                e.deinitializeCommunityHyperedges();

                HEAVY_COARSENING_ASSERT([&] {
                  for (size_t i = e.firstEntry(); i < e.firstInvalidEntry(); ++i) {
                    const HypernodeID& pin = this->_incidence_array[i];
                    if (!hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin)) {
                      LOG << "Hypernode" << pin << "is disabled";
                      return false;
                    }
                  }
                  return true;
                } (), "There are disabled hypernodes in valid part of hyperedge");

                HEAVY_COARSENING_ASSERT([&] {
                  for (size_t i = e.firstInvalidEntry(); i < _hyperedges[he + 1].firstEntry(); ++i) {
                    const HypernodeID& pin = this->_incidence_array[i];
                    if (hypergraph_of_vertex(pin, hypergraphs).nodeIsEnabled(pin)) {
                      LOG << "Hypernode" << pin << "is enabled";
                      return false;
                    }
                  }
                  return true;
                } (), "There are enabled hypernodes in invalid part of hyperedge");
              }
            });
          });
        });
    group.wait();

    CommunityHyperedges community_hyperedges;
    _community_hyperedges = std::move(community_hyperedges);
  }

  // ! Resets the ids of all pins in the incidence array to its original node id
  void resetPinsToOriginalNodeIds(const std::vector<Self>& hypergraphs) {
    tbb::task_group group;
    _arena.execute([&] {
          group.run([&] {
            tbb::parallel_for(0UL, this->_num_pins, [&](const size_t& i) {
              HypernodeID pin = this->_incidence_array[i];
              this->_incidence_array[i] = hypergraph_of_vertex(pin, hypergraphs).hypernode(pin).originalNodeId();
            });
          });
        });
    group.wait();
  }

  // ! Invalidates all disabled hyperedges from the incident nets array of each node
  // ! For further details please take a look at the documentation of uncontraction(...)
  // ! in hypergraph.h
  void invalidateDisabledHyperedgesFromIncidentNets(const std::vector<Self>& hypergraphs) {
    tbb::task_group group;
    _arena.execute([&] {
          group.run([&] {
            tbb::parallel_for(0UL, this->_num_hypernodes, [&](const HypernodeID& id) {
              const HypernodeID hn = get_global_node_id(id);
              invalidateDisabledHyperedgesFromIncidentNets(hn, hypergraphs);
            });
          });
        });
    group.wait();
  }

  // ! Only for assertion
  bool verify_incident_nets_of_hypergraph(const std::vector<Self>& hypergraphs) const {
    for (size_t pos = 0; pos < _incident_nets.size(); ++pos) {
      const HypernodeID& hn = _hypernodes[pos].nodeId();
      for (const HyperedgeID& he : _incident_nets[pos]) {
        const Self& hypergraph_of_he = hypergraph_of_hyperedge(he, hypergraphs);
        const HyperedgeID local_edge_id = get_local_edge_id_of_hyperedge(he);
        ASSERT(local_edge_id < hypergraph_of_he._hyperedges.size());
        const Hyperedge& e = hypergraph_of_he._hyperedges[local_edge_id];
        const auto first = hypergraph_of_he._incidence_array.begin() + e.firstEntry();
        const auto last = first + e.size();
        if (std::find(first, last, hn) == last) {
          LOG << "Hypernode" << hn << "not part of hyperedge" << he << "on numa node" << get_numa_node_of_hyperedge(he);
          hypergraph_of_he.printHyperedgeInfo(he);
          return false;
        }
      }
    }
    return true;
  }

  // ####################### Helper Functions #######################

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_global_node_id(const int node, const size_t node_pos) {
    return (((HyperedgeID)node) << NUMA_NODE_INDENTIFIER) | node_pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_global_edge_id(const int node, const size_t edge_pos) {
    return (((HyperedgeID)node) << NUMA_NODE_INDENTIFIER) | edge_pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_node_id_of_vertex(const HypernodeID u) {
    return ((1UL << NUMA_NODE_INDENTIFIER) - 1) & u;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID u) {
    return (int)(u >> NUMA_NODE_INDENTIFIER);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_hyperedge(const HyperedgeID e) {
    return (int)(e >> NUMA_NODE_INDENTIFIER);
  }

  void printHyperedgeInfo(const HyperedgeID e) const {
    LOG << "Hyperedge:" << e;
    LOG << "Original Size:" << (hyperedge(e + 1).firstEntry() - hyperedge(e).firstEntry());
    if (edgeIsEnabled(e)) {
      LOG << "Current Size:" << hyperedge(e).size();
    }
    for (size_t pos = hyperedge(e).firstEntry(); pos < hyperedge(e + 1).firstEntry(); ++pos) {
      std::cout << _incidence_array[pos] << " ";
    }
    std::cout << std::endl;
  }

 private:
  template <typename HardwareTopology_, typename TBBNumaArena_>
  friend class Hypergraph;

  // ####################### Contract / Uncontract #######################

  // ! Increments the size of an hyperedge to its final size after the batch uncontraction
  void incrementHyperedgeSizeDuringBatchUncontraction(const HyperedgeID e,
                                                      const std::vector<Self>& hypergraphs,
                                                      const kahypar::ds::FastResetFlagArray<>& batch_hypernodes) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if (containsIncidentNet(e)) {
      size_t incidence_array_start = hyperedge(e).firstEntry();
      size_t tmp_size = hyperedge(e).size();
      for ( ; incidence_array_start + tmp_size < hyperedge(e + 1).firstEntry(); ++tmp_size) {
        const size_t pos = incidence_array_start + tmp_size;
        const HypernodeID pin = _incidence_array[pos];
        const HypernodeID original_id = hypergraph_of_vertex(pin, hypergraphs).originalNodeId(pin);
        if (!batch_hypernodes[original_id]) {
          break;
        }
      }
      // Note, since batch_hypernodes is constant within a batch uncontraction, several threads
      // that adapt the size of that hyperedge concurrently should compute the same tmp_size.
      hyperedge(e).setSize(tmp_size);
    }
  }

  // ! Connect hyperedge e to representative hypernode u.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHyperedgeToRepresentative(const HyperedgeID e,
                                                                        const HypernodeID u,
                                                                        Self& hypergraph_of_u) {
    ASSERT(hypergraph_of_u.nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(_node == -1 || _node == get_numa_node_of_hyperedge(e));
    // Hyperedge e does not contain u. Therefore we use the entry of v (i.e. the last entry
    // -- this is ensured by the contract method) in e's edge array to store the information
    // that u is now connected to e and add the edge (u,e) to indicate this conection also from
    // the hypernode's point of view.
    _incidence_array[hyperedge(e).firstInvalidEntry() - 1] = u;
    hypergraph_of_u.incident_nets(u).push_back(e);
  }

  // ! Connect hyperedge e to representative hypernode u.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHyperedgeToRepresentative(const HyperedgeID e,
                                                                        const HypernodeID u,
                                                                        CommunityHyperedge& community_he,
                                                                        Self& hypergraph_of_u) {
    ASSERT(hypergraph_of_u.nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    ASSERT(_node == -1 || _node == get_numa_node_of_hyperedge(e));

    auto& incident_nets_of_u = hypergraph_of_u.incident_nets(u);
    incident_nets_of_u.push_back(e);
    size_t he_pos = incident_nets_of_u.size() - 1;

    // If community hyperedge is not disabled we swap it to valid part of incident nets
    if (!community_he.isDisabled()) {
      size_t invalid_community_nets = hypergraph_of_u.hypernode(u).invalidCommunityNets();
      std::swap(incident_nets_of_u[invalid_community_nets], incident_nets_of_u[he_pos]);
      he_pos = invalid_community_nets;
      hypergraph_of_u.hypernode(u).incrementInvalidCommunityNets();

      // If community hyperedge is a single-pin community hyperedge we swap it to single-pin
      // part of incident nets of u
      if (_remove_single_pin_community_nets && community_he.size() == 1) {
        size_t single_pin_community_nets = hypergraph_of_u.hypernode(u).singlePinCommunityNets();
        std::swap(incident_nets_of_u[single_pin_community_nets], incident_nets_of_u[he_pos]);
        hypergraph_of_u.hypernode(u).incrementSinglePinCommunityNets();
      }
    }
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

    for (int64_t pos = pins_end - 1; pos >= pins_start; --pos) {
      if (_incidence_array[pos] == u) {
        slot_of_u = pos;
        break;
      }
    }

    ASSERT([&] {
          if (slot_of_u == pins_end) {
            LOG << V(pins_start) << V(pins_end) << V(slot_of_u);
            printHyperedgeInfo(he);
            return false;
          }
          return true;
        } (), "Hypernode" << u << "not found in hyperedge" << he);
    ASSERT(_incidence_array[slot_of_u] == u, "Hypernode" << u << "not found in hyperedge" << he);
    _incidence_array[slot_of_u] = v;
  }

  /*!
   * Resets the pin slot containing u back to contain v.
   * If hyperedge he only contained v before the contraction, then the array entry of
   * v in the incidence structure of he is used to store u after the contraction.
   * This method undoes this operation.
   */
  void resetReusedPinSlotToOriginalValue(const HyperedgeID he,
                                         const size_t size,
                                         const HypernodeID u,
                                         const HypernodeID v) {
    ASSERT(hyperedge(he).isDisabled(), "Hyperedge" << he << "is enabled");

    int64_t pins_start = hyperedge(he).firstEntry();
    int64_t pins_end = pins_start + size;
    int64_t slot_of_u = pins_end;

    for (int64_t pos = pins_end - 1; pos >= pins_start; --pos) {
      if (_incidence_array[pos] == u) {
        slot_of_u = pos;
        break;
      }
    }

    ASSERT([&] {
          if (slot_of_u == pins_end) {
            LOG << V(size) << V(pins_start)
                << V(pins_end) << V(slot_of_u);
            printHyperedgeInfo(he);
            return false;
          }
          return true;
        } (), "Hypernode" << u << "not found in hyperedge" << he);
    ASSERT(_incidence_array[slot_of_u] == u, "Hypernode" << u << "not found in hyperedge" << he);
    _incidence_array[slot_of_u] = v;
  }

  // ! Marks all incident nets of hypernode v in bit vector
  // ! Used during uncontraction to detect which hyperedges have to processed (intersection
  // ! of the incident nets of u and v)
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void markAllIncidentNetsOf(const HypernodeID v,
                                                             std::vector<Self>& hypergraphs) {
    // Reset incident nets on each node
    for (size_t node = 0; node < hypergraphs.size(); ++node) {
      hypergraphs[node]._incident_nets_of_v.local().reset();
    }

    for (const HyperedgeID& he : incident_nets(v)) {
      Self& hypergraph_of_he = hypergraph_of_hyperedge(he, hypergraphs);
      HyperedgeID local_he_id = get_local_edge_id_of_hyperedge(he);
      ASSERT(local_he_id < hypergraph_of_he._num_hyperedges);
      hypergraph_of_he._incident_nets_of_v.local().set(local_he_id, true);
    }
  }

  // ! Checks whether hyperedge e is incident to hypernode v
  // ! Note, in order that function works correctly, markAllIncidentNetsOf(v) have to
  // ! called before
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool containsIncidentNet(const HyperedgeID e) {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of this numa node");
    ASSERT(local_id <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    kahypar::ds::FastResetFlagArray<>& local_incident_nets_of_v = _incident_nets_of_v.local();
    return local_incident_nets_of_v[local_id];
  }

  // ! Returns the uncontraction case, if we are in sequential uncontraction mode
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE UncontractionCase get_uncontraction_case(const HyperedgeID he,
                                                                           const size_t size,
                                                                           const HypernodeID v) const {
    size_t incidence_array_start = hyperedge(he).firstEntry();
    if (incidence_array_start + size < hyperedge(he + 1).firstEntry() &&
        _incidence_array[incidence_array_start + size] == v) {
      return UncontractionCase::CASE_1;
    } else {
      return UncontractionCase::CASE_2;
    }
  }

  // ! Returns the uncontraction case, if we are in batch uncontraction mode
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE UncontractionCase get_uncontraction_case(const HyperedgeID he,
                                                                           const size_t size,
                                                                           const HypernodeID v,
                                                                           const std::vector<Self>& hypergraphs,
                                                                           const kahypar::ds::FastResetFlagArray<>& batch_hypernodes) const {
    ASSERT(size > 0);
    const int64_t incidence_array_start = hyperedge(he).firstEntry();
    const int64_t invalid_part_start = incidence_array_start + size;
    const int64_t incidence_array_end = hyperedge(he + 1).firstEntry();
    bool is_forward = invalid_part_start < incidence_array_end &&
                      // Checks wheter the first vertex in the invalid part is part of the batch or not
                      batch_hypernodes[hypergraph_of_vertex(_incidence_array[invalid_part_start],
                                                            hypergraphs).originalNodeId(_incidence_array[invalid_part_start])];

    if (is_forward) {
      // In that case, the size of the hyperedge is the same as before the batch
      // uncontraction.
      for (int64_t pos = invalid_part_start; pos < incidence_array_end; ++pos) {
        const HypernodeID pin = _incidence_array[pos];
        const HypernodeID original_id = hypergraph_of_vertex(pin, hypergraphs).originalNodeId(pin);
        if (!batch_hypernodes[original_id]) {
          return UncontractionCase::CASE_2;
        } else if (pin == v) {
          // Indicates that the batch uncontraction function can adapt the size
          // of the hyperedge to its final size after the batch uncontraction
          return UncontractionCase::CASE_1_FORWARD;
        }
      }
    } else {
      // In that case, the size of the hyperedge is already increased to its size after
      // the batch uncontraction (by an other thread processing the same hyperedge).
      // In that case, we have to look backward from the current size of the hyperedge
      // to decide if we are in uncontraction case 1 or 2.
      for (int64_t pos = invalid_part_start - 1; pos >= incidence_array_start; --pos) {
        const HypernodeID pin = _incidence_array[pos];
        const HypernodeID original_id = hypergraph_of_vertex(pin, hypergraphs).originalNodeId(pin);
        if (!batch_hypernodes[original_id]) {
          return UncontractionCase::CASE_2;
        } else if (pin == v) {
          return UncontractionCase::CASE_1_BACKWARD;
        }
      }
    }

    return UncontractionCase::CASE_2;
  }

  // ####################### Remove / Restore Hyperedges #######################

  // ! Removes hyperedge e from the incident nets of vertex hn
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

  // ! Removes hyperedge e from the incident nets of vertex hn, if community hyperedges are
  // ! initialized.
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
    for ( ; pos >= 0; --pos) {
      if (incident_nets_of_hn[pos] == he) {
        break;
      }
    }
    ASSERT(incident_nets_of_hn[pos] == he, "Hyperedge" << he << "not found");

    // If hyperedge is in single-pin part of incident nets we first swap it
    // to valid part incident nets
    int64_t single_pin_community_nets = hypernode(hn).singlePinCommunityNets();
    ASSERT(pos != -1);
    if (pos < single_pin_community_nets) {
      ASSERT(single_pin_community_nets > 0);
      std::swap(incident_nets_of_hn[pos], incident_nets_of_hn[single_pin_community_nets - 1]);
      pos = single_pin_community_nets - 1;
      hypernode(hn).decrementSinglePinCommunityNets();
    }

    // Afterwards, we swap incident nets to invalid part ...
    std::swap(incident_nets_of_hn[pos], incident_nets_of_hn[invalid_community_nets - 1]);
    pos = invalid_community_nets - 1;
    hypernode(hn).decrementInvalidCommunityNets();

    if (!invalidate_only) {
      // ... and if hyperedge should be removed from incident nets, we swap it to the end
      // and pop back.
      std::swap(incident_nets_of_hn[pos], incident_nets_of_hn.back());
      incident_nets_of_hn.pop_back();
    } else {
      ASSERT(hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled());
    }

    HEAVY_COARSENING_ASSERT([&] {
          size_t invalid_community_nets = hypernode(hn).invalidCommunityNets();
          for (size_t i = 0; i < invalid_community_nets; ++i) {
            const HyperedgeID he = incident_nets_of_hn[i];
            if (hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled()) {
              LOG << "HE" << he << "should be in invalid part of incident nets of HN" << hn;
              return false;
            }
          }
          return true;
        } (), "There is an invalidated community hyperedge in valid part of incident nets");

    HEAVY_COARSENING_ASSERT([&] {
          size_t invalid_community_nets = hypernode(hn).invalidCommunityNets();
          for (size_t i = invalid_community_nets; i < incident_nets_of_hn.size(); ++i) {
            const HyperedgeID he = incident_nets_of_hn[i];
            if (!hypergraph_of_hyperedge(he, hypergraphs).community_hyperedge(he, community_id).isDisabled()) {
              LOG << "HE" << he << "should be in valid part of incident nets of HN" << hn;
              return false;
            }
          }
          return true;
        } (), "There is an valid community hyperedge in invalid part of incident nets");
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID he,
                                                                     const HypernodeID hn) {
    HEAVY_REFINEMENT_ASSERT(std::count(incident_nets(hn).begin() + hypernode(hn).invalidIncidentNets(),
                                    incident_nets(hn).end(), he) == 0,
                        "HN" << hn << "is already connected to HE" << he);
    incident_nets(hn).push_back(he);
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernodeFromInvalidPart(const HyperedgeID he,
                                                                                    const HypernodeID hn) {
    size_t invalid_incident_nets = hypernode(hn).invalidIncidentNets();
    HEAVY_REFINEMENT_ASSERT(std::count(incident_nets(hn).begin() + invalid_incident_nets,
                                       incident_nets(hn).end(), he) == 0,
                            "HN" << hn << "is already connected to HE" << he);

    auto& incident_nets_of_hn = incident_nets(hn);
    size_t slot_of_he = invalid_incident_nets;
    for (size_t pos = 0; pos < invalid_incident_nets; ++pos) {
      if (incident_nets_of_hn[pos] == he) {
        slot_of_he = pos;
        break;
      }
    }

    if (slot_of_he < invalid_incident_nets) {
      ASSERT(incident_nets_of_hn[slot_of_he] == he);
      std::swap(incident_nets_of_hn[slot_of_he],
                incident_nets_of_hn[invalid_incident_nets - 1]);
      hypernode(hn).decrementInvalidIncidentNets();
    } else {
      incident_nets_of_hn.push_back(he);
    }
  }

  // ! Invalidates all disabled hyperedges of the incident nets array of hypernode v
  // ! For a more detailed explanation see hypergraph.h
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void invalidateDisabledHyperedgesFromIncidentNets(const HypernodeID hn,
                                                                                    const std::vector<Self>& hypergraphs) {
    size_t invalid_incident_nets = 0;
    size_t incident_nets_end = incident_nets(hn).size();
    for (size_t pos = 0; pos < incident_nets_end; ++pos) {
      const HyperedgeID he = incident_nets(hn)[pos];
      if (!hypergraph_of_hyperedge(he, hypergraphs).edgeIsEnabled(he)) {
        std::swap(incident_nets(hn)[pos], incident_nets(hn)[invalid_incident_nets++]);
      }
    }
    hypernode(hn).setInvalidIncidentNets(invalid_incident_nets);
  }

  // ####################### Partition Information #######################

  // ! Decrements the number of incident cut hyperedges of hypernode u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void decrementIncidentNumCutHyperedges(const HypernodeID u) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    --_atomic_hn_data[get_local_node_id_of_vertex(u)].num_incident_cut_hes;
  }

  // ! Increments the number of incident cut hyperedges of hypernode u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void incrementIncidentNumCutHyperedges(const HypernodeID u) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    ++_atomic_hn_data[get_local_node_id_of_vertex(u)].num_incident_cut_hes;
  }

  // ! Decrements the number of pins of hyperedge he in block id
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID decrementPinCountInPart(const HyperedgeID he, const PartitionID part_id) {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge" << he << "is disabled");
    /*ASSERT(pinCountInPart(he, id) > 0,
           "HE" << he << ": pin_count[" << id << "]=" << pinCountInPart(he, id)
                << "edgesize=" << edgeSize(he));*/
    ASSERT(part_id < _k && part_id != kInvalidPartition, "Part ID" << part_id << "out of bounds!");
    const HyperedgeID local_id = get_local_edge_id_of_hyperedge(he);
    const size_t offset = local_id * _k + part_id;
    ASSERT(offset < _pins_in_part.size());
    const HypernodeID pin_count_after = --_pins_in_part[offset];
    const bool connectivity_decreased = pin_count_after == 0;
    if (connectivity_decreased) {
      _connectivity_sets.remove(local_id, part_id);
    }
    return pin_count_after;
  }

  // ! Increments the number of pins of hyperedge he in block id
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID incrementPinCountInPart(const HyperedgeID he, const PartitionID part_id) {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge" << he << "is disabled");
    /*ASSERT(pinCountInPart(he, id) <= edgeSize(he),
           "HE" << he << ": pin_count[" << id << "]=" << pinCountInPart(he, id)
                << "edgesize=" << edgeSize(he));
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << "out of bounds!");*/
    const HyperedgeID local_id = get_local_edge_id_of_hyperedge(he);
    const size_t offset = local_id * _k + part_id;
    ASSERT(offset < _pins_in_part.size());
    const HypernodeID pin_count_after = ++_pins_in_part[offset];
    const bool connectivity_increased = pin_count_after == 1;
    if (connectivity_increased) {
      _connectivity_sets.add(local_id, part_id);
    }
    return pin_count_after;
  }

  // ####################### Hypernode Information #######################

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID get_global_node_id() {
    const HypernodeID local_node_id = _next_node_id++;
    const HypernodeID numa_node = static_cast<HypernodeID>(std::max(_node, 0));
    return (numa_node << NUMA_NODE_INDENTIFIER) | local_node_id;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID get_global_node_id(const HypernodeID local_id) const {
    const HypernodeID numa_node = static_cast<HypernodeID>(std::max(_node, 0));
    return (numa_node << NUMA_NODE_INDENTIFIER) | local_id;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Self& hypergraph_of_vertex(const HypernodeID u,
                                                                          const std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_vertex(u);
    ASSERT(node < (int)hypergraph.size());
    return hypergraph[node];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Self& hypergraph_of_vertex(const HypernodeID u,
                                                                    std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_vertex(u);
    ASSERT(node < (int)hypergraph.size());
    return hypergraph[node];
  }

  // ! Accessor for hypernode-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(_node == -1 || get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const StreamingHypergraph&>(*this).hypernode(u));
  }

  // ! Accessor to the incident nets array of hypernode u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(_node == -1 || get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_nets[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) {
    return const_cast<parallel::scalable_vector<HyperedgeID>&>(static_cast<const StreamingHypergraph&>(*this).incident_nets(u));
  }

  // ####################### Hyperedge Information #######################

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeID get_global_edge_id(const size_t edge_pos) const {
    const HyperedgeID numa_node = static_cast<HypernodeID>(std::max(_node, 0));
    return (numa_node << NUMA_NODE_INDENTIFIER) | edge_pos;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_edge_id_of_hyperedge(const HyperedgeID e) {
    return ((1UL << NUMA_NODE_INDENTIFIER) - 1) & e;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Self& hypergraph_of_hyperedge(const HyperedgeID e,
                                                                       std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_hyperedge(e);
    ASSERT(node < (int)hypergraph.size());
    return hypergraph[node];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Self& hypergraph_of_hyperedge(const HyperedgeID e,
                                                                             const std::vector<Self>& hypergraph) {
    int node = get_numa_node_of_hyperedge(e);
    ASSERT(node < (int)hypergraph.size());
    return hypergraph[node];
  }

  // ! Accessor for hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
    // <= instead of < because of sentinel
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const StreamingHypergraph&>(*this).hyperedge(e));
  }

  // ####################### Community Hyperedge Information #######################

  // ! Accessor for community hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const CommunityHyperedge& community_hyperedge(const HyperedgeID e, const PartitionID community_id) const {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hyperedges, "Hyperedge" << e << "does not exist");

    size_t community_hyperedge_position = find_position_of_community_hyperedge(e, community_id);
    ASSERT(community_hyperedge_position < _community_hyperedges[local_id].size(),
           "Community hyperedge" << e << "with community id" << community_id << "not found");
    return _community_hyperedges[local_id][community_hyperedge_position];
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t find_position_of_community_hyperedge(const HyperedgeID e, const PartitionID community_id) const {
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(_node == -1 || get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hyperedges, "Hyperedge" << e << "does not exist");

    size_t pos = 0;
    for ( ; pos < _community_hyperedge_ids[local_id].size(); ++pos) {
      if (_community_hyperedge_ids[local_id][pos] == community_id) {
        break;
      }
    }

    ASSERT(pos < _community_hyperedges[local_id].size(),
           "Community hyperedge" << e << "with community id" << community_id << "not found");
    return pos;
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CommunityHyperedge& community_hyperedge(const HyperedgeID e, const PartitionID community_id) {
    return const_cast<CommunityHyperedge&>(static_cast<const StreamingHypergraph&>(*this).community_hyperedge(e, community_id));
  }

  const int _node;
  // ! number of blocks
  const PartitionID _k;
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
  std::vector<parallel::scalable_vector<PartitionID> > _community_hyperedge_ids;
  // ! Community Hyperedges
  CommunityHyperedges _community_hyperedges;

  // ! Will be used during uncontraction to mark
  // ! all incident nets of contraction partner v
  ThreadLocalFastResetFlagArray _incident_nets_of_v;

  // ! For each hypernode, _atomic_hn_data stores the corresponding block of the vertex
  // ! and the number of incident cut hyperedges
  parallel::scalable_vector<AtomicHypernodeData> _atomic_hn_data;
  // ! For each hyperedge and each block, _pins_in_part stores the number of pins in that block
  parallel::scalable_vector<HypernodeAtomic> _pins_in_part;
  // ! For each hyperedge, _connectivity_sets stores the connectivity and the set of block ids
  // ! which the pins of the hyperedge belongs to
  ConnectivitySets _connectivity_sets;

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
}  // namespace ds
}  // namespace mt_kahypar
