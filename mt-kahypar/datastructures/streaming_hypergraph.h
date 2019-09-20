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

#include "tbb/task_scheduler_observer.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"

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
        _id(0),
        _original_id(0),
        _weight(1),
        _community_id(kInvalidPartition),
        _valid(false) { }

      Hypernode(const HypernodeID id,
                const HypernodeID original_id,
                const HypernodeWeight weight) :
        _id(id),
        _original_id(original_id),
        _weight(weight),
        _community_id(kInvalidPartition),
        _valid(true) { }

      HypernodeID nodeId() const {
        return _id;
      }

      HypernodeID originalNodeId() const {
        return _original_id;
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
        return _community_id;
      }

      void setCommunityID(const PartitionID community_id) {
        _community_id = community_id;
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
      // ! Hypernode weight
      HyperedgeWeight _weight;
      // ! Community id
      PartitionID _community_id;
      // ! Flag indicating whether or not the element is active.
      bool _valid;
  };

  class Hyperedge {

    public:
      using IDType = HyperedgeID;

      Hyperedge() :
        _begin(0),
        _size(0),
        _weight(1),
        _hash(kEdgeHashSeed),
        _valid(false) { }

      Hyperedge(const size_t begin, 
                const size_t size,
                const HyperedgeWeight weight) :
        _begin(begin),
        _size(size),
        _weight(weight),
        _hash(kEdgeHashSeed),
        _valid(true) { }

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

 public:
  explicit StreamingHypergraph(const int node) :
    _node(node),
    _arena(TBBNumaArena::instance().numa_task_arena(node)),
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
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
    _num_hypernodes(other._num_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_pins(other._num_pins),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
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
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).originalNodeId();
  }

  // ! Returns a for-each iterator-pair to loop over the set of incident hyperedges of hypernode u.
  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return std::make_pair(incident_nets(u).cbegin(),
                          incident_nets(u).cend());
  }

  // ! Returns a const reference to the incident net vector of hypernode u
  const parallel::scalable_vector<HyperedgeID>& incidentNets(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return incident_nets(u);
  }

  // ! Returns a for-each iterator-pair to loop over the set pins of hyperedge e.
  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return std::make_pair(_incidence_array.cbegin() + hyperedge(e).firstEntry(),
                          _incidence_array.cbegin() + hyperedge(e).firstInvalidEntry());
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

  bool uncontract(const HypernodeID u, const HypernodeID v, 
                  const HyperedgeID e, const size_t incident_nets_pos,
                  std::vector<Self>& hypergraphs) {
    using std::swap;
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    if ( containsIncidentNet(e) ) {
      // ... then we have to do some kind of restore operation.
      if (_incidence_array[hyperedge(e).firstInvalidEntry()] == v &&
          hyperedge(e).firstInvalidEntry() < hyperedge(e + 1).firstEntry()) {
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
        int node = get_numa_node_of_vertex(u);
        ASSERT(node < (int) hypergraphs.size());
        ASSERT(e == hypergraphs[node].incident_nets(u)[incident_nets_pos]);

        size_t incident_nets_end = hypergraphs[node].incident_nets(u).size();
        swap(hypergraphs[node].incident_nets(u)[incident_nets_pos],
             hypergraphs[node].incident_nets(u)[incident_nets_end - 1]);
        hypergraphs[node].incident_nets(u).pop_back();

        DBG << "resetting reused Pinslot of HE" << e << "from" << u << "to" << v;
        resetReusedPinSlotToOriginalValue(e, u, v);

        // TODO(heuer): Increment pin count in part
        return true;
      }
    }
    return false;
  }

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

  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    hyperedge(e).setWeight(weight);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return incident_nets(u).size();
  }

  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).size();
  }

  size_t edgeHash(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).hash();
  }

  PartitionID communityID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityID();
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

  HypernodeID streamHypernode(HypernodeID original_id, HypernodeWeight weight) {
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(sched_getcpu()) == _node);
    HypernodeID node_id = get_global_node_id();
    _hypernode_stream.stream(node_id, original_id, weight);
    return node_id;
  }

  void streamHyperedge(const std::vector<HypernodeID>& hyperedge, const HyperedgeWeight& weight) {
    int cpu_id = sched_getcpu();
    // Make sure calling process is part of correct numa node
    ASSERT(HardwareTopology::instance().numa_node_of_cpu(cpu_id) == _node,
      "Expected that assigned cpu is on numa node" << _node << ", but was CPU" << cpu_id 
       << "is on node" << HardwareTopology::instance().numa_node_of_cpu(cpu_id));
    _hyperedge_stream.stream(_pin_stream.size(cpu_id), hyperedge.size(), weight);
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
    TBBNumaArena::instance().wait(_node, group);
    // Emplace Back Sentinel
    _hyperedges.emplace_back(_incidence_array.size(), 0, 0);
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
    TBBNumaArena::instance().wait(_node, group);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("sort_and_remap_node_ids", "Sort and Remap Nodes",
      "initialize_numa_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());

    // Compute Total Hypergraph Weight
    start = std::chrono::high_resolution_clock::now();
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
    TBBNumaArena::instance().wait(_node, group);
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
              const int node = get_numa_node_of_vertex(pin);
              ASSERT(node < (int) hypergraphs.size());
              hypergraphs[node].streamIncidentNet(pin, he_id);
              // Initialize edge hash
              he.hash() += kahypar::math::hash(pin);
            }
          }
        });
      });
    });
    TBBNumaArena::instance().wait(_node, group);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("stream_incident_nets", "Stream Incident Nets",
      "initialize_numa_hypernodes", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());


    _hypernode_stream.clear();
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

  // ! Only for assertion
  bool verify_incident_nets_of_hypergraph(const std::vector<Self>& hypergraphs) const {
    for ( size_t pos = 0; pos < _incident_nets.size(); ++pos ) {
      const HypernodeID& hn = _hypernodes[pos].nodeId();
      for ( const HyperedgeID& he : _incident_nets[pos] ) {
        const int node = get_numa_node_of_hyperedge(he);
        const HyperedgeID local_edge_id = get_local_edge_id_of_hyperedge(he);
        ASSERT(node < (int) hypergraphs.size());
        ASSERT(local_edge_id < hypergraphs[node]._hyperedges.size());
        const Hyperedge& e = hypergraphs[node]._hyperedges[local_edge_id];
        const auto first = hypergraphs[node]._incidence_array.begin() + e.firstEntry();
        const auto last = first + e.size();
        if ( std::find(first, last, hn) == last ) {
          LOG << "Hypernode" << hn << "not part of hyperedge" << he << "on numa node" << node;
          return false;
        }
      }
    }
    return true;
  }

  // ! Only for testing
  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
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
    // Hyperedge e does not contain u. Therefore we use the entry of v (i.e. the last entry
    // -- this is ensured by the contract method) in e's edge array to store the information
    // that u is now connected to e and add the edge (u,e) to indicate this conection also from
    // the hypernode's point of view.
    _incidence_array[hyperedge(e).firstInvalidEntry() - 1] = u;
    hypergraph_of_u.incident_nets(u).push_back(e);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void markAllIncidentNetsOf(const HypernodeID v,
                                                             std::vector<Self>& hypergraphs) {
    // Reset incident nets on each node
    for ( size_t node = 0; node < hypergraphs.size(); ++node ) {
      hypergraphs[node]._incident_nets_of_v.reset();
    }

    for ( const HyperedgeID& he : incident_nets(v) ) {
      int node = get_numa_node_of_hyperedge(he);
      HyperedgeID local_he_id = get_local_edge_id_of_hyperedge(he);
      ASSERT(node < (int) hypergraphs.size());
      ASSERT(local_he_id < _num_hyperedges);
      hypergraphs[node]._incident_nets_of_v.set(local_he_id, true);
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
  void resetReusedPinSlotToOriginalValue(const HyperedgeID he, const HypernodeID u, const HypernodeID v) {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge" << he << "is disabled");

    size_t pins_start = hyperedge(he).firstEntry();
    size_t pins_end = hyperedge(he).firstInvalidEntry();
    size_t slot_of_u = pins_end;

    for ( size_t pos = pins_start; pos < pins_end; ++pos ) {
      if ( _incidence_array[pos] == u ) {
        slot_of_u = pos;
        break;
      }
    }
    ASSERT(slot_of_u < pins_end, "Hypernode" << u << "not found in hyperedge" << he);
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

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID he,
                                                                     const HypernodeID hn) {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode" << hn << "is disabled");
    ASSERT(std::count(incident_nets(hn).begin(),
                      incident_nets(hn).end(), he)
            == 0,
            "HN" << hn << "is already connected to HE" << he);
    incident_nets(hn).push_back(he);
    // TODO(heuer): increment pin count in part
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

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID u) {
    return (int) (u >> NUMA_NODE_INDENTIFIER); 
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_node_id_of_vertex(const HypernodeID u) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & u; 
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_hyperedge(const HyperedgeID e) {
    return (int) (e >> NUMA_NODE_INDENTIFIER); 
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_edge_id_of_hyperedge(const HyperedgeID e) {
    return ( (1UL << NUMA_NODE_INDENTIFIER) - 1 ) & e; 
  }

  // ! Accessor for hypernode-related information
  const Hypernode & hypernode(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[local_id];
  }

  // ! Accessor for hyperedge-related information
  const Hyperedge & hyperedge(const HyperedgeID e) const {
    // <= instead of < because of sentinel
    HypernodeID local_id = get_local_edge_id_of_hyperedge(e);
    ASSERT(get_numa_node_of_hyperedge(e) == _node, "Hyperedge" << e << "is not part of numa node" << _node);
    ASSERT(local_id <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[local_id];
  }

  const parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) const {
    HypernodeID local_id = get_local_node_id_of_vertex(u);
    ASSERT(get_numa_node_of_vertex(u) == _node, "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_nets[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  Hypernode & hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const StreamingHypergraph&>(*this).hypernode(u));
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  Hyperedge & hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const StreamingHypergraph&>(*this).hyperedge(e));
  }

  parallel::scalable_vector<HyperedgeID>& incident_nets(const HypernodeID u) {
    return const_cast<parallel::scalable_vector<HyperedgeID>&>(static_cast<const StreamingHypergraph&>(*this).incident_nets(u));
  }

  const int _node;
  // ! task arena for numa node
  tbb::task_arena& _arena;

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