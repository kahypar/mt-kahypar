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

class DynamicHypergraph {

  static constexpr bool debug = false;
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

  struct Uncontraction {
    HypernodeID subtree_size;
    Memento memento;
  };

  struct UncontractionComparator {
    bool operator()(const Uncontraction& lhs, const Uncontraction& rhs) {
      return lhs.subtree_size < rhs.subtree_size ||
        (lhs.subtree_size == rhs.subtree_size && lhs.memento.u < rhs.memento.u) ||
        (lhs.subtree_size == rhs.subtree_size && lhs.memento.u == rhs.memento.u && lhs.memento.v < rhs.memento.v);
    }
  };

  using Batch = parallel::scalable_vector<Memento>;
  using BatchVector = parallel::scalable_vector<Batch>;
  using VersionedBatchVector = parallel::scalable_vector<BatchVector>;

  /**
   * Container used to create a vector of batches in parallel.
   * A batch is a vector of mementos (uncontractions) that are uncontracted in
   * parallel. A batch has a predefined maximum size, but it is also possible
   * that less mementos are part of the batch. Several threads can push mementos
   * into the vector in parallel. Once the batch reaches the maximum batch size
   * we create a new empty batch via an unique lock.
   * Note, the data structure is not implemented in the most efficient way and we do
   * not expect good speedups. However, the number of elements inserted into the data
   * structure is less than the number of hypernodes of the hypergraph and it should
   * suffice to be neglible to the overall algorithm.
   */
  class ConcurrentBatchVector {

    public:
      ConcurrentBatchVector(const size_t max_batch_size) :
        _rw_mutex(),
        _max_batch_size(max_batch_size),
        _active_batch_idx(0),
        _current_idx(0),
        _batches(1, parallel::scalable_vector<Memento>(max_batch_size)) { }

      // ! Returns the final batch vector
      // ! Note, data structure is afterwards not usable any more.
      BatchVector&& get() {
        // In case last batch is smaller than the maximum allowed batch
        // size, we adapt the size of the last batch
        if ( _current_idx < _max_batch_size ) {
          _batches.back().resize(_current_idx);
        }
        std::reverse(_batches.begin(), _batches.end());
        return std::move(_batches);
      }

      // ! Pushes a memento into the current batch. If there is no
      // ! free slot in the current batch an new batch is requested.
      // ! Returns the index of the batch in which we inserted the corresponding
      // ! memento
      size_t push(const Memento& memento) {
        size_t batch_idx = 0;
        {
          std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
          batch_idx = _active_batch_idx;
          const size_t idx = _current_idx++;
          if ( idx < _max_batch_size ) {
            ASSERT(batch_idx < _batches.size());
            _batches[batch_idx][idx] = memento;
            return batch_idx;
          }
        }

        // In case we were unable to insert the memento into the current batch
        // we request a new empty batch and try to insert it again.
        addEmptyBatch(batch_idx);
        return push(memento);
      }

      // ! Returns the index of the current batch
      size_t activeBatchIndex() const {
        return _active_batch_idx;
      }

      // ! Creates a new empty batch via an unique lock
      // ! The expected_batch_idx is the index of the batch on which the calling
      // ! thread was unable to insert its memento.
      void addEmptyBatch(const size_t expected_batch_idx) {
        std::unique_lock<std::shared_timed_mutex> lock(_rw_mutex, std::defer_lock);
        while ( !lock.try_lock() ) {
          // We try to acquire an unique lock and if we fail we check if
          // an other thread added a new batch in the meantime.
          if ( expected_batch_idx < _active_batch_idx ) {
            return;
          }
        }

        // If we acquire the unqiue lock and the expected batch is equal
        // with the current active batch, we add a new empty batch. Otherwise,
        // an other thread already added a new empty batch.
        if ( expected_batch_idx == _active_batch_idx ) {
          // Resize current active batch if it contains less elements
          // than the maximum allowed batch size
          if ( _current_idx < _max_batch_size ) {
            _batches[_active_batch_idx].resize(_current_idx);
          }

          _batches.emplace_back();
          _batches[++_active_batch_idx].resize(_max_batch_size);
          _current_idx = 0;
        }
      }

    private:
      // ! Read-Write lock to protect batch vector, if we add a new empty batch
      std::shared_timed_mutex _rw_mutex;
      // ! Maximum allowed size of a batch
      const size_t _max_batch_size;
      // ! Index of the current active batch
      size_t _active_batch_idx;
      // ! Atomic counter used to assign a memento its index in the current active batch
      std::atomic<size_t> _current_idx;
      // ! Batch vector
      BatchVector _batches;
  };

  using IncidenceArray = Array<HypernodeID>;
  using IncidentNets = parallel::scalable_vector<IncidentNetVector<HyperedgeID>>;
  using OwnershipVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<bool>>;
  using ThreadLocalHyperedgeVector = tbb::enumerable_thread_specific<parallel::scalable_vector<HyperedgeID>>;
  using ThreadLocalBitset = tbb::enumerable_thread_specific<parallel::scalable_vector<bool>>;

 public:
  static constexpr bool is_static_hypergraph = true;
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
    _total_degree(0),
    _total_weight(0),
    _version(0),
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
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _version(other._version),
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
    _community_support(std::move(other._community_support)) { }

  DynamicHypergraph & operator= (DynamicHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _max_edge_size = other._max_edge_size;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _version = other._version;
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
    _removable_single_pin_and_parallel_nets = std::move(_removable_single_pin_and_parallel_nets);
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
    ERROR("numGraphEdges() is not supported in dynamic hypergraph");
    return 0;
  }

  HyperedgeID numNonGraphEdges() const {
    ERROR("numNonGraphEdges() is not supported in dynamic hypergraph");
    return initialNumEdges();
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

  HyperedgeID graphEdgeID(const HyperedgeID) const {
    ERROR("graphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
  }

  HyperedgeID nonGraphEdgeID(const HyperedgeID) const {
    ERROR("nonGraphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
  }

  HypernodeID graphEdgeHead(const HyperedgeID, const HypernodeID) const {
    ERROR("nonGraphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
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
  void contract(const HypernodeID v,
                const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
    ASSERT(_contraction_tree.parent(v) != v, "No contraction registered for hypernode" << v);

    HypernodeID x = _contraction_tree.parent(v);
    HypernodeID y = v;
    ContractionResult res = ContractionResult::CONTRACTED;
    // We perform all contractions registered in the contraction tree
    // as long as there are no pending contractions (_hn_ref_count[y] == 0
    // is equivalent with no pending contractions)
    while ( x != y && res != ContractionResult::PENDING_CONTRACTIONS) {
      // Perform Contraction
      res = contract(x, y, max_node_weight);
      y = x;
      x = _contraction_tree.parent(y);
    }
  }

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   */
  void uncontract(const Batch& batch) {
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
          uncontractHyperedge(memento.u, memento.v, he, removable_incident_nets_of_u);
          releaseHyperedge(he);
        } else {
          failed_hyperedge_uncontractions.push_back(he);
        }
      }

      // Perform uncontractions on which we failed to acquire ownership on the first try
      for ( const HyperedgeID& he : failed_hyperedge_uncontractions ) {
        acquireHyperedge(he);
        uncontractHyperedge(memento.u, memento.v, he, removable_incident_nets_of_u);
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
  VersionedBatchVector createBatchUncontractionHierarchy(const TaskGroupID task_group_id,
                                                         const size_t batch_size,
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
    for ( size_t version = 0; version < num_versions; ++version ) {
      versioned_batches[version] =
        createBatchUncontractionHierarchyForVersion(task_group_id, batch_size, version);
      if ( version > 1 ) {
        batch_sizes_prefix_sum[version] =
          batch_sizes_prefix_sum[version - 1] + versioned_batches[version].size();
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
    return createBatchUncontractionHierarchy(TBBNumaArena::GLOBAL_TASK_GROUP, batch_size, true);
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

    // Adds all restored hyperedges as incident net to its contained pins.
    // In the previous step we inserted all pins together with its hyperedges
    // into a bucket data structure. All hyperedges of the same pin are placed
    // within the same bucket and can be processed sequentially here without
    // locking.
    tbb::parallel_for(0UL, incident_net_map.numBuckets(), [&](const size_t bucket) {
      auto& incident_net_bucket = incident_net_map.getBucket(bucket);
      std::sort(incident_net_bucket.begin(), incident_net_bucket.end(),
        [&](const Memento& lhs, const Memento& rhs) {
          return lhs.u < rhs.u || ( lhs.u == rhs.u || lhs.v < rhs.v);
        });

      // No locking required since vertex u can only occur in one bucket
      for ( const Memento& memento : incident_net_bucket ) {
        const HypernodeID u = memento.u;
        const HyperedgeID he = memento.v;
        _incident_nets[u].push_back(he);
      }

      incident_net_map.free(bucket);
    });

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
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

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
        hypergraph._incident_nets[hn].resize(_incident_nets[hn].size());
        hypergraph._acquired_hns[hn] = _acquired_hns[hn];
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
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

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

    utils::MemoryTreeNode* contraction_tree_node = parent->addChild("Contraction Tree");
    _contraction_tree.memoryConsumption(contraction_tree_node);
    utils::MemoryTreeNode* community_support_node = parent->addChild("Community Support");
    _community_support.memoryConsumption(community_support_node);
  }

 private:
  friend class DynamicHypergraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;

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
    //  3.) Resulting node weight is less or equal than a predefined upper bound
    const bool contraction_partner_valid = nodeIsEnabled(v) && _contraction_tree.pendingContractions(v) == 0;
    const bool less_or_equal_than_max_node_weight =
      hypernode(u).weight() + hypernode(v).weight() <= max_node_weight;
    if ( contraction_partner_valid && less_or_equal_than_max_node_weight ) {
      ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled!");
      hypernode(u).setWeight(nodeWeight(u) + nodeWeight(v));
      hypernode(v).disable();
      releaseHypernode(u);
      releaseHypernode(v);

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

      acquireHypernode(u);
      _contraction_tree.unregisterContraction(u, v);
      releaseHypernode(u);
      return ContractionResult::CONTRACTED;
    } else {
      ContractionResult res = ContractionResult::PENDING_CONTRACTIONS;
      if ( !less_or_equal_than_max_node_weight ) {
        _contraction_tree.unregisterContraction(u, v, true /* failed */);
        res = ContractionResult::WEIGHT_LIMIT_REACHED;
      }
      releaseHypernode(u);
      releaseHypernode(v);
      return res;
    }
  }

  // ! Uncontracts u and v in hyperedge he.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void uncontractHyperedge(const HypernodeID u,
                                                           const HypernodeID v,
                                                           const HyperedgeID he,
                                                           parallel::scalable_vector<bool>& removable_incident_nets_of_u) {
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
      ASSERT(hypernode(pin).batchIndex() <= batch_index);
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

  /**
   * Computes a batch uncontraction hierarchy for a specific version of the hypergraph.
   * A batch is a vector of mementos (uncontractions) that are uncontracted in parallel.
   * Each time we perform single-pin and parallel net detection we create a new version of
   * the hypergraph.
   * A batch of uncontractions must satisfy the condition that all representative of uncontractions
   * are active vertices of the hypergraph. Once a batch is uncontracted all contraction partners
   * of uncontractions within the batch become active. Our algorithm works in two stages. First,
   * we create for each root of the contraction tree a local batch uncontraction order in parallel. Note,
   * that the batches of the local uncontraction hierarchy are not constraint by the maximum batch
   * size. In a second step, we merge the local batch uncontraction orders into one vector of batches
   * that respect the maximum allowed batch size. Note that uncontractions originate from different roots
   * can be interleaved arbitrary as long as we do not add two uncontractions from different local batches
   * of the same root into the same global batch. Our merging step aims to interleave the local batches
   * of different roots in a round-robin fashion. There might be also different merging strategies and
   * our implementation might be not the most efficient parallel implementation. However, the reason
   * of interleaving the batches in a round-robin fashion is that this lead to batches which contains
   * a maximum number of different representatives. This might be beneficial for local searches, because
   * it should lead to mostly indepedent local searches in the hypergraph (meaning that the uncontracted
   * vertices have a certain distance to each other and parallel local searches do not disturb each other).
   */
  BatchVector createBatchUncontractionHierarchyForVersion(const TaskGroupID task_group_id,
                                                          const size_t batch_size,
                                                          const size_t version) {
    utils::Timer::instance().start_timer("compute_root_batches", "Compute Batches For Each Root");
    const parallel::scalable_vector<HypernodeID>& roots = _contraction_tree.roots_of_version(version);
    parallel::scalable_vector<BatchVector> root_batches(roots.size(), BatchVector());
    tbb::parallel_for(0UL, roots.size(), [&](const size_t i) {
      // We create the local uncontraction order for a root of the contraction
      // tree by a simple BFS. Each BFS level induces a new batch. The representative
      // are the vertices contained in that level and the contraction partners are their childs.
      parallel::scalable_queue<HypernodeID> q;
      parallel::scalable_queue<HypernodeID> next_q;
      q.push(roots[i]);

      root_batches[i].emplace_back();
      while ( !q.empty() ) {
        const HypernodeID u = q.front();
        q.pop();

        _contraction_tree.doForEachChildOfVersion(u, version, [&](const HypernodeID v) {
          root_batches[i].back().push_back(Memento { u, v });
          next_q.push(v);
        });

        if ( q.empty() ) {
          std::swap(q, next_q);
          // Sort mementos in increasing order of their subtree size such that heavier
          // vertices are pushed earlier into the global batch vector which leads to that the
          // maximum node weight of the hypergraph degrades faster.
          std::sort(root_batches[i].back().begin(), root_batches[i].back().end(),
            [&](const Memento& lhs, const Memento& rhs) {
            return _contraction_tree.subtreeSize(lhs.v) < _contraction_tree.subtreeSize(rhs.v);
          });
          root_batches[i].emplace_back();
        }
      }

      // Remove empty batches at the end of the batch vector
      while ( root_batches[i].back().empty() ) {
        root_batches[i].pop_back();
      }
      std::reverse(root_batches[i].begin(), root_batches[i].end());
    });
    utils::Timer::instance().stop_timer("compute_root_batches");


    utils::Timer::instance().start_timer("merge_root_batches", "Merge Root Batches");
    // Push root indices into a global task queue
    tbb::concurrent_queue<size_t> root_index_queue;
    tbb::parallel_for(0UL, roots.size(), [&](const size_t i) {
      root_index_queue.push(i);
    });

    ConcurrentBatchVector tmp_batch_vec(batch_size);
    parallel::scalable_vector<uint8_t> active_tasks(std::thread::hardware_concurrency(), false);
    TBBNumaArena::instance().execute_task_on_each_thread(task_group_id,
      [&](const int, const int, const int) {
      int slot_idx = tbb::this_task_arena::current_thread_index();
      active_tasks[slot_idx] = true;
      // Indices of the roots this thread is repsonsible for
      parallel::scalable_vector<size_t> root_indices;
      // Last active batch index in which we merged a vertex of the corresponding root
      parallel::scalable_vector<size_t> last_batch_indices;

      // Try to pop roots from the global task queue. On success the current
      // thread is responsible for merging the local uncontraction order of
      // that root into the global batch vector
      size_t idx = 0;
      while ( root_index_queue.try_pop(idx) ) {
        root_indices.push_back(idx);
        last_batch_indices.push_back(0);
      }

      // Number of roots with unmerged batches
      size_t remaining_roots = root_indices.size();
      // Last index of the batch we merged into
      size_t last_merge_batch_idx = 0;
      while ( remaining_roots > 0 ) {
        bool still_active = false;

        // We visit each active root and merge at most one memento of their current batch
        // into the current global active batch (round-robin fashion)
        for ( size_t i = 0; i < remaining_roots; ++i ) {
          const size_t root_idx = root_indices[i];
          parallel::scalable_vector<Batch>& batches = root_batches[root_idx];
          ASSERT(batches.size() > 0);

          // If the current batch of the corresponding root is empty, we have to wait until the index
          // of the current active global batch is greater than the last global batch index in which
          // we merged a memento of the corresponding root. It ensures that no two mementos from different
          // local batches are merged into same global batch.
          if ( batches.back().empty() && tmp_batch_vec.activeBatchIndex() > last_batch_indices[i] ) {
            ASSERT(batches.size() > 1);
            batches.pop_back();
          }

          if ( !batches.back().empty() ) {
            const Memento memento = batches.back().back();
            still_active = true;
            batches.back().pop_back();
            last_merge_batch_idx = tmp_batch_vec.push(memento);
            last_batch_indices[i] = last_merge_batch_idx;

            // If we finish processing all batches of the current root,
            // we remove it from our local task list.
            if ( batches.back().empty() && batches.size() == 1 ) {
              batches.pop_back();
              std::swap(root_indices[i], root_indices[--remaining_roots]);
              std::swap(last_batch_indices[i--], last_batch_indices[remaining_roots]);
              root_indices.pop_back();
              last_batch_indices.pop_back();
            }
          }
        }

        // If we were not able to merge at least one memento into the global batch vector,
        // we perform busy-waiting until a new batch is available and inside the busy-waiting
        // loop we check if all active threads are idle. If all threads are idle, we force the
        // creation of new batch in the global batch vector.
        if ( !still_active && remaining_roots > 0 ) {
          active_tasks[slot_idx] = false;
          while ( last_merge_batch_idx == tmp_batch_vec.activeBatchIndex() ) {
            bool all_idle = true;
            for ( uint32_t i = 0; i < std::thread::hardware_concurrency(); ++i ) {
              all_idle &= !active_tasks[i];
            }

            if ( all_idle ) {
              tmp_batch_vec.addEmptyBatch(last_merge_batch_idx);
            }
          }
          active_tasks[slot_idx] = true;
        }
      }

      active_tasks[slot_idx] = false;
    });

    BatchVector batches = std::move(tmp_batch_vec.get());
    utils::Timer::instance().stop_timer("merge_root_batches");

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
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;
  // ! Version of the hypergraph, each time we remove a single-pin and parallel nets,
  // ! we create a new version
  size_t _version;

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

  // ! Community Information and Stats
  CommunitySupport<DynamicHypergraph> _community_support;

};

} // namespace ds
} // namespace mt_kahypar