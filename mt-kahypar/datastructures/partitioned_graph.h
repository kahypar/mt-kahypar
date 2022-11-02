/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <atomic>
#include <type_traits>
#include <mutex>

#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
class PartitionedGraph {
private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  enum class EdgeLockState : uint8_t {
    UNCONTENDED = 0,
    MOVING_SMALLER_ID_NODE = 1,
    MOVING_LARGER_ID_NODE = 2,
    LOCKED = 3
  };

  // ! Multi-state lock to synchronize moves
  struct EdgeLock {
    EdgeLock() :
      state(0),
      move_target(kInvalidPartition) {
    }

    CAtomic<uint32_t> state;
    PartitionID move_target;
  };

  class ConnectivityIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PartitionID;
    using reference = PartitionID&;
    using pointer = const PartitionID*;
    using difference_type = std::ptrdiff_t;

    /*!
     * Constructs a connectivity iterator based on a pin iterator
     */
    ConnectivityIterator(PartitionID first, PartitionID second, unsigned int count) :
      _first(first),
      _second(second),
      _iteration_count(count) {
        if (_first == _second) {
          ++_iteration_count;
        }
        if (_first == kInvalidPartition) {
          ++_iteration_count;
        } else if (_second == kInvalidPartition) {
          ++_iteration_count;
          _second = _first;
        }
        _iteration_count = std::min<unsigned int>(_iteration_count, 2);
    }

    // ! Returns the current partiton id.
    PartitionID operator* () const {
      ASSERT(_iteration_count < 2);
      return _iteration_count == 0 ? _first : _second;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    ConnectivityIterator & operator++ () {
      ASSERT(_iteration_count < 2);
      ++_iteration_count;
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    ConnectivityIterator operator++ (int) {
      ConnectivityIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const ConnectivityIterator& rhs) {
      return _first != rhs._first || _second != rhs._second ||
             _iteration_count != rhs._iteration_count;
    }

    bool operator== (const ConnectivityIterator& rhs) {
      return _first == rhs._first && _second == rhs._second &&
             _iteration_count == rhs._iteration_count;
    }


   private:
    PartitionID _first = 0;
    PartitionID _second = 0;
    // state of the iterator
    unsigned int _iteration_count = 0;
  };

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_partitioned = true;
  static constexpr bool supports_connectivity_set = true;

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);
  static constexpr size_t SIZE_OF_EDGE_LOCK = sizeof(EdgeLock);

  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;

  PartitionedGraph() = default;

  explicit PartitionedGraph(const PartitionID k,
                            Hypergraph& hypergraph) :
    _is_gain_cache_initialized(false),
    _top_level_num_nodes(hypergraph.initialNumNodes()),
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(
      "Refinement", "part_ids", hypergraph.initialNumNodes(), false, false),
    _incident_weight_in_part(),
    _edge_locks(
      "Refinement", "edge_locks", hypergraph.maxUniqueID(), false, false),
    _edge_markers(Hypergraph::is_static_hypergraph ? 0 : hypergraph.maxUniqueID()) {
    _part_ids.assign(hypergraph.initialNumNodes(), CAtomic<PartitionID>(kInvalidPartition), false);
  }

  explicit PartitionedGraph(const PartitionID k,
                            Hypergraph& hypergraph,
                            parallel_tag_t) :
    _is_gain_cache_initialized(false),
    _top_level_num_nodes(hypergraph.initialNumNodes()),
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(),
    _incident_weight_in_part(),
    _edge_locks(),
    _edge_markers() {
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "part_ids", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), CAtomic<PartitionID>(kInvalidPartition));
    }, [&] {
      _edge_locks.resize(
        "Refinement", "edge_locks", static_cast<size_t>(hypergraph.maxUniqueID()));
      _edge_locks.assign(hypergraph.maxUniqueID(), EdgeLock());
    }, [&] {
      if (!Hypergraph::is_static_hypergraph) {
        _edge_markers.setSize(hypergraph.maxUniqueID());
      }
    });
  }

  PartitionedGraph(const PartitionedGraph&) = delete;
  PartitionedGraph & operator= (const PartitionedGraph &) = delete;

  PartitionedGraph(PartitionedGraph&& other) = default;
  PartitionedGraph & operator= (PartitionedGraph&& other) = default;

  ~PartitionedGraph() {
    freeInternalData();
  }

  void resetData() {
    _is_gain_cache_initialized = false;
    resetMoveState();
    tbb::parallel_invoke([&] {
    }, [&] {
      _part_ids.assign(_part_ids.size(), CAtomic<PartitionID>(kInvalidPartition));
    }, [&] {
      _incident_weight_in_part.assign(_incident_weight_in_part.size(),  CAtomic<HyperedgeWeight>(0));
    }, [&] {
      for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);
    });
  }

  // ####################### General Hypergraph Stats ######################

  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
  }

  void setHypergraph(Hypergraph& hypergraph) {
    _hg = &hypergraph;
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _hg->initialNumNodes();
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _hg->numRemovedHypernodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _hg->initialNumEdges();
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _hg->initialNumPins();
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _hg->initialTotalVertexDegree();
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _hg->totalWeight();
  }

  // ! Number of blocks this hypergraph is partitioned into
  PartitionID k() const {
    return _k;
  }


  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    _hg->doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    _hg->doParallelForAllEdges(f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return _hg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    return _hg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    return _hg->pins(e);
  }

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivityIterator> connectivitySet(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    PartitionID first = partID(edgeSource(e));
    PartitionID second = partID(edgeTarget(e));
    return IteratorRange<ConnectivityIterator>(
      ConnectivityIterator(first, second, 0),
      ConnectivityIterator(first, second, 2));
  }

  // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return _hg->nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    const PartitionID block = partID(u);
    if ( block != kInvalidPartition ) {
      ASSERT(block < _k);
      const HypernodeWeight delta = weight - _hg->nodeWeight(u);
      _part_weights[block] += delta;
    }
    _hg->setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return _hg->nodeDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return _hg->nodeIsEnabled(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u, const PartitionID to) {
    _hg->restoreDegreeZeroHypernode(u);
    setNodePart(u, to);
  }

  // ####################### Hyperedge Information #######################

  // ! Target of an edge
  HypernodeID edgeTarget(const HyperedgeID e) const {
    return _hg->edgeTarget(e);
  }

  // ! Source of an edge
  HypernodeID edgeSource(const HyperedgeID e) const {
    return _hg->edgeSource(e);
  }

  // ! Whether the edge is a single pin edge
  bool isSinglePin(const HyperedgeID e) const {
    return _hg->isSinglePin(e);
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
  }

  // ! Unique id of a hyperedge
  HyperedgeID uniqueEdgeID(const HyperedgeID e) const {
    return _hg->uniqueEdgeID(e);
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    _hg->setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return _hg->edgeSize(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return _hg->edgeIsEnabled(e);
  }

  // ####################### Uncontraction #######################

  void uncontract(const Batch& batch) {
    resetMoveState();
    // Set block ids of contraction partners
    tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
      const Memento& memento = batch[i];
      ASSERT(nodeIsEnabled(memento.u));
      ASSERT(!nodeIsEnabled(memento.v));
      const PartitionID part_id = partID(memento.u);
      ASSERT(part_id != kInvalidPartition && part_id < _k);
      setOnlyNodePart(memento.v, part_id);
    });

    _hg->uncontract(batch,
      [&](const HyperedgeID e) { return !_edge_markers.compare_and_set_to_true(uniqueEdgeID(e)); },
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID e) {
        // In this case, e was a single pin edge before uncontraction
        if ( _is_gain_cache_initialized ) {
          // the edge weight is added to u and v
          const PartitionID block = partID(u);
          const HyperedgeWeight we = edgeWeight(e);
          _incident_weight_in_part[incident_weight_index(u, block)].fetch_add(we, std::memory_order_relaxed);
          _incident_weight_in_part[incident_weight_index(v, block)].fetch_add(we, std::memory_order_relaxed);
        }
      },
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID e) {
        // In this case, u is replaced by v in e
        if ( _is_gain_cache_initialized ) {
          // the edge weight shifts from u to v
          const PartitionID targetBlock = partID(edgeTarget(e));
          const HyperedgeWeight we = edgeWeight(e);
          _incident_weight_in_part[incident_weight_index(u, targetBlock)].fetch_sub(we, std::memory_order_relaxed);
          _incident_weight_in_part[incident_weight_index(v, targetBlock)].fetch_add(we, std::memory_order_relaxed);
        }
      });
  }

  // ####################### Restore Hyperedges #######################

  void restoreLargeEdge(const HyperedgeID& he) {
    _hg->restoreLargeEdge(he);
  }

  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<typename Hypergraph::ParallelHyperedge>& hes_to_restore) {
    _edge_markers.reset();
    _hg->restoreSinglePinAndParallelNets(hes_to_restore);
  }

  // ####################### Partition Information #######################

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    return _part_ids[u].load(std::memory_order_relaxed);
  }

  void extractPartIDs(Array<CAtomic<PartitionID>>& part_ids) {
    std::swap(_part_ids, part_ids);
  }


  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u].load() == kInvalidPartition);
    moveAssertions();
    _part_ids[u].store(p, std::memory_order_relaxed);
  }

  void setNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(_part_ids[u].load() == kInvalidPartition);
    setOnlyNodePart(u, p);
    _part_weights[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      DeltaFunc&& delta_func) {
    return changeNodePartImpl<true>(u, from, to, max_weight_to, report_success, delta_func);
  }

  // overload
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func) {
    return changeNodePart(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func);
  }

  // overload for case that does not require locking
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to) {
    return changeNodePartImpl<false>(u, from, to, std::numeric_limits<HypernodeWeight>::max(),
                                     []{}, NOOP_FUNC);
  }

  // Make sure not to call phg.gainCacheUpdate(..) in delta_func for changeNodePartWithGainCacheUpdate
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         PartitionID from,
                                         PartitionID to,
                                         HypernodeWeight max_weight_to,
                                         SuccessFunc&& report_success,
                                         DeltaFunc&& delta_func) {
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    auto my_delta_func = [&](const HyperedgeID he, const HyperedgeWeight edge_weight, const HypernodeID edge_size,
            const HypernodeID pin_count_in_from_part_after, const HypernodeID pin_count_in_to_part_after) {
      delta_func(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
      gainCacheUpdate(he, edge_weight, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
    };
    return changeNodePart(u, from, to, max_weight_to, report_success, my_delta_func);
  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u, PartitionID from, PartitionID to) {
    return changeNodePartWithGainCacheUpdate(u, from, to, std::numeric_limits<HypernodeWeight>::max(), [] { },
                                             NoOpDeltaFunc());
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_weights[p].load(std::memory_order_relaxed);
  }

  // ! Returns whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    const PartitionID part_id = partID(u);
    if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD ) {
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        if ( partID(edgeTarget(he)) != part_id ) {
          return true;
        }
      }
    }
    return false;
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    const PartitionID part_id = partID(u);
    HypernodeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( partID(edgeTarget(he)) != part_id ) {
        ++num_incident_cut_hyperedges;
      }
    }
    return num_incident_cut_hyperedges;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    const PartitionID source_id = partID(edgeSource(e));
    const PartitionID target_id = partID(edgeTarget(e));
    PartitionID sum = 0;
    if (source_id != kInvalidPartition) {
      ++sum;
    }
    if (target_id != kInvalidPartition && target_id != source_id) {
      ++sum;
    }
    return sum;
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    HypernodeID count = 0;
    if (p == partID(edgeSource(e))) {
      count++;
    }
    if (!isSinglePin(e) && p == partID(edgeTarget(e))) {
      count++;
    }
    return count;
  }

  HyperedgeWeight moveFromPenalty(const HypernodeID u) const {
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _incident_weight_in_part[incident_weight_index(u, partID(u))].load(std::memory_order_relaxed);
  }

  HyperedgeWeight moveToBenefit(const HypernodeID u, PartitionID p) const {
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _incident_weight_in_part[incident_weight_index(u, p)].load(std::memory_order_relaxed);
  }

  HyperedgeWeight incidentWeightInPart(const HypernodeID u, PartitionID p) const {
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _incident_weight_in_part[incident_weight_index(u, p)].load(std::memory_order_relaxed);
  }

  void initializeGainCacheEntry(const HypernodeID u, parallel::scalable_vector<Gain>& benefit_aggregator) {
    for (HyperedgeID e : incidentEdges(u)) {
      if (!isSinglePin(e)) {
        benefit_aggregator[partID(edgeTarget(e))] += edgeWeight(e);
      }
    }

    for (PartitionID i = 0; i < _k; ++i) {
      _incident_weight_in_part[incident_weight_index(u, i)].store(
        benefit_aggregator[i], std::memory_order_relaxed);
      benefit_aggregator[i] = 0;
    }
  }

  HyperedgeWeight km1Gain(const HypernodeID u, PartitionID from, PartitionID to) const {
    unused(from);
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveToBenefit(u, to) - moveFromPenalty(u);
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, block weights must be initialized explicitly here.
  void initializePartition() {
    initializeBlockWeights();
  }

  bool isGainCacheInitialized() const {
    return _is_gain_cache_initialized;
  }

  // ! Initialize gain cache
  void initializeGainCache() {
    allocateGainTableIfNecessary();

    // assert that part has been initialized
    ASSERT(std::none_of(nodes().begin(), nodes().end(),
                            [&](HypernodeID u) { return partID(u) == kInvalidPartition || partID(u) > k(); }) );
    // assert that current gain values are zero
    ASSERT(!_is_gain_cache_initialized
           && std::none_of(_incident_weight_in_part.begin(), _incident_weight_in_part.end(),
                           [&](const auto& weight) { return weight.load() != 0; }));

    // Calculate gain in parallel over all edges. Note that because the edges
    // are grouped by source node, this is still cache-efficient.
    doParallelForAllEdges([&](const HyperedgeID e) {
      const HypernodeID node = edgeSource(e);
      if (nodeIsEnabled(node) && !isSinglePin(e)) {
        size_t index = incident_weight_index(node, partID(edgeTarget(e)));
        _incident_weight_in_part[index].fetch_add(edgeWeight(e), std::memory_order_relaxed);
      }
    });

    _is_gain_cache_initialized = true;
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    _part_ids.assign(_part_ids.size(), CAtomic<PartitionID>(kInvalidPartition), false);
    _incident_weight_in_part.assign(_incident_weight_in_part.size(),  CAtomic<HyperedgeWeight>(0), false);
    for (auto& weight : _part_weights) {
      weight.store(0, std::memory_order_relaxed);
    }
  }

  // ! Resets the edge locks. Should be called e.g. after a rollback.
  // !
  // ! More precisely, this needs to be called if (in this order):
  // ! 1. Nodes are moved via changeNodePart(), involving a delta function
  // ! 2. Nodes are reassigned (via setOnlyNodePart() or changeNodePart() without a delta function)
  // ! 3. Nodes are moved again
  // ! Then, resetMoveState() must be called between steps 2 and 3.
  void resetMoveState() {
    if (_lock_treshold > std::numeric_limits<decltype(_lock_treshold)>::max() - 3) {
      tbb::parallel_for(ID(0), _hg->maxUniqueID(), [&](const HyperedgeID id) {
        _edge_locks[id].state.store(0, std::memory_order_relaxed);
      });
      _lock_treshold = 0;
    } else {
      _lock_treshold += 3;
    }
    moveAssertions();
  }

  // ! Only for testing
  void recomputePartWeights() {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].store(0);
    }

    for (HypernodeID u : nodes()) {
      _part_weights[ partID(u) ] += nodeWeight(u);
    }
  }

  void allocateGainTableIfNecessary() {
    if (_incident_weight_in_part.size() == 0) {
      _incident_weight_in_part.resize("Refinement", "incident_weight_in_part", _top_level_num_nodes * size_t(_k), true);
    }
  }


  // ! Only for testing
  HyperedgeWeight moveFromPenaltyRecomputed(const HypernodeID u) const {
    PartitionID part_id = partID(u);
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (!isSinglePin(e) && partID(edgeTarget(e)) == part_id) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  // ! Only for testing
  HyperedgeWeight moveToBenefitRecomputed(const HypernodeID u, PartitionID p) const {
    PartitionID part_id = partID(u);
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (!isSinglePin(e) && partID(edgeTarget(e)) == p) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  void recomputeMoveFromPenalty(const HypernodeID u) {
    // Nothing to do here
  }

  // ! Only for testing
  bool checkTrackedPartitionInformation() {
    bool success = true;

    for (HyperedgeID e : edges()) {
      PartitionID expected_connectivity = 0;
      for (PartitionID i = 0; i < k(); ++i) {
        expected_connectivity += (pinCountInPart(e, i) > 0);
      }
      if ( expected_connectivity != connectivity(e) ) {
        LOG << "Connectivity of hyperedge" << e << "=>" <<
            "Expected:" << V(expected_connectivity)  << "," <<
            "Actual:" << V(connectivity(e));
        success = false;
      }
    }

    if ( _is_gain_cache_initialized ) {
      for (HypernodeID u : nodes()) {
        if ( moveFromPenalty(u) != moveFromPenaltyRecomputed(u) ) {
          LOG << "Move from benefit of hypernode" << u << "=>" <<
              "Expected:" << V(moveFromPenaltyRecomputed(u)) << ", " <<
              "Actual:" <<  V(moveFromPenalty(u));
          success = false;
        }

        for (PartitionID i = 0; i < k(); ++i) {
          if (partID(u) != i) {
            if ( moveToBenefit(u, i) != moveToBenefitRecomputed(u, i) ) {
              LOG << "Move to penalty of hypernode" << u << "in block" << i << "=>" <<
                  "Expected:" << V(moveToBenefitRecomputed(u, i)) << ", " <<
                  "Actual:" <<  V(moveToBenefit(u, i));
              success = false;
            }
          }
        }
      }
    }
    return success;
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);

    parent->addChild("Part Weights", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Incident Weight in Part", sizeof(CAtomic<HyperedgeWeight>) * _incident_weight_in_part.size());
  }

  // ####################### Extract Block #######################

  // ! Extracts a block of a partition as separate graph.
  // ! It also returns a vertex-mapping from the original graph to the sub-graph.
  std::pair<Hypergraph, parallel::scalable_vector<HypernodeID> > extract(
    PartitionID block,
    bool /*cut_net_splitting*/,
    bool stable_construction_of_incident_edges
  ) {
    ASSERT(block != kInvalidPartition && block < _k);

    // Compactify vertex ids
    parallel::scalable_vector<HypernodeID> node_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    parallel::scalable_vector<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    HypernodeID num_nodes = 0;
    HypernodeID num_edges = 0;
    tbb::parallel_invoke([&] {
      for (const HypernodeID& node : nodes()) {
        if (partID(node) == block) {
          node_mapping[node] = num_nodes++;
        }
      }
    }, [&] {
      for (const HyperedgeID& edge : edges()) {
        const HypernodeID source = edgeSource(edge);
        const HypernodeID target = edgeTarget(edge);
        if (partID(source) == block && partID(target) == block && source < target) {
          he_mapping[edge] = num_edges++;
        }
      }
    });

    // Extract plain hypergraph data for corresponding block
    using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
    EdgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> edge_weight;
    parallel::scalable_vector<HypernodeWeight> node_weight;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_edges);
      edge_weight.resize(num_edges);
      doParallelForAllEdges([&](const HyperedgeID edge) {
        const HypernodeID source = edgeSource(edge);
        const HypernodeID target = edgeTarget(edge);
        if (partID(source) == block && partID(target) == block && source < target) {
          ASSERT(he_mapping[edge] < num_edges);
          edge_weight[he_mapping[edge]] = edgeWeight(edge);
          for (const HypernodeID& pin : pins(edge)) {
            edge_vector[he_mapping[edge]] = {node_mapping[source], node_mapping[target]};
          }
        }
      });
    }, [&] {
      node_weight.resize(num_nodes);
      doParallelForAllNodes([&](const HypernodeID node) {
        if (partID(node) == block) {
          node_weight[node_mapping[node]] = nodeWeight(node);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_graph = HypergraphFactory::construct_from_graph_edges(
               num_nodes, num_edges, edge_vector, edge_weight.data(), node_weight.data(),
               stable_construction_of_incident_edges);

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& node) {
      if (partID(node) == block) {
        const HypernodeID extracted_node = node_mapping[node];
        extracted_graph.setCommunityID(extracted_node, _hg->communityID(node));
      }
    });
    return std::make_pair(std::move(extracted_graph), std::move(node_mapping));
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      parallel::parallel_free(_part_ids, _incident_weight_in_part, _edge_locks);
    }
    _k = 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID he, const HyperedgeWeight we,
                       const PartitionID from, const HypernodeID /*pin_count_in_from_part_after*/,
                       const PartitionID to, const HypernodeID /*pin_count_in_to_part_after*/) {
    const HypernodeID target = edgeTarget(he);
    const size_t index_in_from_part = incident_weight_index(target, from);
    _incident_weight_in_part[index_in_from_part].fetch_sub(we, std::memory_order_relaxed);
    const size_t index_in_to_part = incident_weight_index(target, to);
    _incident_weight_in_part[index_in_to_part].fetch_add(we, std::memory_order_relaxed);
  }

 private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t incident_weight_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * _k  + p;
  }

  template<bool HandleLocks, typename SuccessFunc, typename DeltaFunc>
  bool changeNodePartImpl(const HypernodeID u,
                          PartitionID from,
                          PartitionID to,
                          HypernodeWeight max_weight_to,
                          SuccessFunc&& report_success,
                          DeltaFunc&& delta_func) {
    ASSERT(partID(u) == from);
    ASSERT(from != to);
    const HypernodeWeight weight = nodeWeight(u);
    const HypernodeWeight to_weight_after = _part_weights[to].add_fetch(weight, std::memory_order_relaxed);
    if (to_weight_after <= max_weight_to) {
      _part_weights[from].fetch_sub(weight, std::memory_order_relaxed);
      report_success();
      if (HandleLocks) {
        DBG << "<<< Start changing node part: " << V(u) << " - " << V(from) << " - " << V(to);
        parallel::scalable_vector<std::pair<HyperedgeID, PartitionID>> locks_to_restore;
        for (const HyperedgeID edge : incidentEdges(u)) {
          if (!isSinglePin(edge)) {
            const PartitionID target_part = targetPartWithLockSynchronization(u, to, edge, locks_to_restore);
            const HypernodeID pin_count_in_from_part_after = target_part == from ? 1 : 0;
            const HypernodeID pin_count_in_to_part_after = target_part == to ? 2 : 1;
            delta_func(edge, edgeWeight(edge), edgeSize(edge), pin_count_in_from_part_after, pin_count_in_to_part_after);
          }
        }
        _part_ids[u].store(to, std::memory_order_relaxed);
        DBG << "Done changing node part: " << V(u) << " >>>";
        restoreLockInformation(u, std::move(locks_to_restore));
      } else {
        moveAssertions();
        _part_ids[u].store(to, std::memory_order_relaxed);
      }
      return true;
    } else {
      _part_weights[to].fetch_sub(weight, std::memory_order_relaxed);
      return false;
    }
  }

  // ! Determines partition id of v, synchronized by the edge lock of {u, v}.
  // ! Information that needs to be restored is pushed to locks_to_restore.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  PartitionID targetPartWithLockSynchronization(const HypernodeID u,
                                                const PartitionID to,
                                                const HyperedgeID edge,
                                                parallel::scalable_vector<std::pair<HyperedgeID, PartitionID>>& locks_to_restore) {
    const HypernodeID v = edgeTarget(edge);
    ASSERT(u != v);
    const bool is_smaller_id = u < v;
    EdgeLock& e_lock = _edge_locks[_hg->uniqueEdgeID(edge)];
    const auto [state, part_id] = lock(e_lock, is_smaller_id);
    PartitionID target_part;
    if (state == EdgeLockState::UNCONTENDED) {
      target_part = partID(v);
      DBG << "EdgeLockState::UNCONTENDED: " << V(u) << " - " << V(v);
    } else if (state == moveState(is_smaller_id)) {
      // this state means u was already moved previously
      ASSERT(part_id == partID(u), "Lock in invalid state: " << V(part_id) << " - "
             << V(partID(u)) << " - " << V(to) << " - was resetMoveState() called properly?");
      target_part = partID(v);
      DBG << "EdgeLock::moveState(is_smaller_id): " << V(u) << " - " << V(v);
    } else if (state == moveState(!is_smaller_id)) {
      // this state means v might be moved concurrently
      target_part = part_id;
      DBG << "EdgeLock::moveState(!is_smaller_id): " << V(u) << " - " << V(v);

      if (part_id != partID(v)) {
        // v is moved concurrently, and this information must be restored afterwards
        locks_to_restore.push_back({edge, part_id});
        DBG << "Lock needs to be restored: " << V(u) << " - " << V(v);
      }
    } else {
      ALWAYS_ASSERT(false, "unreachable");
    }
    unlockMoving(e_lock, is_smaller_id, to);
    return target_part;
  }

  void restoreLockInformation(const HypernodeID u, parallel::scalable_vector<std::pair<HyperedgeID, PartitionID>>&& locks_to_restore) {
      for (const auto& [edge, part_id] : locks_to_restore) {
        const HypernodeID v = edgeTarget(edge);
        const bool is_smaller_id = u < v;
        EdgeLock& e_lock = _edge_locks[_hg->uniqueEdgeID(edge)];
        ASSERT(!isUncontended(e_lock));
        const auto [success, lock_part_id] = tryLock(e_lock, moveState(is_smaller_id));
        if (success) {
          ASSERT(lock_part_id == partID(u));
          if (partID(v) != part_id) {
            DBG << "Reset to moving: " << V(u) << " - " << V(v);
            // reset lock to moving state for v
            unlockMoving(e_lock, !is_smaller_id, part_id);
          } else {
            DBG << "Reset to uncontended: " << V(u) << " - " << V(v);
            // reset lock to uncontended state
            unlock(e_lock);
          }
        }
      }
  }

  void initializeBlockWeights() {
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
      [&](tbb::blocked_range<HypernodeID>& r) {
        // this is not enumerable_thread_specific because of the static partitioner
        parallel::scalable_vector<HypernodeWeight> part_weight_deltas(_k, 0);
        for (HypernodeID node = r.begin(); node < r.end(); ++node) {
          if (nodeIsEnabled(node)) {
            part_weight_deltas[partID(node)] += nodeWeight(node);
          }
        }
        for (PartitionID p = 0; p < _k; ++p) {
          _part_weights[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
        }
      },
      tbb::static_partitioner()
    );
  }

  void moveAssertions() {
    HEAVY_REFINEMENT_ASSERT(
      [&]{
        for (size_t id = 0; id < _hg->maxUniqueID(); ++id) {
          if (!isUncontended(_edge_locks[id])) {
            return false;
          }
        }
        return true;
      }(),
      "Some lock is in invalid state - was resetMoveState() called properly?"
    );
  }

  // ####################### Edge Locks #######################

  // ! Returns whether the lock is in UNCONTENDED state.
  bool isUncontended(const EdgeLock& e_lock) const {
    return lockState(e_lock.state.load(std::memory_order_relaxed)) == EdgeLockState::UNCONTENDED;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  EdgeLockState lockState(uint32_t value) const {
    if (value <= _lock_treshold) {
      return EdgeLockState::UNCONTENDED;
    }
    return static_cast<EdgeLockState>(value - _lock_treshold);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  uint32_t lockValue(EdgeLockState state) const {
    return _lock_treshold + static_cast<uint32_t>(state);
  }

  // ! Returns whether locking was successful and the partition id.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<bool, PartitionID> tryLock(EdgeLock& e_lock, EdgeLockState expected) {
    const uint32_t locked = lockValue(EdgeLockState::LOCKED);
    uint32_t expected_val = lockValue(expected);
    if (e_lock.state.compare_exchange_weak(expected_val, locked, std::memory_order_acquire)) {
      return {true, e_lock.move_target};
    } else {
      return {false, kInvalidPartition};
    }
  }

  // ! Returns the previous state and the partition id.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<EdgeLockState, PartitionID> lock(EdgeLock& e_lock, bool smaller_id) {
    const uint32_t locked = lockValue(EdgeLockState::LOCKED);
    uint32_t expected_val = lockValue(moveState(!smaller_id));
    while (!e_lock.state.compare_exchange_weak(expected_val, locked, std::memory_order_acquire)) {
      ASSERT(expected_val < _lock_treshold + 4);
      if (expected_val == locked) {
        expected_val = lockValue(moveState(!smaller_id));
      }
    }
    return {lockState(expected_val), e_lock.move_target};
  }

  // ! Unlocks and sets the state to uncontended.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void unlock(EdgeLock& e_lock) {
    ASSERT(e_lock.state.load() == lockValue(EdgeLockState::LOCKED));
    e_lock.move_target = kInvalidPartition;
    e_lock.state.store(lockValue(EdgeLockState::UNCONTENDED), std::memory_order_release);
  }

  // ! Unlocks, sets the state to moving and sets the partition id.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void unlockMoving(EdgeLock& e_lock, bool smaller_id, PartitionID target) {
    ASSERT(e_lock.state.load() == lockValue(EdgeLockState::LOCKED));
    e_lock.move_target = target;
    e_lock.state.store(lockValue(moveState(smaller_id)), std::memory_order_release);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  static EdgeLockState moveState(bool smaller_id) {
    return smaller_id ? EdgeLockState::MOVING_SMALLER_ID_NODE : EdgeLockState::MOVING_LARGER_ID_NODE;
  }

  // ! Indicate whether gain cache is initialized
  bool _is_gain_cache_initialized;

  size_t _top_level_num_nodes = 0;

  // ! Number of blocks
  PartitionID _k = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  parallel::scalable_vector< CAtomic<HypernodeWeight> > _part_weights;

  // ! Current block IDs of the vertices
  Array< CAtomic<PartitionID> > _part_ids;

  // ! For each node and block, the sum of incident edge weights where the target is in that part
  Array< CAtomic<HyperedgeWeight> > _incident_weight_in_part;

  // ! For each edge we use an atomic lock to synchronize moves
  Array< EdgeLock > _edge_locks;

  // ! We need to synchronize uncontractions via atomic markers
  ThreadSafeFastResetFlagArray<uint8_t> _edge_markers;

  // ! Fast reset threshold for edge locks
  uint32_t _lock_treshold;
};

} // namespace ds
} // namespace mt_kahypar
