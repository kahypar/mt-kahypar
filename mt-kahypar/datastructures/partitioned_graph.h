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

template <typename Hypergraph = Mandatory>
class PartitionedGraph {
private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  // Factory
  using HypergraphFactory = typename Hypergraph::Factory;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

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

  struct EdgeMove {
    EdgeMove() :
      lock(),
      u(kInvalidHypernode),
      to(kInvalidPartition) { }

    SpinLock lock;
    HypernodeID u;
    PartitionID to;
  };

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_graph = Hypergraph::is_graph;
  static constexpr bool is_partitioned = true;
  static constexpr bool supports_connectivity_set = true;
  static constexpr mt_kahypar_partition_type_t TYPE = PartitionedGraphType<Hypergraph>::TYPE;

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);
  static constexpr size_t SIZE_OF_EDGE_LOCK = sizeof(EdgeMove);

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
    _edge_sync(
      "Refinement", "edge_sync", hypergraph.maxUniqueID(), false, false),
    _edge_markers(Hypergraph::is_static_hypergraph ? 0 : hypergraph.maxUniqueID()) {
    _part_ids.assign(hypergraph.initialNumNodes(), CAtomic<PartitionID>(kInvalidPartition), false);
    _edge_sync.assign(hypergraph.maxUniqueID(), EdgeMove(), false);
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
    _edge_sync(),
    _edge_markers() {
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "part_ids", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), CAtomic<PartitionID>(kInvalidPartition));
    }, [&] {
      _edge_sync.resize(
        "Refinement", "edge_sync", static_cast<size_t>(hypergraph.maxUniqueID()));
      _edge_sync.assign(hypergraph.maxUniqueID(), EdgeMove());
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
    tbb::parallel_invoke([&] {
    }, [&] {
      _part_ids.assign(_part_ids.size(), CAtomic<PartitionID>(kInvalidPartition));
    }, [&] {
      if ( _is_gain_cache_initialized ) {
        _incident_weight_in_part.assign(_incident_weight_in_part.size(),  CAtomic<HyperedgeWeight>(0));
      }
    }, [&] {
      for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);
    }, [&] {
      _edge_sync.assign(_hg->maxUniqueID(), EdgeMove());
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
    // Set block ids of contraction partners
    tbb::parallel_for(UL(0), batch.size(), [&](const size_t i) {
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
    // If we pass the input hypergraph to initial partitioning, then initial partitioning
    // will pass an part ID vector of size |V'|, where V' are the number of nodes of
    // smallest hypergraph, while the _part_ids vector of the input hypergraph is initialized
    // with the original number of nodes. This can cause segmentation fault when we simply swap them
    // during main uncoarsening.
    if ( _part_ids.size() == part_ids.size() ) {
      std::swap(_part_ids, part_ids);
    } else {
      ASSERT(part_ids.size() <= _part_ids.size());
      tbb::parallel_for(UL(0), part_ids.size(), [&](const size_t i) {
        part_ids[i].store(_part_ids[i], std::memory_order_relaxed);
      });
    }
  }


  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u].load() == kInvalidPartition);
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
    return changeNodePartImpl(u, from, to, max_weight_to, report_success, delta_func);
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
    return changeNodePartImpl(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, NOOP_FUNC);
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
    return changeNodePartWithGainCacheUpdate(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), [] { }, NoOpDeltaFunc());
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
    _edge_sync.assign(_hg->maxUniqueID(), EdgeMove(), false);
    for (auto& weight : _part_weights) {
      weight.store(0, std::memory_order_relaxed);
    }
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

  void recomputeMoveFromPenalty(const HypernodeID) {
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
    parent->addChild("Part Weights", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Incident Block Weights", sizeof(CAtomic<HyperedgeWeight>) * _incident_weight_in_part.size());
    parent->addChild("Edge Synchronization", sizeof(CAtomic<PartitionID>) * _edge_sync.size());
    parent->addChild("Edge Markers", sizeof(uint8_t) * _edge_markers.size());
  }

  // ####################### Extract Block #######################

  std::pair<Hypergraph, vec<HypernodeID> > extract(
          const PartitionID block,
          bool cut_net_splitting,
          bool stable_construction_of_incident_edges) {
    ASSERT(block != kInvalidPartition && block < _k);
    vec<HypernodeID> node_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    Hypergraph block_graph = extract(block, node_mapping,
      cut_net_splitting, stable_construction_of_incident_edges);
    return std::make_pair(std::move(block_graph), std::move(node_mapping));
  }

  // ! Extracts a block of a partition as separate graph.
  // ! It also returns a vertex-mapping from the original graph to the sub-graph.
  Hypergraph extract(const PartitionID block,
                     vec<HypernodeID>& node_mapping,
                     bool /*cut_net_splitting*/,
                     bool stable_construction_of_incident_edges) {
    ASSERT(block != kInvalidPartition && block < _k);
    ASSERT(_hg->initialNumNodes() == static_cast<HypernodeID>(node_mapping.size()));

    // Compactify vertex ids
    vec<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
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
    using EdgeVector = vec<std::pair<HypernodeID, HypernodeID>>;
    EdgeVector edge_vector;
    vec<HyperedgeWeight> edge_weight;
    vec<HypernodeWeight> node_weight;
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
    return extracted_graph;
  }

  std::pair<vec<Hypergraph>, vec<HypernodeID>> extractAllBlocks(const PartitionID k,
                                                                const bool /*cut_net_splitting*/,
                                                                const bool stable_construction_of_incident_edges) {
    ASSERT(k <= _k);

    // Compactify node and edge ids
    vec<HypernodeID> hn_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    vec<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    vec<parallel::AtomicWrapper<HypernodeID>> nodes_cnt(
      k, parallel::AtomicWrapper<HypernodeID>(0));
    vec<parallel::AtomicWrapper<HyperedgeID>> edges_cnt(
      k, parallel::AtomicWrapper<HyperedgeID>(0));
    if ( stable_construction_of_incident_edges ) {
      // Stable construction for deterministic behavior requires
      // to determine node and edge IDs sequentially
      tbb::parallel_invoke([&] {
        for ( const HypernodeID& hn : nodes() ) {
          const PartitionID block = partID(hn);
          if ( block < k ) {
            hn_mapping[hn] = nodes_cnt[block]++;
          }
        }
      }, [&] {
        for ( const HyperedgeID& he : edges() ) {
          const HypernodeID source = edgeSource(he);
          const HypernodeID target = edgeTarget(he);
          const PartitionID sourceBlock = partID(source);
          const PartitionID targetBlock = partID(target);
          if (source < target && sourceBlock == targetBlock && sourceBlock < k) {
            he_mapping[he] = edges_cnt[sourceBlock]++;
          }
        }
      });
    } else {
      tbb::parallel_invoke([&] {
        doParallelForAllNodes([&](const HypernodeID& hn) {
          const PartitionID block = partID(hn);
          if ( block < k ) {
            hn_mapping[hn] = nodes_cnt[block]++;
          }
        });
      }, [&] {
        doParallelForAllEdges([&](const HyperedgeID& he) {
          const HypernodeID source = edgeSource(he);
          const HypernodeID target = edgeTarget(he);
          const PartitionID sourceBlock = partID(source);
          const PartitionID targetBlock = partID(target);
          if (source < target && sourceBlock == targetBlock && sourceBlock < k) {
            he_mapping[he] = edges_cnt[sourceBlock]++;
          }
        });
      });
    }

    using EdgeVector = vec<std::pair<HypernodeID, HypernodeID>>;
    vec<EdgeVector> edge_vector(k);
    vec<vec<HyperedgeWeight>> edge_weight(k);
    vec<vec<HypernodeWeight>> node_weight(k);
    // Allocate auxilliary graph data structures
    tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
      const HypernodeID num_nodes = nodes_cnt[p];
      const HyperedgeID num_edges = edges_cnt[p];
      tbb::parallel_invoke([&] {
        edge_vector[p].resize(num_edges);
      }, [&] {
        edge_weight[p].resize(num_edges);
      }, [&] {
        node_weight[p].resize(num_nodes);
      });
    });

    // Write blocks to auxilliary graph data structure
    tbb::parallel_invoke([&] {
      doParallelForAllEdges([&](const HyperedgeID& he) {
        const HyperedgeID mapped_he = he_mapping[he];
        const HypernodeID source = edgeSource(he);
        const HypernodeID target = edgeTarget(he);
        const PartitionID sourceBlock = partID(source);
        const PartitionID targetBlock = partID(target);
        if (source < target && sourceBlock == targetBlock && sourceBlock < k) {
          ASSERT(UL(mapped_he) < edge_weight[sourceBlock].size());
          edge_weight[sourceBlock][mapped_he] = edgeWeight(he);
          edge_vector[sourceBlock][mapped_he] =
            { hn_mapping[edgeSource(he)], hn_mapping[edgeTarget(he)] };
        }
      });
    }, [&] {
      doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID block = partID(hn);
        const HypernodeID mapped_hn = hn_mapping[hn];
        if ( block < k ) {
          ASSERT(UL(mapped_hn) < node_weight[block].size());
          node_weight[block][mapped_hn] = nodeWeight(hn);
        }
      });
    });

    // Construct graph of each block
    vec<Hypergraph> extracted_graphs(k);
    tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
      const HypernodeID num_nodes = nodes_cnt[p];
      const HyperedgeID num_edges = edges_cnt[p];
      extracted_graphs[p] = HypergraphFactory::construct_from_graph_edges(
        num_nodes, num_edges, edge_vector[p], edge_weight[p].data(), node_weight[p].data(),
        stable_construction_of_incident_edges);
    });

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID block = partID(hn);
      if ( block < k ) {
        extracted_graphs[block].setCommunityID(hn_mapping[hn], _hg->communityID(hn));
      }
    });

    parallel::parallel_free(edge_vector);
    parallel::parallel_free(edge_weight, node_weight);

    return std::make_pair(std::move(extracted_graphs), std::move(hn_mapping));
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      parallel::parallel_free(_part_ids, _incident_weight_in_part, _edge_sync);
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

  template<typename SuccessFunc, typename DeltaFunc>
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
      DBG << "<<< Start changing node part: " << V(u) << " - " << V(from) << " - " << V(to);
      for (const HyperedgeID edge : incidentEdges(u)) {
        if (!isSinglePin(edge)) {
          const PartitionID block_of_target_node = synchronizeMoveOnEdge(edge, u, to);
          const HypernodeID pin_count_in_from_part_after = block_of_target_node == from ? 1 : 0;
          const HypernodeID pin_count_in_to_part_after = block_of_target_node == to ? 2 : 1;
          delta_func(edge, edgeWeight(edge), edgeSize(edge), pin_count_in_from_part_after, pin_count_in_to_part_after);
        }
      }
      _part_ids[u].store(to, std::memory_order_relaxed);
      DBG << "Done changing node part: " << V(u) << " >>>";
      return true;
    } else {
      _part_weights[to].fetch_sub(weight, std::memory_order_relaxed);
      return false;
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

  // ####################### Edge Locks #######################

  // This function synchronizes a move on an edge and returns the block ID
  // of the target node of the corresponding edge. The function assumes that
  // node u is moved to the block 'to'.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  PartitionID synchronizeMoveOnEdge(const HyperedgeID edge,
                                    const HypernodeID u,
                                    const PartitionID to) {
    const HypernodeID v = edgeTarget(edge);
    PartitionID block_of_v = partID(v);
    EdgeMove& edge_move = _edge_sync[uniqueEdgeID(edge)];
    edge_move.lock.lock();
    if ( edge_move.u == v ) {
      ASSERT(edge_move.to < _k && edge_move.to != kInvalidPartition);
      block_of_v = edge_move.to;
    }
    edge_move.u = u;
    edge_move.to = to;
    edge_move.lock.unlock();
    return block_of_v;
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
  Array< EdgeMove > _edge_sync;

  // ! We need to synchronize uncontractions via atomic markers
  ThreadSafeFastResetFlagArray<uint8_t> _edge_markers;
};

} // namespace ds
} // namespace mt_kahypar
