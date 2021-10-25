/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <atomic>
#include <type_traits>
#include <mutex>

#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/datastructures/async/async_common.h"
#include "mt-kahypar/datastructures/gain_cache.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"


namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory,
          typename GainCacheT = Mandatory>
class PartitionedHypergraph {
private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  struct NoOpGainCacheUpdateFunc {
      template<typename PinIteratorT>
      void operator() (const HyperedgeID, const HyperedgeWeight, IteratorRange<PinIteratorT>, const PartitionID, const HypernodeID, const PartitionID, const HypernodeID) {}
  };

  // REVIEW NOTE: Can't we use a lambda in changeNodePart. And write a second function that calls the first with a lambda that does nothing.
  // Then we could guarantee inlining
  // This would also reduce the code/documentation copy-pasta for with or without gain updates

  static constexpr bool enable_heavy_assert = true;

  using PHG = PartitionedHypergraph<Hypergraph, HypergraphFactory, GainCacheT>;

 public:

  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_partitioned = true;
  static constexpr bool supports_connectivity_set = true;

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;

  using GainCache = GainCacheT;

  PartitionedHypergraph() = default;

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
          _is_gain_cache_initialized(false),
          _k(k),
          _hg(&hypergraph),
          _part_weights(k, CAtomic<HypernodeWeight>(0)),
          _part_ids(
        "Refinement", "part_ids", hypergraph.initialNumNodes(), false, false),
          _pins_in_part(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize(), false),
          _connectivity_set(hypergraph.initialNumEdges(), k, false),
          _gain_cache(hypergraph.initialNumNodes(), k, &_part_ids, &_pins_in_part, &_connectivity_set),
          _pin_count_update_ownership(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true, false),
          _direct_conn_set_snapshots(_k, kInvalidPartition),
          _direct_parts_with_one_pin_snapshots(_k, kInvalidPartition),
          _pin_count_in_part_bitcopy_snapshots(),
          _connectivity_set_bitcopy_snapshots(),
          _pins_snapshots(_hg->maxEdgeSize(),invalidNode),
          _num_stable_pins_seen(0),
          _num_volatile_pins_seen(0){
    _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition, false);
  }

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph,
                                 parallel_tag_t) :
          _is_gain_cache_initialized(false),
          _k(k),
          _hg(&hypergraph),
          _part_weights(k, CAtomic<HypernodeWeight>(0)),
          _part_ids(),
          _pins_in_part(),
          _connectivity_set(0, 0),
          _gain_cache(&_part_ids, &_pins_in_part, &_connectivity_set, parallel_tag_t()),
          _pin_count_update_ownership(),
          _direct_conn_set_snapshots(_k, kInvalidPartition),
          _direct_parts_with_one_pin_snapshots(_k, kInvalidPartition),
          _pin_count_in_part_bitcopy_snapshots(),
          _connectivity_set_bitcopy_snapshots(),
          _pins_snapshots(_hg->maxEdgeSize(), invalidNode),
          _num_stable_pins_seen(0),
          _num_volatile_pins_seen(0)  {
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "vertex_part_info", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition);
    }, [&] {
      _pins_in_part.initialize(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize());
    }, [&] {
      _connectivity_set = ConnectivitySets(hypergraph.initialNumEdges(), k);
    }, [&] {
      _gain_cache.resize("Refinement", "gain_cache", hypergraph.initialNumNodes(), k);
    }, [&] {
      _pin_count_update_ownership.resize(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true);
    });
  }

  // REVIEW NOTE why do we delete copy assignment/construction? wouldn't it be useful to make a copy, e.g. for initial partitioning
  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other)  noexcept :
          _is_gain_cache_initialized(other._is_gain_cache_initialized),
          _k(other._k),
          _hg(std::move(other._hg)),
          _part_weights(std::move(other._part_weights)),
          _part_ids(std::move(other._part_ids)),
          _pins_in_part(std::move(other._pins_in_part)),
          _connectivity_set(std::move(other._connectivity_set)),
          _gain_cache(std::move(other._gain_cache)),
          _pin_count_update_ownership(std::move(other._pin_count_update_ownership)),
          _direct_conn_set_snapshots(std::move(other._direct_conn_set_snapshots)),
          _direct_parts_with_one_pin_snapshots(std::move(other._direct_parts_with_one_pin_snapshots)),
          _pin_count_in_part_bitcopy_snapshots(std::move(other._pin_count_in_part_bitcopy_snapshots)),
          _connectivity_set_bitcopy_snapshots(std::move(other._connectivity_set_bitcopy_snapshots)),
          _pins_snapshots(std::move(other._pins_snapshots)),
          _num_stable_pins_seen(std::move(other._num_stable_pins_seen)),
          _num_volatile_pins_seen(std::move(other._num_volatile_pins_seen))  {

      // Reset query functions for GainCache (so references point to functions in new PHG)
      _gain_cache.assignQueryObjects(&_part_ids, &_pins_in_part, &_connectivity_set);
  }

  PartitionedHypergraph & operator= (PartitionedHypergraph&& other)  noexcept {
      _is_gain_cache_initialized = other._is_gain_cache_initialized;
      _k = other._k;
      _hg = std::move(other._hg);
      _part_weights = std::move(other._part_weights);
      _part_ids = std::move(other._part_ids);
      _pins_in_part = std::move(other._pins_in_part);
      _connectivity_set = std::move(other._connectivity_set);
      _gain_cache = std::move(other._gain_cache);
      _pin_count_update_ownership = std::move(other._pin_count_update_ownership);
    _direct_conn_set_snapshots = std::move(other._direct_conn_set_snapshots);
    _direct_parts_with_one_pin_snapshots = std::move(other._direct_parts_with_one_pin_snapshots);
    _pin_count_in_part_bitcopy_snapshots = std::move(other._pin_count_in_part_bitcopy_snapshots);
    _connectivity_set_bitcopy_snapshots = std::move(other._connectivity_set_bitcopy_snapshots);
      _pins_snapshots = std::move(other._pins_snapshots);
      _num_stable_pins_seen = std::move(other._num_stable_pins_seen);
      _num_volatile_pins_seen = std::move(other._num_volatile_pins_seen);

      _gain_cache.assignQueryObjects(&_part_ids, &_pins_in_part, &_connectivity_set);

      return *this;
  }

  ~PartitionedHypergraph() {
    freeInternalData();
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

  // ! Version of the underlying hypergraph
  size_t version() const {
      return _hg->version();
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    _hg->doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllEdges(f);
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
  IteratorRange<ConnectivitySets::Iterator> connectivitySet(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    return _connectivity_set.connectivitySet(e);
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

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    _hg->enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    _hg->disableHypernode(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u, const PartitionID to) {
    _hg->restoreDegreeZeroHypernode(u);
    setNodePart(u, to);
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
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

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    _hg->enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    _hg->disableHyperedge(e);
  }

  bool isGraphEdge(const HyperedgeID e) const {
    return _hg->isGraphEdge(e);
  }

  HyperedgeID graphEdgeID(const HyperedgeID e) const {
    return _hg->graphEdgeID(e);
  }

  HyperedgeID nonGraphEdgeID(const HyperedgeID e) const {
    return _hg->nonGraphEdgeID(e);
  }

  HypernodeID graphEdgeHead(const HyperedgeID e, const HypernodeID tail) const {
    return _hg->graphEdgeHead(e, tail);
  }

  // ! Sets snapshot edge size threshold. Not to be used during uncoarsening!
  void setSnapshotEdgeSizeThreshold(const size_t threshold) {
      _hg->setSnapshotEdgeSizeThreshold(threshold);
  }

  // ####################### Uncontraction #######################

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   */
  void uncontract(const Batch& batch) {
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
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        // In this case, u and v are incident to hyperedge he after uncontraction
        const PartitionID block = partID(u);
        const HypernodeID pin_count_in_part_after = incrementPinCountInPartWithoutGainUpdate(he, block);
        ASSERT(pin_count_in_part_after > 1, V(u) << V(v) << V(he));

        if ( _is_gain_cache_initialized ) {
          // If u was the only pin of hyperedge he in its block before then moving out vertex u
          // of hyperedge he does not decrease the connectivity any more after the
          // uncontraction => b(u) -= w(he)
          const HyperedgeWeight edge_weight = edgeWeight(he);
            _gain_cache.syncUpdateForUncontractCaseOne(he, edge_weight, v, block, pin_count_in_part_after, pins(he));
        }
      },
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        // In this case, u is replaced by v in hyperedge he
        // => Pin counts of hyperedge he does not change
        if ( _is_gain_cache_initialized ) {
          const HyperedgeWeight edge_weight = edgeWeight(he);
          _gain_cache.syncUpdateForUncontractCaseTwo(he, edge_weight, u, v);
        }
      });
  }

//  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void takeConnectivitySetSnapshotsFromBitcopySnapshots(const PinCountInPart::Snapshot& pcip_snapshot, const ConnectivitySets::Snapshot& cs_snapshot,
                                                        CompressedConnectivitySetSnapshot& conn_set_snapshot, CompressedConnectivitySetSnapshot& parts_with_one_pin_snapshot) {
      ASSERT(conn_set_snapshot.empty());
      ASSERT(parts_with_one_pin_snapshot.empty());
      for (auto block : cs_snapshot.snapshottedConnectivitySet()) {
        auto pcip = pcip_snapshot.pinCountInPart(block);
        if (pcip >= 1) {
          conn_set_snapshot.push_back(block);
          if (pcip == 1) {
            parts_with_one_pin_snapshot.push_back(block);
          }
        }
      }
  }

  void takeConnectivitySnapshotDirectly(const HyperedgeID he, CompressedConnectivitySetSnapshot& conn_set_snapshot,
                                        CompressedConnectivitySetSnapshot& parts_with_one_pin_snapshot) {
    ASSERT(conn_set_snapshot.empty());
    ASSERT(parts_with_one_pin_snapshot.empty());
    for (PartitionID block : _connectivity_set.connectivitySet(he)) {
      conn_set_snapshot.push_back(block);
      if (pinCountInPart(he, block) == 1) {
        parts_with_one_pin_snapshot.push_back(block);
      }
    }
  }

//  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  IteratorRange<HypernodeID*> takeVolatilePinsSnapshot(const HyperedgeID he) {
    std::vector<HypernodeID>& pins_snapshot = _pins_snapshots.local();
    auto pin_range = _hg->volatile_pins(he);
    HypernodeID* ptr_to_begin = &(*pin_range.begin());
    HypernodeID num_pins = pin_range.end() - pin_range.begin();
    ASSERT(num_pins <= pins_snapshot.size());
    std::memcpy(pins_snapshot.data(), ptr_to_begin, static_cast<size_t>(num_pins) * sizeof(HypernodeID));

    return IteratorRange(pins_snapshot.data(), pins_snapshot.data() + num_pins);
  }

//  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  IteratorRange<IncidenceIterator> takeStablePinsSnapshot(const HyperedgeID he) const {
    return _hg->stable_pins(he);
  }

//  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  IteratorRange<PinSnapshotIterator> takePinsSnapshot(const HyperedgeID he) {

    IteratorRange<IncidenceIterator> stable_pins = takeStablePinsSnapshot(he);
    IteratorRange<HypernodeID*> volatile_pins = takeVolatilePinsSnapshot(he);

    size_t num_stable = stable_pins.end() - stable_pins.begin();
    size_t num_volatile = volatile_pins.end() - volatile_pins.begin();
    _num_stable_pins_seen.fetch_add(num_stable, std::memory_order_relaxed);
    _num_volatile_pins_seen.fetch_add(num_volatile, std::memory_order_relaxed);

    return PinSnapshotIterator::stitchPinIterators(stable_pins, volatile_pins);
  }

  size_t getNumStablePinsSeen() const {
    return _num_stable_pins_seen.load(std::memory_order_relaxed);
  }

  size_t getNumVolatilePinsSeen() const {
    return _num_volatile_pins_seen.load(std::memory_order_relaxed);
  }

  // ! Debug purposes
  [[maybe_unused]] void lockHyperedgePinCountLock(const HyperedgeID he) {
      _pin_count_update_ownership[he].lock();
  }

  // ! Debug purposes
  [[maybe_unused]] void unlockHyperedgePinCountLock(const HyperedgeID he) {
      _pin_count_update_ownership[he].unlock();
  }

  void uncontract(const ContractionGroup& group,
                  const ContractionGroupID groupID,
                  bool useBitcopySnapshots = false) {

      // Set block ids of contraction partners
     for (auto& memento : group) {
          ASSERT(nodeIsEnabled(memento.u));
          ASSERT(!nodeIsEnabled(memento.v));
          const PartitionID part_id = partID(memento.u);
          ASSERT(part_id != kInvalidPartition && part_id < _k);
          setOnlyNodePart(memento.v, part_id);
     }

    CompressedConnectivitySetSnapshot& conn_set_snapshot = _direct_conn_set_snapshots.local();
    CompressedConnectivitySetSnapshot& parts_with_one_pin_snapshot = _direct_parts_with_one_pin_snapshots.local();
    std::unique_ptr<PinCountInPart::Snapshot>& pcip_bitcopy_snapshot = _pin_count_in_part_bitcopy_snapshots.local();
    if (!pcip_bitcopy_snapshot) {
      pcip_bitcopy_snapshot = std::make_unique<PinCountInPart::Snapshot>(_pins_in_part.generateEmptySnapshot());
    }
    std::unique_ptr<ConnectivitySets::Snapshot>& cs_bitcopy_snapshot = _connectivity_set_bitcopy_snapshots.local();
    if (!cs_bitcopy_snapshot) {
      cs_bitcopy_snapshot = std::make_unique<ConnectivitySets::Snapshot>(_connectivity_set.generateEmptySnapshot());
    }

     _hg->uncontract(group,
                     groupID,
                      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
                            // (This is always called while pin count update ownership for he has been locked by _hg->uncontract();
                            // release that lock here at right point)
                            // In this case, u and v are incident to hyperedge he after uncontraction

                          // Update pin count and take snapshot of pin state for gain cache update. This way the gain
                          // cache does not have to be updated within the lock.
                          const PartitionID block = partID(u);
                          const HypernodeID pin_count_in_part_after = incrementPinCountInPartWithoutGainUpdate(
                                  he, block);
                          ASSERT(pin_count_in_part_after > 1, V(u) << V(v) << V(he));
                          if (_is_gain_cache_initialized) {

                            const HyperedgeWeight edge_weight = edgeWeight(he);
                            if (_hg->totalSize(he) < _hg->snapshotEdgeSizeThreshold() || !_gain_cache.isPinCountThatTriggersUpdateForAllPinsForUncontractCaseOne(pin_count_in_part_after)) {
                                // Gain cache update within lock, no snapshots
                                _gain_cache.syncUpdateForUncontractCaseOne(he, edge_weight, v, block, pin_count_in_part_after, pins(he));
                                _pin_count_update_ownership[he].unlock();
                              } else {
                                // If e is large enough and update to all pins of e is necessary, perform gain cache update outside of lock:
                                // Snapshot pins and connectivity set in lock, gain cache update outside of lock

                                auto pins_snapshot = std::make_unique<IteratorRange<PinSnapshotIterator>>(takePinsSnapshot(he));
                                if (useBitcopySnapshots) {
                                  _connectivity_set.takeBitcopySnapshotForHyperedge(he, *cs_bitcopy_snapshot);
                                  _pins_in_part.takeBitcopySnapshotForHyperedge(he, *pcip_bitcopy_snapshot);
                                  _pin_count_update_ownership[he].unlock();
                                  _gain_cache.asyncUpdateForUncontractCaseOne(edge_weight, v, block,
                                                                              pin_count_in_part_after, *pins_snapshot,
                                                                              *cs_bitcopy_snapshot, *pcip_bitcopy_snapshot);

                                } else {
                                  takeConnectivitySnapshotDirectly(he, conn_set_snapshot, parts_with_one_pin_snapshot);
                                  _pin_count_update_ownership[he].unlock();
                                  _gain_cache.asyncUpdateForUncontractCaseOne(edge_weight, v, block,
                                                                              pin_count_in_part_after, *pins_snapshot,
                                                                              conn_set_snapshot, parts_with_one_pin_snapshot);
                                  conn_set_snapshot.clear();
                                  parts_with_one_pin_snapshot.clear();
                                }
                              }
                          } else {
                              _pin_count_update_ownership[he].unlock();
                          }
                          },
                      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
                          // (This is always called while pin count update ownership for he has been locked by _hg->uncontract();
                          // release that lock here at right point)
                          // In this case, u is replaced by v in hyperedge he
                          // => Pin counts of hyperedge he does not change
                          if (_is_gain_cache_initialized) {

//                            if (_hg->totalSize(he) < _hg->snapshotEdgeSizeThreshold()) {
                              // Gain cache update within lock, no snapshots
                              const HyperedgeWeight edge_weight = edgeWeight(he);
                              _gain_cache.syncUpdateForUncontractCaseTwo(he, edge_weight, u, v);
                              _pin_count_update_ownership[he].unlock();
//                            } else {
//
//                              // In this case only v was part of hyperedge e before and
//                              // u must be replaced by v in hyperedge e
//                              const HyperedgeWeight edge_weight = edgeWeight(he);
//
//                              if (useBitcopySnapshots) {
//                                _connectivity_set.takeBitcopySnapshotForHyperedge(he, *cs_bitcopy_snapshot);
//                                _pins_in_part.takeBitcopySnapshotForHyperedge(he, *pcip_bitcopy_snapshot);
//                                _pin_count_update_ownership[he].unlock();
//                                _gain_cache.asyncUpdateForUncontractCaseTwo(edge_weight, u, v, *cs_bitcopy_snapshot, *pcip_bitcopy_snapshot);
//
//                              } else {
//                                takeConnectivitySnapshotDirectly(he, conn_set_snapshot, parts_with_one_pin_snapshot);
//                                _pin_count_update_ownership[he].unlock();
//                                _gain_cache.asyncUpdateForUncontractCaseTwo(edge_weight, u, v, conn_set_snapshot, parts_with_one_pin_snapshot);
//                                conn_set_snapshot.clear();
//                                parts_with_one_pin_snapshot.clear();
//                              }
//                            }
                          } else {
                              _pin_count_update_ownership[he].unlock();
                          }
                      },
                      [&](const HyperedgeID he) {
                            return _pin_count_update_ownership[he].tryLock();
     });
  }


    // ####################### Restore Hyperedges #######################

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    _hg->restoreLargeEdge(he);

    // Recalculate pin count in parts
    const size_t incidence_array_start = _hg->hyperedge(he).firstEntry();
    const size_t incidence_array_end = _hg->hyperedge(he).firstInvalidEntry();
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _hg->_incidence_array[pos];
      const PartitionID block = partID(pin);
      ++ets_pin_count_in_part.local()[block];
    });

    // Aggregate local pin count for each block
    for ( PartitionID block = 0; block < _k; ++block ) {
      HypernodeID pin_count_in_part = 0;
      for ( const vec<HypernodeID>& local_pin_count : ets_pin_count_in_part ) {
        pin_count_in_part += local_pin_count[block];
      }

      if ( pin_count_in_part > 0 ) {
        _pins_in_part.setPinCountInPart(he, block, pin_count_in_part);
        _connectivity_set.add(he, block);
      }
    }
  }

  /**
   * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
   * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
   */
  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore) {
    // Restore hyperedges in hypergraph
    _hg->restoreSinglePinAndParallelNets(hes_to_restore);

    // Compute pin counts of restored hyperedges and gain cache values of vertices contained
    // single-pin hyperedges. Note, that restoring parallel hyperedges does not change any
    // value in the gain cache, since it already contributes to the gain via its representative.
    utils::Timer::instance().start_timer("update_pin_counts_and_gain_cache", "Update Pin Counts And Gain Cache");
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    tbb::parallel_for(0UL, hes_to_restore.size(), [&](const size_t i) {
      const HyperedgeID he = hes_to_restore[i].removed_hyperedge;
      const HyperedgeID representative = hes_to_restore[i].representative;
      ASSERT(edgeIsEnabled(he));
      const bool is_single_pin_he = edgeSize(he) == 1;
      if ( is_single_pin_he ) {
        // Restore single-pin net
        HypernodeID single_vertex_of_he = kInvalidHypernode;
        for ( const HypernodeID& pin : pins(he) ) {
          single_vertex_of_he = pin;
        }
        ASSERT(single_vertex_of_he != kInvalidHypernode);

        const PartitionID block_of_single_pin = partID(single_vertex_of_he);
        _connectivity_set.add(he, block_of_single_pin);
        _pins_in_part.setPinCountInPart(he, block_of_single_pin, 1);

        if ( _is_gain_cache_initialized ) {
          const HyperedgeWeight edge_weight = edgeWeight(he);
          _gain_cache.updateForRestoringSinglePinNet(edge_weight, single_vertex_of_he, block_of_single_pin);
        }
      } else {
        // Restore parallel net => pin count information given by representative
        ASSERT(edgeIsEnabled(representative));
        for ( const PartitionID& block : connectivitySet(representative) ) {
          _connectivity_set.add(he, block);
          _pins_in_part.setPinCountInPart(he, block, pinCountInPart(representative, block));
        }

        HEAVY_REFINEMENT_ASSERT([&] {
          for ( PartitionID block = 0; block < _k; ++block ) {
            if ( pinCountInPart(he, block) != pinCountInPartRecomputed(he, block) ) {
              LOG << "Pin count in part of hyperedge" << he << "in block" << block
                  << "is" << pinCountInPart(he, block) << ", but should be"
                  << pinCountInPartRecomputed(he, block);
              return false;
            }
          }
          return true;
        }());
      }
    });
    utils::Timer::instance().stop_timer("update_pin_counts_and_gain_cache");
  }

  // ####################### Partition Information #######################

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    return _part_ids[u];
  }

  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u] == kInvalidPartition);
    _part_ids[u] = p;
  }

  void setNodePart(const HypernodeID u, PartitionID p) {
    setOnlyNodePart(u, p);
    _part_weights[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
    for (HyperedgeID he : incidentEdges(u)) {
      _pin_count_update_ownership[he].lock();
      incrementPinCountInPartWithoutGainUpdate(he, p);
      _pin_count_update_ownership[he].unlock();
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.

  template<typename SuccessFunc, typename DeltaFunc, typename GainCacheUpdateFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      DeltaFunc&& delta_func,
                      GainCacheUpdateFunc&& gain_cache_update_func,
                      bool concurrent_uncontractions = false) {
      assert(nodeIsEnabled(u));
      assert(partID(u) == from);
      assert(from != to);
      ASSERT(from != kInvalidPartition);
      ASSERT(to != kInvalidPartition);
      const HypernodeWeight wu = nodeWeight(u);
      const HypernodeWeight to_weight_after = _part_weights[to].add_fetch(wu, std::memory_order_relaxed);
      const HypernodeWeight from_weight_after = _part_weights[from].fetch_sub(wu, std::memory_order_relaxed);
    if ( to_weight_after <= max_weight_to && from_weight_after > 0 ) {
      _part_ids[u] = to;
      report_success();

      if (concurrent_uncontractions) {
        std::vector<HyperedgeID> retry_edges;
        for ( const HyperedgeID he : incidentEdges(u) ) {
          bool success = asyncUpdatePinCountOfHyperedge(he, from, to, delta_func, gain_cache_update_func);
          if (!success) retry_edges.push_back(he);
        }

        size_t i = 0;
        while (!retry_edges.empty()) {
          const HyperedgeID he = retry_edges[i];
          bool success = asyncUpdatePinCountOfHyperedge(he, from, to, delta_func, gain_cache_update_func);
          if (success) {
            retry_edges[i] = retry_edges.back();
            retry_edges.pop_back();
          } else {
            ++i;
          }
          if (i >= retry_edges.size()) {
            i = 0;
          }
        }
      } else {
        for ( const HyperedgeID he : incidentEdges(u) ) {
          updatePinCountOfHyperedge(he, from, to, delta_func, gain_cache_update_func);
        }
      }

      return true;
    } else {
      _part_weights[to].fetch_sub(wu, std::memory_order_relaxed);
      _part_weights[from].fetch_add(wu, std::memory_order_relaxed);
      return false;
    }
  }

  // curry
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      DeltaFunc&& delta_func,
                      bool concurrent_uncontractions = false) {
     return changeNodePart(u, from, to, max_weight_to, report_success, delta_func, NoOpGainCacheUpdateFunc(), concurrent_uncontractions);
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC,
                      bool concurrent_uncontractions = false) {
    return changeNodePart(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func, concurrent_uncontractions);
  }


  template <typename PinIteratorT>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID, const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                         const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                         const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
    _gain_cache.updateForMove(we, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
  }

  struct GainCacheUpdateFuncWrapper {

        explicit GainCacheUpdateFuncWrapper(PartitionedHypergraph& caller) : _parent_phg(caller) {}

        PartitionedHypergraph& _parent_phg;

        template<typename PinIteratorT>
        void operator() (const HyperedgeID he, const HyperedgeWeight edge_weight, IteratorRange<PinIteratorT> pins,
                         const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                         const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
          _parent_phg.gainCacheUpdate(he, edge_weight, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
        }
  };

  // Make sure not to call phg.gainCacheUpdate(..) in delta_func for changeNodePartWithGainCacheUpdate
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         PartitionID from,
                                         PartitionID to,
                                         HypernodeWeight max_weight_to,
                                         SuccessFunc&& report_success,
                                         DeltaFunc&& delta_func,
                                         bool concurrent_uncontractions = false) {
    ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
//
//    auto gain_cache_update = [&](const HyperedgeID he, const HyperedgeWeight edge_weight, IteratorRange<std::vector<HypernodeID>::const_iterator> pins,
//                        const PartitionID from, const HypernodeID pin_count_in_from_part_after,
//                        const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
//        gainCacheUpdate(he, edge_weight, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
//    };

    return changeNodePart(u, from, to, max_weight_to, report_success, delta_func, GainCacheUpdateFuncWrapper(*this), concurrent_uncontractions);

  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u, PartitionID from, PartitionID to,
                                         bool concurrent_uncontractions = false) {
    return changeNodePartWithGainCacheUpdate(u, from, to, std::numeric_limits<HypernodeWeight>::max(), [] { },
                                             NoOpDeltaFunc(), concurrent_uncontractions);
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_weights[p].load(std::memory_order_relaxed);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD ) {
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        if ( connectivity(he) > 1 ) {
          return true;
        }
      }
      return false;
    } else {
      // TODO maybe we should allow these in label propagation? definitely not in FM
      // In case u is a high degree vertex, we omit the border node check and
      // and return false. Assumption is that it is very unlikely that such a
      // vertex can change its block.
      return false;
    }
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    HypernodeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        ++num_incident_cut_hyperedges;
      }
    }
    return num_incident_cut_hyperedges;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    return _connectivity_set.connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    return _pins_in_part.pinCountInPart(e, p);
  }

  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _gain_cache.moveFromBenefit(u);
  }

  HyperedgeWeight moveFromBenefit(const HypernodeID u, const PartitionID p) const {
      return _gain_cache.moveFromBenefit(u, p);
  }

  HyperedgeWeight moveToPenalty(const HypernodeID u, PartitionID p) const {
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _gain_cache.moveToPenalty(u, p);
  }

  void initializeGainCacheEntry(const HypernodeID u, vec<Gain>& benefit_aggregator, vec<Gain>& penalty_aggregator) {
    PartitionID pu = partID(u);
    Gain incident_edges_weight = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      HyperedgeWeight ew = edgeWeight(e);
      if (pinCountInPart(e, pu) == 1) {
        benefit_aggregator[pu] += ew;
      }
      for (PartitionID i : connectivitySet(e)) {
        penalty_aggregator[i] += ew;
      }
      incident_edges_weight += ew;
    }

    _gain_cache.initializeEntry(u, benefit_aggregator, incident_edges_weight, penalty_aggregator);
  }

  HyperedgeWeight km1Gain(const HypernodeID u, PartitionID from, PartitionID to) const {
    unused(from);
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveFromBenefit(u) - moveToPenalty(u, to);
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, block weights and pin counts in part for
  // ! each hyperedge must be initialized explicitly here.
  void initializePartition() {
    tbb::parallel_invoke(
            [&] { initializeBlockWeights(); },
            [&] { initializePinCountInPart(); }
    );
  }

  bool isGainCacheInitialized() const {
    return _is_gain_cache_initialized;
  }

  // ! Initialize gain cache
  // ! NOTE: Requires that pin counts are already initialized and reflect the
  // ! current state of the partition
  void initializeGainCache() {
    // check whether part has been initialized
    ASSERT( _part_ids.size() == initialNumNodes()
            && std::none_of(nodes().begin(), nodes().end(),
                            [&](HypernodeID u) { return partID(u) == kInvalidPartition || partID(u) > k(); }) );


    auto aggregate_contribution_of_he_for_vertex =
      [&](const PartitionID block_of_u,
          const HyperedgeID he,
          vec<HyperedgeWeight>& l_move_from_benefit,
          HyperedgeWeight& incident_edges_weight,
          vec<HyperedgeWeight>& l_move_to_penalty) {
      HyperedgeWeight edge_weight = edgeWeight(he);
      for (const PartitionID block : connectivitySet(he)) {
        if (pinCountInPart(he, block) == 1) {
          l_move_from_benefit[block] += edge_weight;
        }
        l_move_to_penalty[block] += edge_weight;
      }
      incident_edges_weight += edge_weight;
    };

    // Gain calculation consist of two stages
    //  1. Compute gain of all low degree vertices sequential (with a parallel for over all vertices)
    //  2. Compute gain of all high degree vertices parallel (with a sequential for over all high degree vertices)
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mtp(_k, 0);
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mfb(_k, 0);
    std::mutex high_degree_vertex_mutex;
    parallel::scalable_vector<HypernodeID> high_degree_vertices;

    // Compute gain of all low degree vertices sequential (parallel for over all vertices)
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
      [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& l_move_to_penalty = ets_mtp.local();
        vec<HyperedgeWeight>& l_move_from_benefit = ets_mfb.local();
        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
          if ( nodeIsEnabled(u)) {
            if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD) {
              const PartitionID from = partID(u);
              HyperedgeWeight incident_edges_weight = 0;
              for (HyperedgeID he : incidentEdges(u)) {
                aggregate_contribution_of_he_for_vertex(from, he,
                  l_move_from_benefit, incident_edges_weight, l_move_to_penalty);
              }

              // Call to initializeEntry also resets l_move_from_benefit and l_move_to_penalty entries to zero
              _gain_cache.initializeEntry(u,l_move_from_benefit,incident_edges_weight,l_move_to_penalty);
            } else {
              // Collect high degree vertex for subsequent parallel gain computation
              std::lock_guard<std::mutex> lock(high_degree_vertex_mutex);
              high_degree_vertices.push_back(u);
            }
          }
        }
      });

    auto check_mfb_zero = [&](){
        for (const auto& mfbs : ets_mfb) {
            for (const auto&  mfb : mfbs) {
                if (mfb != 0) return false;
            }
        }
        return true;
    };
    unused(check_mfb_zero);
    HEAVY_REFINEMENT_ASSERT(check_mfb_zero());

    // Compute gain of all high degree vertices parallel (sequential for over all high degree vertices)
      vec<HyperedgeWeight> aggregated_penalties(_k,0);
      vec<HyperedgeWeight> aggregated_benefits(_k,0);
    for ( const HypernodeID& u : high_degree_vertices ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> ets_iew(0);
      const PartitionID from = partID(u);
      const HypernodeID degree_of_u = _hg->nodeDegree(u);
      tbb::parallel_for(tbb::blocked_range<HypernodeID>(ID(0), degree_of_u),
        [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& l_move_to_penalty = ets_mtp.local();
        vec<HyperedgeWeight>& l_move_from_benefit = ets_mfb.local();
        HyperedgeWeight& l_incident_edges_weight = ets_iew.local();
        size_t current_pos = r.begin();
        for ( const HyperedgeID& he : _hg->incident_nets_of(u, r.begin()) ) {
          aggregate_contribution_of_he_for_vertex(from, he,
            l_move_from_benefit, l_incident_edges_weight, l_move_to_penalty);
          ++current_pos;
          if ( current_pos == r.end() ) {
            break;
          }
        }
      });

      // Aggregate thread locals to compute overall gain of the high degree vertex
      const HyperedgeWeight incident_edges_weight = ets_iew.combine(std::plus<>());
      for (PartitionID p = 0; p < _k; ++p) {
          for (auto& l_move_from_benefit : ets_mfb) {
              aggregated_benefits[p] += l_move_from_benefit[p];
              l_move_from_benefit[p] = 0;
          }
          for (auto& l_move_to_penalty : ets_mfb) {
              aggregated_penalties[p] += l_move_to_penalty[p];
              l_move_to_penalty[p] = 0;
          }
      }

      // Call to initializeEntry internally resets aggregated_benefits and aggregated_penalties to 0
      _gain_cache.initializeEntry(u, aggregated_benefits, incident_edges_weight, aggregated_penalties);
    }

    _is_gain_cache_initialized = true;
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    _part_ids.assign(_part_ids.size(), kInvalidPartition, false);
    for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      for ( const PartitionID& block : connectivitySet(he) ) {
        _pins_in_part.setPinCountInPart(he, block, 0);
      }
      _connectivity_set.clear(he);
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

  // ! Recomputes the benefits for node u for each block (if a full benefit cache is used)
  vec<HyperedgeWeight> moveFromBenefitsRecomputed(const HypernodeID u) const {
      vec<HyperedgeWeight> benefits(_k,0);
      for (HyperedgeID e : incidentEdges(u)) {
          for (PartitionID p : _connectivity_set.connectivitySet(e)) {
              if (pinCountInPart(e, p) == 1) {
                  benefits[p] += edgeWeight(e);
              }
          }
      }
      return benefits;
  }

  // ! Recomputes the aggregated benefit for a block if the aggregated benefit-cache is used
  HyperedgeWeight moveFromBenefitRecomputedAggregated(const HypernodeID u) const {
    const PartitionID p = partID(u);
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) == 1) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  // ! Only for testing
  HyperedgeWeight moveToPenaltyRecomputed(const HypernodeID u, PartitionID p) const {
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) == 0) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  void recomputeMoveFromBenefit(const HypernodeID u) {
    auto recomputed_benefits = moveFromBenefitsRecomputed(u);
    _gain_cache.storeRecomputedMoveFromBenefits(u, recomputed_benefits);
  }

  // ! Only for testing
  bool checkTrackedPartitionInformation() {
    bool success = true;

    for (HyperedgeID e : edges()) {
      PartitionID expected_connectivity = 0;
      for (PartitionID i = 0; i < k(); ++i) {
        const HypernodeID actual_pin_count_in_part = pinCountInPart(e, i);
        if ( actual_pin_count_in_part != pinCountInPartRecomputed(e, i) ) {
          LOG << "Pin count of hyperedge" << e << "in block" << i << "=>" <<
              "Expected:" << V(pinCountInPartRecomputed(e, i)) << "," <<
              "Actual:" <<  V(pinCountInPart(e, i));
          success = false;
        }
        expected_connectivity += (actual_pin_count_in_part > 0);
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
        if (moveFromBenefit(u) != moveFromBenefitRecomputedAggregated(u) ) {
          LOG << "Move from benefit of hypernode" << u << "=>" <<
              "Expected:" << V(moveFromBenefitRecomputedAggregated(u)) << ", " <<
              "Actual:" <<  V(moveFromBenefit(u));
          success = false;
        }

        for (PartitionID i = 0; i < k(); ++i) {
          if (partID(u) != i) {
            if ( moveToPenalty(u, i) != moveToPenaltyRecomputed(u, i) ) {
              LOG << "Move to penalty of hypernode" << u << "in block" << i << "=>" <<
                  "Expected:" << V(moveToPenaltyRecomputed(u, i)) << ", " <<
                  "Actual:" <<  V(moveToPenalty(u, i));
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
    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    _connectivity_set.memoryConsumption(connectivity_set_node);

    parent->addChild("Part Weights", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", _pins_in_part.size_in_bytes());
    parent->addChild("Gain Cache", _gain_cache.size_in_bytes());
    parent->addChild("HE Ownership", sizeof(SpinLock) * _hg->initialNumNodes());

    size_t size_of_snapshots = 0;
    for (const auto& snapshot : _direct_conn_set_snapshots) {
      size_of_snapshots += snapshot.size_in_bytes();
    }
    for (const auto& snapshot : _pin_count_in_part_bitcopy_snapshots) {
      size_of_snapshots += snapshot->size_in_bytes();
    }
    parent->addChild("Connectivity Set Snapshots", size_of_snapshots);
  }

  // ####################### Extract Block #######################

  // ! Extracts a block of a partition as separate hypergraph.
  // ! It also returns a vertex-mapping from the original hypergraph to the sub-hypergraph.
  // ! If cut_net_splitting is activated, hyperedges that span more than one block (cut nets) are split, which is used for the connectivity metric.
  // ! Otherwise cut nets are discarded (cut metric).
  std::pair<Hypergraph, parallel::scalable_vector<HypernodeID> > extract(
          PartitionID block,
          bool cut_net_splitting,
          bool stable_construction_of_incident_edges) {
    ASSERT(block != kInvalidPartition && block < _k);

    // Compactify vertex ids
    parallel::scalable_vector<HypernodeID> hn_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    parallel::scalable_vector<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    HypernodeID num_hypernodes = 0;
    HypernodeID num_hyperedges = 0;
    tbb::parallel_invoke([&] {
      for ( const HypernodeID& hn : nodes() ) {
        if ( partID(hn) == block ) {
          hn_mapping[hn] = num_hypernodes++;
        }
      }
    }, [&] {
      for ( const HyperedgeID& he : edges() ) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          he_mapping[he] = num_hyperedges++;
        }
      }
    });

    // Extract plain hypergraph data for corresponding block
    using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
    HyperedgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.resize(num_hyperedges);
      doParallelForAllEdges([&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          ASSERT(he_mapping[he] < num_hyperedges);
          hyperedge_weight[he_mapping[he]] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[he_mapping[he]].push_back(hn_mapping[pin]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      doParallelForAllNodes([&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          hypernode_weight[hn_mapping[hn]] = nodeWeight(hn);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_hypergraph = HypergraphFactory::construct(num_hypernodes, num_hyperedges,
            edge_vector, hyperedge_weight.data(), hypernode_weight.data(), stable_construction_of_incident_edges);

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( partID(hn) == block ) {
        const HypernodeID extracted_hn = hn_mapping[hn];
        extracted_hypergraph.setCommunityID(extracted_hn, _hg->communityID(hn));
      }
    });
    return std::make_pair(std::move(extracted_hypergraph), std::move(hn_mapping));
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      tbb::parallel_invoke( [&] {
        parallel::parallel_free(_part_ids, _pin_count_update_ownership);
      }, [&] {
        parallel::free(_pins_in_part.data());
      }, [&] {
        _connectivity_set.freeInternalData();
      } );
    }
    _k = 0;
  }

 private:

  void applyPartWeightUpdates(vec<HypernodeWeight>& part_weight_deltas) {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  void initializeBlockWeights() {
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HypernodeWeight> pws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID pu = partID( u );
          const HypernodeWeight wu = nodeWeight( u );
          pws[pu] += wu;
        }
      }
      applyPartWeightUpdates(pws);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  void initializePinCountInPart() {
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if ( edgeIsEnabled(he) ) {
          for (const HypernodeID& pin : pins(he)) {
            ++pin_counts[partID(pin)];
          }

          for (PartitionID p = 0; p < _k; ++p) {
            assert(pinCountInPart(he, p) == 0);
            if (pin_counts[p] > 0) {
              _connectivity_set.add(he, p);
              _pins_in_part.setPinCountInPart(he, p, pin_counts[p]);
            }
            pin_counts[p] = 0;
          }
        }
      }
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()), assign);
  }

  HypernodeID pinCountInPartRecomputed(const HyperedgeID e, PartitionID p) const {
    HypernodeID pcip = 0;
    for (HypernodeID u : pins(e)) {
      if (partID(u) == p) {
        pcip++;
      }
    }
    return pcip;
  }

  // ! Updates pin count in part using a spinlock.
  template<typename DeltaFunc, typename GainCacheUpdateFunc>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID he,
                                                                    const PartitionID from,
                                                                    const PartitionID to,
                                                                    DeltaFunc&& delta_func,
                                                                    GainCacheUpdateFunc&& gain_cache_update_func) {
    ASSERT(he < _pin_count_update_ownership.size());
    _pin_count_update_ownership[he].lock();
    const HypernodeID pin_count_in_from_part_after = decrementPinCountInPartWithoutGainUpdate(he, from);
    const HypernodeID pin_count_in_to_part_after = incrementPinCountInPartWithoutGainUpdate(he, to);
    _pin_count_update_ownership[he].unlock();
    delta_func(he, edgeWeight(he), edgeSize(he), pin_count_in_from_part_after, pin_count_in_to_part_after);
    gain_cache_update_func(he, edgeWeight(he), pins(he), from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
  }

    // ! Updates pin count in part using a spinlock. In this variant for asynchronous uncoarsening concurrent
    // ! uncontractions and moves are allowed. For small edges, the gain cache update is performed inside the HE lock
    // ! and for large edges the pins of the hyperedge are stored within the pin count update lock in order to make the
    // ! gain cache update work on the right pins outside of the lock. Returns true if successful and false if the pin
    // ! count update lock could not be acquired.
  template<typename DeltaFunc, typename GainCacheUpdateFunc>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool asyncUpdatePinCountOfHyperedge(const HyperedgeID he,
                                                                         const PartitionID from,
                                                                         const PartitionID to,
                                                                         DeltaFunc&& delta_func,
                                                                         GainCacheUpdateFunc&& gain_cache_update_func) {
      ASSERT(he < _pin_count_update_ownership.size());
      if (!_pin_count_update_ownership[he].tryLock()) return false;
      const HypernodeID pin_count_in_from_part_after = decrementPinCountInPartWithoutGainUpdate(he, from);
      const HypernodeID pin_count_in_to_part_after = incrementPinCountInPartWithoutGainUpdate(he, to);
      size_t edge_size = edgeSize(he);

      if (_hg->totalSize(he) < _hg->snapshotEdgeSizeThreshold()) {
        gain_cache_update_func(he, edgeWeight(he), pins(he), from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
        _pin_count_update_ownership[he].unlock();
      } else {

        std::unique_ptr<IteratorRange<PinSnapshotIterator>> pins_snapshot = std::make_unique<IteratorRange<PinSnapshotIterator>>(PinSnapshotIterator::emptyPinIteratorRange());
        if (_gain_cache.arePinCountsAfterMoveThatTriggerUpdateForAllPins(pin_count_in_from_part_after, pin_count_in_to_part_after)) {
          pins_snapshot = std::make_unique<IteratorRange<PinSnapshotIterator>>(takePinsSnapshot(he));
        }
        _pin_count_update_ownership[he].unlock();
        gain_cache_update_func(he, edgeWeight(he), *pins_snapshot, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
      }

      delta_func(he, edgeWeight(he), edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
      return true;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.decrementPinCountInPart(e, p);
    if ( pin_count_after == 0 ) {
      _connectivity_set.remove(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.incrementPinCountInPart(e, p);
    if ( pin_count_after == 1 ) {
      _connectivity_set.add(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void lockHyperedge(const HyperedgeID he) {
//      _he_ownership[he].lock();
      _hg->acquireHyperedge(he);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void unlockHyperedge(const HyperedgeID he) {
//      _he_ownership[he].unlock();
      _hg->releaseHyperedge(he);
  }

  // ! Indicate wheater gain cache is initialized
  bool _is_gain_cache_initialized;

  // ! Number of blocks
  PartitionID _k = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_weights;

  // ! Current block IDs of the vertices
  Array< PartitionID > _part_ids;

  // ! For each hyperedge and each block, _pins_in_part stores the
  // ! number of pins in that block
  PinCountInPart _pins_in_part;

  // ! For each hyperedge, _connectivity_set stores the set of blocks that the hyperedge spans
  ConnectivitySets _connectivity_set;

  // ! Gain-cache that is kept up to date by uncontraction and move operations (with eventual correctness) and is used
  // ! to determine the best move
  GainCache _gain_cache;

  // ! In order to update the pin count of a hyperedge thread-safe, a thread must acquire
  // ! the ownership of a hyperedge via a CAS operation.
  Array<SpinLock> _pin_count_update_ownership;

  // ! Thread-local PartitionBitSets used for snapshots of connectivity sets in asynchronous uncoarsening
  tbb::enumerable_thread_specific<CompressedConnectivitySetSnapshot> _direct_conn_set_snapshots;
  tbb::enumerable_thread_specific<CompressedConnectivitySetSnapshot> _direct_parts_with_one_pin_snapshots;
  tbb::enumerable_thread_specific<std::unique_ptr<PinCountInPart::Snapshot>> _pin_count_in_part_bitcopy_snapshots;
    tbb::enumerable_thread_specific<std::unique_ptr<ConnectivitySets::Snapshot>> _connectivity_set_bitcopy_snapshots;
  tbb::enumerable_thread_specific<std::vector<HypernodeID>> _pins_snapshots;

  // ! Stats counters
  CAtomic<size_t> _num_stable_pins_seen;
  CAtomic<size_t> _num_volatile_pins_seen;

};

} // namespace ds
} // namespace mt_kahypar