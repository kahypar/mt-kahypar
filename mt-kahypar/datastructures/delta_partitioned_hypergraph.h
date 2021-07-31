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

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Special hypergraph data structure for the FM algorithm.
 * In order to perform frequent moves on the hypergraph without affecting other
 * local searches one can wrap this hypergraph around a partitioned hypergraph and
 * perform those moves locally on this hypergraph. The changes are not reflected
 * in the global partitioned hypergraph. Instead, the data structure stores deltas
 * relative to the global partitioned hypergraph for each internal member such that
 * by applying those deltas to the original partitioned hypergraph it would reflect the
 * global state of the partitioned hypergraph after all local moves would be applied
 * to it.
 * The rationale behind this is that the majority of local searches to not yield to
 * an improvement and are immediately reverted. However, applying them directly to global
 * partitioned hypergraph would affect other local searches running concurrently, which build
 * upon that state. This special partitioned hypergraph allows a local search to hide its
 * current search state from other searches in a space efficient manner.
 */
template <typename PartitionedHypergraph = Mandatory>
class DeltaPartitionedHypergraph {
 private:

  using HypernodeIterator = typename PartitionedHypergraph::HypernodeIterator;
  using HyperedgeIterator = typename PartitionedHypergraph::HyperedgeIterator;
  using IncidenceIterator = typename PartitionedHypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename PartitionedHypergraph::IncidentNetsIterator;

  using GainCacheDelta = typename PartitionedHypergraph::GainCacheDelta;

 public:
  static constexpr bool supports_connectivity_set = false;
  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = PartitionedHypergraph::HIGH_DEGREE_THRESHOLD;

  explicit DeltaPartitionedHypergraph(const PartitionID k) :
    _k(k),
    _phg(nullptr),
    _part_weights_delta(k, 0),
    _part_ids_delta(),
    _pins_in_part_delta(),
    _gain_cache_deltas(k, [&](const HypernodeID hn){return partID(hn);}) { }

  DeltaPartitionedHypergraph(const DeltaPartitionedHypergraph&) = delete;
  DeltaPartitionedHypergraph & operator= (const DeltaPartitionedHypergraph &) = delete;

  DeltaPartitionedHypergraph(DeltaPartitionedHypergraph&& other)  noexcept = default;
  DeltaPartitionedHypergraph & operator= (DeltaPartitionedHypergraph&& other)  noexcept = default;

  ~DeltaPartitionedHypergraph() = default;

  void setPartitionedHypergraph(PartitionedHypergraph* phg) {
    ASSERT(phg);
    if (!phg) LOG << "phg is nullptr in DeltaPHG::setPartitionedHypergraph()!";

    if constexpr (GainCacheDelta::notify_about_updates_on_phg) {
      if (_phg != phg) {
        phg->addDeltaObserver(this);
      }
      _gain_cache_deltas.resize(phg->initialNumNodes(), phg->initialNumEdges());
    }

    _phg = phg;
  }

  // ####################### Iterators #######################

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    ASSERT(_phg);
    return _phg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    ASSERT(_phg);
    return _phg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->pins(e);
  }

  // ####################### Hypernode Information #######################

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->nodeWeight(u);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->nodeDegree(u);
  }

  bool nodeIsEnabled(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->nodeIsEnabled(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->edgeSize(e);
  }

  HyperedgeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->edgeWeight(e);
  }

  // ####################### Partition Information #######################

  // ! Changes the block of hypernode u from 'from' to 'to'.
  // ! Move is successful, if it is not violating the balance
  // ! constraint specified by 'max_weight_to'.
  template<typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to,
                      DeltaFunc&& delta_func,
                      const bool async_context = false) {
    ASSERT(_phg);
    assert(partID(u) == from);
    assert(from != to);
    ASSERT(_phg->nodeIsEnabled(u));
    const HypernodeWeight wu = _phg->nodeWeight(u);
    if ( partWeight(to) + wu <= max_weight_to ) {
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += wu;
      _part_weights_delta[from] -= wu;
      auto inc_edges = _phg->incidentEdges(u);
      for ( const HyperedgeID& he : inc_edges ) {
        if (async_context) _phg->lockHyperedgePinCountLock(he);
        const HypernodeID pin_count_in_from_part_after = decrementPinCountInPart(he, from);
        const HypernodeID pin_count_in_to_part_after = incrementPinCountInPart(he, to);
        delta_func(he, _phg->edgeWeight(he), _phg->edgeSize(he), pin_count_in_from_part_after, pin_count_in_to_part_after);
        if (async_context) _phg->unlockHyperedgePinCountLock(he);
      }
      return true;
    } else {
      return false;
    }
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to,
                      const bool async_context = false) {
    return changeNodePart(u, from, to, max_weight_to, NoOpDeltaFunc(), async_context);
  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         const PartitionID from,
                                         const PartitionID to,
                                         const HypernodeWeight max_weight_to,
                                         const bool async_context = false) {
    auto delta_gain_func = [&]( HyperedgeID he, HyperedgeWeight edge_weight,
                                HypernodeID, HypernodeID pcip_from, HypernodeID pcip_to ) {
      gainCacheUpdate(he, edge_weight, pins(he), from, pcip_from, to, pcip_to);
    };
    return changeNodePart(u, from, to, max_weight_to, delta_gain_func, async_context);
  }

  template <typename PinIteratorT>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID he, const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                         const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                         const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
    _gain_cache_deltas.updateForDeltaMove(he, we, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
  }

  // ! Returns the block of hypernode u
  PartitionID partID(const HypernodeID u) const {
    ASSERT(_phg);
    const PartitionID* part_id = _part_ids_delta.get_if_contained(u);
    return part_id ? *part_id : _phg->partID(u);
  }

  // ! Returns the total weight of block p
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    return _phg->partWeight(p) + _part_weights_delta[p];
  }

  // ! Returns the number of pins of hyperedge e in block p
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    const int32_t* pin_count_delta = _pins_in_part_delta.get_if_contained(e * _k + p);
    return std::max(static_cast<int32_t>(_phg->pinCountInPart(e, p)) +
      ( pin_count_delta ? *pin_count_delta : 0 ), 0);
  }

  // ! Returns the sum of all edges incident to u, where u is the last remaining
  // ! pin in its block
  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->moveFromBenefit(u) + _gain_cache_deltas.benefitDelta(u);
  }

  // ! Returns the sum of all edges incident to u, where p is not part of
  // ! their connectivity set.
  HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    return _phg->moveToPenalty(u, p) + _gain_cache_deltas.penaltyDelta(u, p);
  }

  Gain km1Gain(const HypernodeID u, const PartitionID from, const PartitionID to) const {
    unused(from);
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveFromBenefit(u) - moveToPenalty(u, to);
  }

  void initializeGainCacheEntry(const HypernodeID u, vec<Gain>& benefit_aggregator, vec<Gain>& penalty_aggregator) {
    _phg->initializeGainCacheEntry(u, benefit_aggregator, penalty_aggregator);
  }

  // ! Clears all deltas applied to the partitioned hypergraph
  void clear() {
    // O(k)
    _part_weights_delta.assign(_k, 0);
    // Constant Time
    _part_ids_delta.clear();
    _pins_in_part_delta.clear();
    _gain_cache_deltas.clear();
  }

  void dropMemory() {
    if (!_memory_dropped) {
      _memory_dropped = true;
      _part_ids_delta.freeInternalData();
      _pins_in_part_delta.freeInternalData();
      _gain_cache_deltas.dropMemory();
    }
  }

  size_t combinedMemoryConsumption() const {
    return _pins_in_part_delta.size_in_bytes()
           + _gain_cache_deltas.benefit_size_in_bytes()
           + _gain_cache_deltas.penalty_size_in_bytes()
           + _part_ids_delta.size_in_bytes();
  }

  PartitionID k() const {
    return _k;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* delta_phg_node = parent->addChild("Delta Partitioned Hypergraph");
    utils::MemoryTreeNode* part_weights_node = delta_phg_node->addChild("Delta Part Weights");
    part_weights_node->updateSize(_part_weights_delta.capacity() * sizeof(HypernodeWeight));
    utils::MemoryTreeNode* part_ids_node = delta_phg_node->addChild("Delta Part IDs");
    part_ids_node->updateSize(_part_ids_delta.size_in_bytes());
    utils::MemoryTreeNode* pins_in_part_node = delta_phg_node->addChild("Delta Pins In Part");
    pins_in_part_node->updateSize(_pins_in_part_delta.size_in_bytes());
    utils::MemoryTreeNode* move_from_benefit_node = delta_phg_node->addChild("Delta Move From Benefit");
    move_from_benefit_node->updateSize(_gain_cache_deltas.benefit_size_in_bytes());
    utils::MemoryTreeNode* move_to_penalty_node = delta_phg_node->addChild("Delta Move To Penalty");
    move_to_penalty_node->updateSize(_gain_cache_deltas.penalty_size_in_bytes());
  }

  void markActive() {
    ASSERT(!_active);
   _active = true;
  }

  void markInactive() {
    ASSERT(_active);
    _active = false;
  }

  template<typename PinIteratorT>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void updateDeltasForMoveOnUnderlyingPHG(const HyperedgeID he, const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                                                    const PartitionID from, const HypernodeID pin_count_in_from_part_after_in_delta, const HypernodeID pin_count_in_from_part_after_in_underlying,
                                                    const PartitionID to, const HypernodeID pin_count_in_to_part_after_in_delta, const HypernodeID pin_count_in_to_part_after_in_underlying) {
    if (!_active) return;
    _gain_cache_deltas.updateForMoveOnUnderlyingPHG(he, we, pins,
                                                    from, pin_count_in_from_part_after_in_delta, pin_count_in_from_part_after_in_underlying,
                                                    to, pin_count_in_to_part_after_in_delta, pin_count_in_to_part_after_in_underlying);
  }

    // ! Variant for call in pin count update edge lock so pin count in delta part is queried on demand
    template<typename PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateDeltasForMoveOnUnderlyingPHG(const HyperedgeID he, const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                                            const PartitionID from,  const HypernodeID pin_count_in_from_part_after_in_underlying,
                                            const PartitionID to,  const HypernodeID pin_count_in_to_part_after_in_underlying) {
      if (!_active) return;
      _gain_cache_deltas.updateForMoveOnUnderlyingPHG(he, we, pins,
                                                      from, pinCountInPart(he, from), pin_count_in_from_part_after_in_underlying,
                                                      to, pinCountInPart(he, to), pin_count_in_to_part_after_in_underlying);
    }

    template<typename PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateDeltasForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, HypernodeID v, IteratorRange<PinIteratorT> pins,
                                    const PartitionID block, const HypernodeID pin_count_in_part_after_in_delta,
                                    const HypernodeID pin_count_in_part_after_in_underlying) {
      if (!_active) return;
      _gain_cache_deltas.updateForUncontractCaseOne(he, we, v, pins, block, pin_count_in_part_after_in_delta, pin_count_in_part_after_in_underlying);
    }

    // ! Variant for call in pin count update edge lock so pin count in delta part is queried on demand
    template<typename PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateDeltasForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, HypernodeID v, IteratorRange<PinIteratorT> pins,
                                          const PartitionID block,
                                          const HypernodeID pin_count_in_part_after_in_underlying) {
      if (!_active) return;
      _gain_cache_deltas.updateForUncontractCaseOne(he, we, v, pins, block, pinCountInPart(he, block), pin_count_in_part_after_in_underlying);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateDeltasForUncontractCaseTwo(const HyperedgeID he, const HypernodeID u, const HypernodeID v,
                                    const PartitionID block) {
      if (!_active) return;
      _gain_cache_deltas.updateForUncontractCaseTwo(he, u, v, block);
    }

 private:

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    return std::max(static_cast<int32_t>(
      _phg->pinCountInPart(e, p)) + --_pins_in_part_delta[e * _k + p], static_cast<int32_t>(0));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    return std::max(static_cast<int32_t>(
      _phg->pinCountInPart(e, p)) + ++_pins_in_part_delta[e * _k + p], static_cast<int32_t>(0));
  }

  bool _memory_dropped = false;
  bool _active = false;

  // ! Number of blocks
  const PartitionID _k;

  // ! Partitioned hypergraph where all deltas are stored relative to
  PartitionedHypergraph* _phg;

  // ! Delta for block weights
  vec< HypernodeWeight > _part_weights_delta;

  // ! Stores for each locally moved node, its new block id
  DynamicSparseMap<HypernodeID, PartitionID> _part_ids_delta;

  // ! Stores the delta of each locally touched pin count entry
  // ! relative to the _pins_in_part member in '_phg'
  DynamicSparseMap<size_t, int32_t> _pins_in_part_delta;

  // ! Stores the delta of each locally touched move from benefit entry
  // ! relative to the _move_from_benefit member in '_phg'
  GainCacheDelta _gain_cache_deltas;
};


} // namespace ds
} // namespace mt_kahypar