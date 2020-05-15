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
 * an improvement and are immediatly reverted. However, applying them directly to global
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

 public:
  DeltaPartitionedHypergraph(const PartitionID k) :
    _k(k),
    _phg(nullptr),
    _part_weights_delta(k, 0),
    _part_ids_delta(),
    _pins_in_part_delta(),
    _move_to_penalty_delta(),
    _move_from_benefit_delta() { }

  // REVIEW NOTE why do we delete copy assignment/construction? wouldn't it be useful to make a copy, e.g. for initial partitioning
  DeltaPartitionedHypergraph(const DeltaPartitionedHypergraph&) = delete;
  DeltaPartitionedHypergraph & operator= (const DeltaPartitionedHypergraph &) = delete;

  DeltaPartitionedHypergraph(DeltaPartitionedHypergraph&& other) = default;
  DeltaPartitionedHypergraph & operator= (DeltaPartitionedHypergraph&& other) = default;

  ~DeltaPartitionedHypergraph() = default;

  void setPartitionedHypergraph(PartitionedHypergraph* phg) {
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

  // ####################### Hyperedge Information #######################

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->edgeSize(e);
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
                      DeltaFunc&& delta_func) {
    ASSERT(_phg);
    assert(partID(u) == from);
    assert(from != to);
    const HypernodeWeight wu = _phg->nodeWeight(u);
    if ( partWeight(to) + wu <= max_weight_to ) {
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += wu;
      _part_weights_delta[from] -= wu;
      for ( const HyperedgeID& he : _phg->incidentEdges(u) ) {
        const HypernodeID pin_count_in_from_part_after =
          decrementPinCountInPartWithGainUpdate(he, from);
        const HypernodeID pin_count_in_to_part_after =
          incrementPinCountInPartWithGainUpdate(he, to);
        delta_func(he, _phg->edgeWeight(he), _phg->edgeSize(he),
          pin_count_in_from_part_after, pin_count_in_to_part_after);
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
                      const HypernodeWeight max_weight_to) {
    return changeNodePart(u, from, to, max_weight_to, NoOpDeltaFunc());
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
    const uint64_t* pin_count_delta = _pins_in_part_delta.get_if_contained(e * _k + p);
    return std::max(static_cast<uint64_t>(_phg->pinCountInPart(e, p)) +
      ( pin_count_delta ? *pin_count_delta : 0 ), static_cast<uint64_t>(0));
  }

  // ! Returns the sum of all edges incident to u, where u is the last remaining
  // ! pin in its block
  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    ASSERT(_phg);
    const HyperedgeWeight* move_from_benefit_delta =
      _move_from_benefit_delta.get_if_contained(u);
    return _phg->moveFromBenefit(u) + ( move_from_benefit_delta ? *move_from_benefit_delta : 0 );
  }

  // ! Returns the sum of all edges incident to u, where p is not part of
  // ! their connectivity set.
  HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight* move_to_penalty_delta =
      _move_to_penalty_delta.get_if_contained(u * _k + p);
    return _phg->moveToPenalty(u, p) + ( move_to_penalty_delta ? *move_to_penalty_delta : 0 );
  }

  // ! Clears all deltas applied to the partitioned hypergraph
  void clear() {
    // O(k)
    _part_weights_delta.assign(_k, 0);
    // Constant Time
    _part_ids_delta.clear();
    _pins_in_part_delta.clear();
    _move_to_penalty_delta.clear();
    _move_from_benefit_delta.clear();
  }

  PartitionID k() const {
    return _k;
  }

 private:


  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = std::max(static_cast<uint64_t>(
      _phg->pinCountInPart(e, p)) + --_pins_in_part_delta[e * _k + p], static_cast<uint64_t>(0));
    if (pin_count_after == 1) {
      const HyperedgeWeight we = _phg->edgeWeight(e);
      for (HypernodeID u : _phg->pins(e)) {
        if (partID(u) == p) {
          _move_from_benefit_delta[u] += we;
        }
      }
    } else if (pin_count_after == 0) {
      const HyperedgeWeight we = _phg->edgeWeight(e);
      for (HypernodeID u : _phg->pins(e)) {
        _move_to_penalty_delta[u * _k + p] += we;
      }
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = std::max(static_cast<uint64_t>(
      _phg->pinCountInPart(e, p)) + ++_pins_in_part_delta[e * _k + p], static_cast<uint64_t>(0));
    if (pin_count_after == 1) {
      const HyperedgeWeight we = _phg->edgeWeight(e);
      for (HypernodeID u : _phg->pins(e)) {
        _move_to_penalty_delta[u * _k + p] -= we;
      }
    } else if (pin_count_after == 2) {
      const HyperedgeWeight we = _phg->edgeWeight(e);
      for (HypernodeID u : _phg->pins(e)) {
        if (partID(u) == p) {
          _move_from_benefit_delta[u] -= we;
        }
      }
    }
    return pin_count_after;
  }

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
  DynamicSparseMap<size_t, uint64_t> _pins_in_part_delta;

  // ! Stores the delta of each locally touched move to penalty entry
  // ! relative to the _move_to_penalty member in '_phg'
  DynamicSparseMap<size_t, HyperedgeWeight> _move_to_penalty_delta;

  // ! Stores the delta of each locally touched move from benefit entry
  // ! relative to the _move_from_benefit member in '_phg'
  DynamicSparseMap<HypernodeID, HyperedgeWeight> _move_from_benefit_delta;
};

} // namespace ds
} // namespace mt_kahypar