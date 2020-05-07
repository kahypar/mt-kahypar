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


template <typename PartitionedHypergraph = Mandatory>
class DeltaPartitionedHypergraph {
private:

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

  template<typename F>
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to,
                      F&& report_success) {
    ASSERT(_phg);
    assert(partID(u) == from);
    assert(from != to);
    const HypernodeWeight wu = _phg->nodeWeight(u);
    if ( partWeight(to) + wu <= max_weight_to ) {
      report_success();
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += wu;
      _part_weights_delta[from] -= wu;
      for ( const HyperedgeID& he : _phg->incidentEdges(u) ) {
        decrementPinCountInPartWithGainUpdate(he, from);
        incrementPinCountInPartWithGainUpdate(he, to);
      }
      return true;
    } else {
      return false;
    }
  }

  PartitionID partID(const HypernodeID u) const {
    ASSERT(_phg);
    if ( _part_ids_delta.contains(u) ) {
      return _part_ids_delta.get(u);
    } else {
      return _phg->partID(u);
    }
  }

  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    return _phg->partWeight(p) + _part_weights_delta[p];
  }

  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    uint64_t pin_count_delta = 0;
    if ( _pins_in_part_delta.contains(e * _k + p) ) {
      pin_count_delta = _pins_in_part_delta.get(e * _k + p);
    }
    return std::max(static_cast<uint64_t>(_phg->pinCountInPart(e, p)) +
      pin_count_delta, static_cast<uint64_t>(0));
  }

  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    ASSERT(_phg);
    uint64_t move_from_benefit_delta = 0;
    if ( _move_from_benefit_delta.contains(u) ) {
      move_from_benefit_delta = _move_from_benefit_delta.get(u);
    }
    return _phg->moveFromBenefit(u) + move_from_benefit_delta;
  }

  HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    uint64_t move_to_penalty_delta = 0;
    if ( _move_to_penalty_delta.contains(u * _k + p) ) {
      move_to_penalty_delta = _move_to_penalty_delta.get(u * _k + p);
    }
    return _phg->moveToPenalty(u, p) + move_to_penalty_delta;
  }

  void clear() {
    _part_weights_delta.assign(_k, 0);
    _part_ids_delta.clear();
    _pins_in_part_delta.clear();
    _move_to_penalty_delta.clear();
    _move_from_benefit_delta.clear();
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

  const PartitionID _k;
  PartitionedHypergraph* _phg;
  vec< HypernodeWeight > _part_weights_delta;
  DynamicSparseMap<HypernodeID, PartitionID> _part_ids_delta;
  DynamicSparseMap<size_t, uint64_t> _pins_in_part_delta;
  DynamicSparseMap<size_t, HyperedgeWeight> _move_to_penalty_delta;
  DynamicSparseMap<HypernodeID, HyperedgeWeight> _move_from_benefit_delta;
};

} // namespace ds
} // namespace mt_kahypar