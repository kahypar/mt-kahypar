/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/partition/context.h>

namespace mt_kahypar {
class JudiciousGainCache final {
public:
  using PriorityQueue = ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>>;

  explicit JudiciousGainCache(const Context& context, HypernodeID num_nodes) :
    _context(context),
    _toPQs(static_cast<size_t>(context.partition.k), PriorityQueue(num_nodes)),
    _blockPQ(static_cast<size_t>(context.partition.k)),
    _target_parts(num_nodes, kInvalidPartition),
    _parts(context.partition.k),
    _blocks_enabled(context.partition.k, true) {
    std::iota(_parts.begin(), _parts.end(), 0);
  }

  void insert(const PartitionedHypergraph& phg, const HypernodeID v) {
    const PartitionID pv = phg.partID(v);
    auto [target, gain] = computeBestTargetBlock(phg, v, pv, _parts);
    _toPQs[target].insert(v, gain);
    _target_parts[v] = target;
  }

  void initBlockPQ() {
    ASSERT(_blockPQ.empty());
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, _toPQs[i].topKey());
      }
    }
  }

  void setActivePart(const PartitionID part_id) {
    _active_part = part_id;
  }

  void updateGain(const PartitionedHypergraph& phg, const HypernodeID v, const Move& move) {
    const PartitionID pv = phg.partID(v);
    const PartitionID designatedTargetV = _target_parts[v];
    ASSERT(_toPQs[designatedTargetV].contains(v));
    if (designatedTargetV == kInvalidPartition) return;
    Gain gain = 0;
    PartitionID newTarget = kInvalidPartition;

    if (phg.k() < 4 || designatedTargetV == move.from || designatedTargetV == move.to || _rebalancing) {
      // moveToPenalty of designatedTargetV is affected.
      // and may now be greater than that of other blocks --> recompute full
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, v, pv, _parts);
    } else {
      // moveToPenalty of designatedTargetV is not affected.
      // only move.from and move.to may be better
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, v, pv, { designatedTargetV, move.to });
    }

    if (designatedTargetV == newTarget) {
      _toPQs[designatedTargetV].adjustKey(v, gain);
    } else {
      _toPQs[designatedTargetV].remove(v);
      _toPQs[newTarget].insert(v, gain);
      _target_parts[v] = newTarget;
    }
  }

  bool findNextMove(const PartitionedHypergraph& phg, Move& m) {
    ASSERT(!_blockPQ.contains(_active_part));
    if (!updatePQs()) {
      return false;
    }
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(to != _active_part);
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    ASSERT(phg.partID(u) == _active_part);
    auto [_to, gain] = computeBestTargetBlock(phg, u, phg.partID(u), _parts);
    ASSERT(_to != _active_part);
    ASSERT(blockIsEnabled(_to));
    m.node = u;
    m.from = phg.partID(u);
    m.to = _to;
    m.gain = gain;
    _toPQs[to].deleteTop();
    return true;
  }

  void resetGainCache() {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _toPQs[i].clear();
      _blocks_enabled[i] = true;
    }
    _blockPQ.clear();
    _rebalancing = false;
    _target_parts.assign(_target_parts.size(), kInvalidPartition);
  }

  void setOnlyEnabledBlock(PartitionID p) {
    _blocks_enabled.assign(_blocks_enabled.size(), false);
    _blocks_enabled[p] = true;
    _rebalancing = true;
  }

private:

  bool blockIsEnabled(PartitionID p) const {
    return _blocks_enabled[p];
  }

  bool updatePQs() {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i);
    }
    return !_blockPQ.empty();
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i) {
      if (!_toPQs[i].empty()) {
        if (blockIsEnabled(i)) {
          _blockPQ.insertOrAdjustKey(i, _toPQs[i].topKey());
        }
      } else if (_blockPQ.contains(i)) {
        _blockPQ.remove(i);
      }
  }

  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PartitionedHypergraph& phg,
                                                                 const HypernodeID u,
                                                                 const PartitionID from,
                                                                 const vec<PartitionID> parts) {
    const HyperedgeWeight from_load = phg.partLoad(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_load = std::numeric_limits<HyperedgeWeight>::max();
    HyperedgeWeight to_load_after = std::numeric_limits<HyperedgeWeight>::max();
    for (PartitionID i : parts) {
      if (i != from && blockIsEnabled(i)) {
        const HyperedgeWeight load = phg.partLoad(i);
        const HyperedgeWeight load_after = load + phg.moveToPenalty(u, i);
        if (load_after < to_load_after || (load_after == to_load_after && load < to_load)) {
          to = i;
          to_load = load;
          to_load_after = load_after;
        }
      }
    }
    ASSERT(to != kInvalidPartition);
    to_load_after += phg.weightOfDisabledEdges(u);
    Gain benefit = phg.moveFromBenefit(u) + phg.weightOfDisabledEdges(u);
    HyperedgeWeight from_load_after = from_load - benefit;
    const Gain gain = to_load_after < from_load_after ? benefit :
                                                        from_load - to_load_after;

    return std::make_pair(to, gain);
  }

  const Context& _context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  vec<PartitionID> _target_parts;
  vec<PartitionID> _parts;
  vec<bool> _blocks_enabled;
  bool _rebalancing = false;
  // only used for assertions
  PartitionID _active_part;
};
}
