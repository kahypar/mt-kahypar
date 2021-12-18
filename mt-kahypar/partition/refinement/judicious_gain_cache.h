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

  /*! \enum pqStatus
   *
   *  Return status of findNextMove()
   *  ok: move found
   *  rollback: disabled blocks have been enabled and caller should perform a rollback
   *  empty: no moves left
   */
  enum pqStatus {
    ok = 0,
    rollback = 1,
    empty = 2,
  };

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

  void initBlockPQ(PartitionedHypergraph& phg, HypernodeWeight from_load) {
    ASSERT(_blockPQ.empty());
    vec<std::pair<HyperedgeWeight, PartitionID>> remove_candidates;
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      const HyperedgeWeight load = phg.partLoad(i);
      if (static_cast<double>(load) / from_load > _context.refinement.judicious.block_disable_factor) {
        remove_candidates.push_back(std::make_pair(load, i));
        _blocks_enabled[i] = false;
      } else if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, _toPQs[i].topKey());
      }
    }
    if (_blockPQ.size() > static_cast<size_t>(_context.partition.k / 2)) {
      return;
    }
    // if too many blocks would be disabled, add some of them to the blockPQ by increasing load
    std::sort(remove_candidates.begin(), remove_candidates.end());
    for (size_t i = 0; i < remove_candidates.size() &&
         _blockPQ.size() <= static_cast<size_t>(_context.partition.k / 2); ++i) {
      const PartitionID candidate = remove_candidates[i].second;
      if (!_toPQs[candidate].empty()) {
        _blockPQ.insert(candidate, _toPQs[candidate].topKey());
        _blocks_enabled[candidate] = true;
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

    if (phg.k() < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
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

  pqStatus findNextMove(const PartitionedHypergraph& phg, Move& m) {
    ASSERT(!_blockPQ.contains(_active_part));
    const auto status = updatePQs();
    if (status != pqStatus::ok) {
      return status;
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
    return status;
  }

  void resetGainCache() {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _toPQs[i].clear();
      _blocks_enabled[i] = true;
    }
    _blockPQ.clear();
    _enable_all_blocks = false;
    _rebalancing = false;
    _target_parts.assign(_target_parts.size(), kInvalidPartition);
  }

  PartitionID getTargetPart(const HypernodeID v) const {
    return _target_parts[v];
  }

  void updateEnabledBlocks(const PartitionID to, const HyperedgeWeight from_load, const HyperedgeWeight to_load) {
    // only consider disabling block if not all blocks should be enabled, at least half the blocks are enabled and the block is larger enough
    if (!_enable_all_blocks
        && _blockPQ.size() > static_cast<size_t>(_context.partition.k / 2)
        && static_cast<double>(to_load) / from_load > _context.refinement.judicious.block_disable_factor) {
      _blocks_enabled[to] = false;
      if (_blockPQ.contains(to)) {
        _blockPQ.remove(to);
      }
    }
  }

  void setOnlyEnabledBlock(PartitionID p) {
    _enable_all_blocks = false;
    _blocks_enabled.assign(_blocks_enabled.size(), false);
    _blocks_enabled[p] = true;
    _rebalancing = true;
  }

private:

  bool blockIsEnabled(PartitionID p) {
    return _blocks_enabled[p] || _enable_all_blocks;
  }

  pqStatus updatePQs() {
    // first update the blockPQ
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i);
    }
    // then if it is empty try enabling all blocks
    if (_blockPQ.empty() && !_enable_all_blocks && !_rebalancing) {
      _enable_all_blocks = true;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
          updateOrRemoveToPQFromBlocks(i);
      }
      return _blockPQ.empty() ? pqStatus::empty : pqStatus::rollback;
    }
    return _blockPQ.empty() ? pqStatus::empty : pqStatus::ok;
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
/*
    if (to == kInvalidPartition) {
      return std::make_pair(to, std::numeric_limits<HyperedgeWeight>::min());
    }
*/
    ASSERT(to != kInvalidPartition);
    to_load_after += phg.weightOfDisabledEdges(u);
    Gain benefit = phg.moveFromBenefit(u) + phg.weightOfDisabledEdges(u);
    HyperedgeWeight from_load_after = from_load - benefit;
    // (Review Note) If to block is light enough, we only care about benefit.
    // light enough means part_load[to] + penalty < part_load[from], I think
    // but penalty and load can be a good tie-breaker because otherwise all moves look the same
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
  bool _enable_all_blocks = false;
  bool _rebalancing = false;
  PartitionID _active_part;
};
}
