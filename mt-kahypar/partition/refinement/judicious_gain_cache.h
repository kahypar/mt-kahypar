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

  void updateGain(const PartitionedHypergraph& phg, const HypernodeID v, const Move& move) {
    const PartitionID pv = phg.partID(v);
    const PartitionID designatedTargetV = _target_parts[v];
    if (designatedTargetV == kInvalidPartition || designatedTargetV == pv) { // no localized searches for now
      return;
    }
    Gain gain = 0;
    PartitionID newTarget = kInvalidPartition;

    if (phg.k() < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
      // moveToPenalty of designatedTargetV is affected.
      // and may now be greater than that of other blocks --> recompute full
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, v, pv, _parts);
    } else {
      // moveToPenalty of designatedTargetV is not affected.
      // only move.from and move.to may be better
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, v, pv, { designatedTargetV, move.from, move.to });
    }

    if (designatedTargetV == newTarget) {
      _toPQs[designatedTargetV].adjustKey(v, gain);
    } else if (designatedTargetV != kInvalidPartition) {
      _toPQs[designatedTargetV].remove(v);
      _toPQs[newTarget].insert(v, gain);
      _target_parts[v] = newTarget;

    }
  }

  bool findNextMove(const PartitionedHypergraph& phg, Move& m, bool& should_refiner_perform_rollback) {
    if (!updatePQs(should_refiner_perform_rollback)) {
      return false;
    }
    // this is pretty hacky... maybe come up with something better
    if (should_refiner_perform_rollback) {
      m.to = kInvalidPartition;
      return true;
    }

    const PartitionID to = _blockPQ.top();
    const HypernodeID u = _toPQs[to].top();
    /*const Gain estimated_gain = _toPQs[to].topKey();*/
    auto [_to, gain] = computeBestTargetBlock(phg, u, phg.partID(u), _parts);
    ASSERT(_blocks_enabled[to] || _enable_all_blocks);
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    _toPQs[to].deleteTop();  // blockPQ updates are done later, collectively.
    return true;
  }

  void resetGainCache() {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _toPQs[i].clear();
      _blocks_enabled[i] = true;
    }
    _blockPQ.clear();
    _enable_all_blocks = false;
    _target_parts.assign(_target_parts.size(), kInvalidPartition);
  }

  PartitionID getTargetPart(const HypernodeID v) const {
    return _target_parts[v];
  }

  void updateEnabledBlocks(const PartitionID to, const HyperedgeWeight from_load, const HyperedgeWeight to_load) {
    // only consider disabling block if not all blocks should be enabled, at least half the blocks are enabled and the block is larger enough
    if (!_enable_all_blocks
        && _blockPQ.size() > static_cast<size_t>(_context.partition.k / 2)
        && 1.f * to_load / from_load > _block_disable_factor) {
      _blocks_enabled[to] = false;
      _blockPQ.remove(to);
    }
  }

private:
  bool updatePQs(bool& should_refiner_perform_rollback) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (_blocks_enabled[i] || _enable_all_blocks) {
        updateOrRemoveToPQFromBlocks(i);
      }
    }
    if (_enable_all_blocks) {
      return !_blockPQ.empty();
    }
    // if not all blocks were enabled and no block is left, retry with all blocks
    if (_blockPQ.empty()) {
      _enable_all_blocks = true;
      should_refiner_perform_rollback = true;
      return updatePQs(should_refiner_perform_rollback);
    }
    return true;
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insertOrAdjustKey(i, _toPQs[i].topKey());
      } else if (_blockPQ.contains(i)) {
        _blockPQ.remove(i);
      }
  }

  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PartitionedHypergraph& phg,
                                                                 const HypernodeID u,
                                                                 const PartitionID from,
                                                                 const vec<PartitionID> parts) {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i : parts) {
      if (i != from) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
        if ((penalty < to_penalty || (penalty == to_penalty && to_weight < best_to_weight))) {
          to_penalty = penalty;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    const Gain gain = to != kInvalidPartition ? phg.moveFromBenefit(u) - to_penalty
                                              : std::numeric_limits<HyperedgeWeight>::min();
    return std::make_pair(to, gain);
  }

  const Context& _context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  vec<PartitionID> _target_parts;
  vec<PartitionID> _parts;
  vec<bool> _blocks_enabled;
  bool _enable_all_blocks = false;
  const double _block_disable_factor = 0.9;
};
}
