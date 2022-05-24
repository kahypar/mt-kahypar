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
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/partition/context.h>

namespace mt_kahypar {
class JudiciousGainCache final {
public:
  using PriorityQueue =
      ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>>;

  explicit JudiciousGainCache(const Context &context,
                              const HypernodeID num_nodes)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)),
        _target_parts(num_nodes, kInvalidPartition) {}

  void init(const PartitionedHypergraph &phg,
            const vec<HypernodeID> &refinement_nodes, const PartitionID from) {
    _from = from;
    for (const HypernodeID v : refinement_nodes) {
      ASSERT(phg.partID(v) == from);
      auto [target, gain] = computeBestTargetBlock(phg, v);
      _toPQs[target].insert(v, gain);
      _target_parts[v] = target;
    }
    initBlockPQ(phg);
  }

  void updateGain(const PartitionedHypergraph &phg, const HypernodeID v) {
    const PartitionID designatedTargetV = _target_parts[v];
    ASSERT(_toPQs[designatedTargetV].contains(v));
    ASSERT(designatedTargetV != kInvalidPartition);
    auto [newTarget, gain] = computeBestTargetBlock(phg, v);
    if (designatedTargetV == newTarget) {
      _toPQs[designatedTargetV].adjustKey(v, gain);
    } else {
      _toPQs[designatedTargetV].remove(v);
      _toPQs[newTarget].insert(v, gain);
      _target_parts[v] = newTarget;
    }
  }

  bool findNextMove(const PartitionedHypergraph &phg, Move &m) {
    ASSERT(!_blockPQ.contains(_from));
    if (!updatePQs(phg)) {
      return false;
    }
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    ASSERT(phg.partID(u) == _from);
    const Gain gain = _blockPQ.topKey();
    ASSERT(gain == blockGain(phg, to));
    m.node = u;
    m.from = _from;
    m.to = to;
    m.gain = gain;
    _toPQs[to].deleteTop();
    return true;
  }

  void reset() {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _toPQs[i].clear();
    }
    _blockPQ.clear();
    _target_parts.assign(_target_parts.size(), kInvalidPartition);
  }

private:
  void initBlockPQ(const PartitionedHypergraph &phg) {
    ASSERT(_blockPQ.empty());
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, blockGain(phg, i));
      }
    }
  }

  bool updatePQs(const PartitionedHypergraph &phg) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(phg, i);
    }
    return !_blockPQ.empty();
  }

  void updateOrRemoveToPQFromBlocks(const PartitionedHypergraph &phg,
                                    const PartitionID i) {
    if (!_toPQs[i].empty()) {
      _blockPQ.insertOrAdjustKey(i, blockGain(phg, i));
    } else if (_blockPQ.contains(i)) {
      _blockPQ.remove(i);
    }
  }

  std::pair<PartitionID, HyperedgeWeight>
  computeBestTargetBlock(const PartitionedHypergraph &phg,
                         const HypernodeID u) {
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_load = std::numeric_limits<HyperedgeWeight>::max();
    HyperedgeWeight to_load_after = std::numeric_limits<HyperedgeWeight>::max();
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (i != _from) {
        const HyperedgeWeight load = phg.partLoad(i);
        const HyperedgeWeight load_after = load + phg.moveToPenalty(u, i);
        if (load_after < to_load_after ||
            (load_after == to_load_after && load < to_load)) {
          to = i;
          to_load = load;
          to_load_after = load_after;
        }
      }
    }
    ASSERT(to != kInvalidPartition);
    const Gain gain = phg.moveFromBenefit(u) - phg.moveToPenalty(u, to);

    return std::make_pair(to, gain);
  }

  // optimistic gain of the top node of the block; independent of other block
  // loads
  Gain blockGain(const PartitionedHypergraph &phg, const PartitionID to) {
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    const Gain benefit = phg.moveFromBenefit(u) + phg.weightOfDisabledEdges(u);
    const Gain penalty =
        phg.moveToPenalty(u, to) + phg.weightOfDisabledEdges(u);
    const HyperedgeWeight from_load_after = phg.partLoad(_from) - benefit;
    const HyperedgeWeight to_load_after = phg.partLoad(to) + penalty;
    return to_load_after <= from_load_after
               ? benefit
               : phg.partLoad(_from) - to_load_after;
  }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  vec<PartitionID> _target_parts;
  PartitionID _from;
};
} // namespace mt_kahypar
