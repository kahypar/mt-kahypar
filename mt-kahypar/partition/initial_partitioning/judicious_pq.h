/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2022 Noah Wahl <noah.wahl@student.kit.edu>
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
#include <mt-kahypar/partition/metrics.h>

namespace mt_kahypar {
class JudiciousPQ final {
public:
  using PriorityQueue = ds::ExclusiveHandleHeap<
      ds::Heap<HypernodeWeight, PartitionID, std::greater<>>>;

  explicit JudiciousPQ(const Context &context, HypernodeID num_nodes)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)),
        _disabled_blocks(context.partition.k, false) { }

  void insert(const PartitionedHypergraph &phg, const HypernodeID v) {
    const PartitionID pv = phg.partID(v);
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (i == pv)
        continue;
      Gain gain = 0;
      if (pv == kInvalidPartition) {
        gain = computeGainForInvalidFrom(phg, v, i);
      } else {
        gain = computeGainForValidFrom(phg, v, pv, i);
      }
      _toPQs[i].insert(v, gain);
    }
  }

  void initBlockPQ(const PartitionedHypergraph &phg, const PartitionID default_block) {
    ASSERT(_blockPQ.empty());
    _judicious_load = default_block == kInvalidPartition ? 0 : phg.partLoad(default_block);
    // ASSERT(_judicious_load == metrics::judiciousLoad(phg));
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, blockGain(phg, i));
      }
    }
  }

  void increaseGain(const PartitionedHypergraph &phg, const HypernodeID v,
                  const HyperedgeID he, const PartitionID to) {
    _toPQs[to].increaseKey(v, _toPQs[to].keyOf(v) - phg.edgeWeight(he));
  }

  void decreaseGain(const PartitionedHypergraph &phg, const HypernodeID v,
                  const HyperedgeID he, const PartitionID to) {
    _toPQs[to].decreaseKey(v, _toPQs[to].keyOf(v) + phg.edgeWeight(he));
  }

  bool findNextMove(const PartitionedHypergraph &phg, Move &m) {
    if (!updatePQs(phg)) {
      return false;
    }
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    const Gain gain = -_blockPQ.topKey();
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    for (auto &pq : _toPQs) {
      if (pq.contains(u)) {
        pq.remove(u);
      }
    }
    return true;
  }

  void updateJudiciousLoad(const PartitionedHypergraph& phg, const PartitionID from, const PartitionID to) {
    const HyperedgeWeight from_load = from != kInvalidPartition ? phg.partLoad(from) : 0;
    _judicious_load = std::max(_judicious_load, std::max(from_load, phg.partLoad(to)));
  }

  void disableBlock(const PartitionID p) {
    _disabled_blocks[p] = true;
  }

private:
  bool updatePQs(const PartitionedHypergraph &phg) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i, phg);
    }
    return !_blockPQ.empty();
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i,
                                    const PartitionedHypergraph &phg) {
    if (!_toPQs[i].empty()) {
      _blockPQ.insertOrAdjustKey(i, blockGain(phg, i));
    } else if ((_toPQs[i].empty() || _disabled_blocks[i]) && _blockPQ.contains(i)) {
      _blockPQ.remove(i);
    }
  }

  Gain computeGainForInvalidFrom(const PartitionedHypergraph &phg,
                                 const HypernodeID u, const PartitionID to) {
    Gain penalty = 0;
    for (const auto &he : phg.incidentEdges(u)) {
      const HyperedgeWeight pin_count_in_to = phg.pinCountInPart(he, to);
      if (pin_count_in_to == 0) {
        penalty += phg.edgeWeight(he);
      }
    }
    penalty += phg.weightOfDisabledEdges(u);
    return penalty;
  }

  Gain computeGainForValidFrom(const PartitionedHypergraph &phg,
                                 const HypernodeID u, const PartitionID from, const PartitionID to) {
    Gain penalty = 0;
    for (const auto &he : phg.incidentEdges(u)) {
      const HyperedgeWeight pin_count_in_to = phg.pinCountInPart(he, to);
      const HyperedgeWeight pin_count_in_from = phg.pinCountInPart(he, from);
      if (pin_count_in_to == 0) {
        penalty += phg.edgeWeight(he);
      }
      if (pin_count_in_from == 1) {
        penalty -= phg.edgeWeight(he);
      }
    }
    penalty += phg.weightOfDisabledEdges(u);
    return penalty;
  }

  Gain blockGain(const PartitionedHypergraph& phg, const PartitionID p) {
    if (_context.initial_partitioning.use_judicious_increase) {
      return std::max(phg.partLoad(p) + _toPQs[p].topKey() - _judicious_load, 0);
    } else if (_context.initial_partitioning.use_block_load_only) {
      return phg.partLoad(p);
    } else {
      return _toPQs[p].topKey();
    }
  }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  HyperedgeWeight _judicious_load;
  vec<bool> _disabled_blocks;
};
} // namespace mt_kahypar
