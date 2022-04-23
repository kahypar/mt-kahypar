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
#include <algorithm>
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

namespace mt_kahypar {

  static constexpr bool debug = true;

  struct GreedyJudiciousInitialPartitionerStats {
    size_t num_moved_nodes = 0;
    vec<Gain> gain_sequence;

    void print() {
      if (!debug) {
        return;
      }
      ASSERT(num_moved_nodes == gain_sequence.size());
      LOG << V(num_moved_nodes);
      for (const auto i : gain_sequence) {
        LLOG << i;
      }
    }
  };

class JudiciousPQ final {
public:
  using PriorityQueue = ds::ExclusiveHandleHeap<
      ds::Heap<std::pair<HypernodeWeight, size_t>, PartitionID, std::greater<>>>;

  explicit JudiciousPQ(const Context &context, const HypernodeID num_nodes, const size_t seed, GreedyJudiciousInitialPartitionerStats& stats)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)),
        _part_loads(static_cast<size_t>(context.partition.k)),
        _disabled_blocks(context.partition.k, false),
        _g(seed),
        _stats(stats) { }

  void init(const PartitionedHypergraph& phg, const PartitionID default_part) {
    vec<Gain> penalties(phg.initialNumNodes(), 0);
    for (const auto& he : phg.edges()) {
      for (const auto& v : phg.pins(he)) {
        penalties[v] += phg.edgeWeight(he);
      }
    }
    for (const auto &v : phg.nodes()) {
      const Gain penalty = penalties[v] + phg.weightOfDisabledEdges(v);
      const size_t tag = _context.initial_partitioning.random_selection ? _g() : penalty;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        if (i == default_part)
          continue;
        _toPQs[i].insert(v, std::make_pair(penalty, tag));
      }
    }
    initBlockPQ(phg);
  }

  void increaseGain(const PartitionedHypergraph &phg, const HypernodeID v,
                  const HyperedgeID he, const PartitionID to) {
    auto key = _toPQs[to].keyOf(v);
    _toPQs[to].increaseKey(v, std::make_pair(key.first - phg.edgeWeight(he), key.second));
  }

  void decreaseGain(const PartitionedHypergraph &phg, const HypernodeID v,
                  const HyperedgeID he, const PartitionID to) {
    auto key = _toPQs[to].keyOf(v);
    _toPQs[to].decreaseKey(v, std::make_pair(key.first + phg.edgeWeight(he), key.second));
  }

  bool getNextMove(const PartitionedHypergraph& phg, Move &move, const PartitionID default_part) {
    while (findNextMove(phg, move)) {
      if (phg.partID(move.node) == default_part) {
        return true;
      }
    }
    return false;
  }

  void updateJudiciousLoad(const PartitionedHypergraph& phg, const PartitionID from, const PartitionID to) {
    const HyperedgeWeight judicious_load_before = _part_loads.topKey();
    _part_loads.adjustKey(to, phg.partLoad(to));
    if (from != kInvalidPartition) {
      _part_loads.adjustKey(from, phg.partLoad(from));
    }
    const HyperedgeWeight judicious_load_after = _part_loads.topKey();
    ASSERT([&]() {
      HyperedgeWeight max_load = 0;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        max_load = std::max(max_load, phg.partLoad(i));
      }
      return max_load == _part_loads.topKey();
    }());
    _stats.gain_sequence.push_back(judicious_load_before - judicious_load_after);
  }

  void disableBlock(const PartitionID p) {
    _disabled_blocks[p] = true;
    _num_disabled_blocks++;
  }

private:
  bool findNextMove(const PartitionedHypergraph &phg, Move &m) {
    if (!updatePQs(phg)) {
      return false;
    }
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(!_toPQs[to].empty());
    const HypernodeID u = _toPQs[to].top();
    const Gain gain = _context.initial_partitioning.use_judicious_increase ? -_blockPQ.topKey().first : std::min(_part_loads.topKey() - (phg.partLoad(to) + _toPQs[to].topKey().first), 0);
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    _toPQs[to].deleteTop();
    return true;
  }

  void initBlockPQ(const PartitionedHypergraph &phg) {
    ASSERT(_blockPQ.empty());
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, phg.partLoad(i));
    }
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (!_toPQs[i].empty()) {
        _blockPQ.insert(i, blockGain(phg, i));
      }
    }
  }

  bool updatePQs(const PartitionedHypergraph &phg) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i, phg);
    }
    return !_blockPQ.empty()
      && _num_disabled_blocks < static_cast<size_t>(_context.partition.k + 1 / 2);
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i,
                                    const PartitionedHypergraph &phg) {
    if (!_toPQs[i].empty() && !_disabled_blocks[i]) {
      _blockPQ.insertOrAdjustKey(i, blockGain(phg, i));
    } else if (_blockPQ.contains(i)) {
      _blockPQ.remove(i);
    }
  }

  std::pair<Gain, size_t> blockGain(const PartitionedHypergraph& phg, const PartitionID p) {
    Gain gain = 0;
    if (_context.initial_partitioning.use_judicious_increase) {
      gain = std::max(phg.partLoad(p) + _toPQs[p].topKey().first - _part_loads.topKey(), 0);
    } else if (_context.initial_partitioning.use_block_load_only) {
      gain = phg.partLoad(p);
    } else {
      gain = _toPQs[p].topKey().first;
    }
    return std::make_pair(gain, _toPQs[p].topKey().second);
  }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>> _part_loads;
  vec<bool> _disabled_blocks;
  size_t _num_disabled_blocks = 0;
  std::mt19937 _g;
  GreedyJudiciousInitialPartitionerStats& _stats;
};
} // namespace mt_kahypar
