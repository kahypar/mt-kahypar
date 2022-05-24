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
#include "mt-kahypar/datastructures/judicious_partitioned_hypergraph.h"
#include "mt-kahypar/partition/initial_partitioning/judicious_ip_commons.h"
#include <algorithm>
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

namespace mt_kahypar {

static constexpr bool debug = false;

class JudiciousPQ final {
public:
  using PriorityQueue =
      ds::ExclusiveHandleHeap<ds::Heap<std::pair<HypernodeWeight, size_t>,
                                       PartitionID, std::greater<>, 16>>;
  using JudiciousPartitionedHypergraph = ds::JudiciousPartitionedHypergraph;

  explicit JudiciousPQ(const Context &context, const HypernodeID num_nodes,
                       const size_t seed,
                       GreedyJudiciousInitialPartitionerStats &stats,
                       const GreedyJudiciousInitialPartitionerConfig &config)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)),
        _part_loads(static_cast<size_t>(context.partition.k)),
        _disabled_blocks(context.partition.k, false), _g(seed), _stats(stats),
        _config(config) {
    if (_config.preassign_nodes) {
      _move_from_benefit.resize(num_nodes, 0);
    }
  }

  void init(const JudiciousPartitionedHypergraph &phg,
            const PartitionID default_part) {
    HighResClockTimepoint init_start =
        std::chrono::high_resolution_clock::now();
    vec<Gain> penalties(phg.initialNumNodes(), 0);
    for (const auto &he : phg.edges()) {
      if (_config.preassign_nodes && phg.edgeSize(he) == 1) {
        for (const auto &v : phg.pins(he)) {
          _move_from_benefit[v] += phg.edgeWeight(he);
        }
      }
      for (const auto &v : phg.pins(he)) {
        penalties[v] += phg.edgeWeight(he);
      }
    }
    for (const auto &v : phg.nodes()) {
      const Gain penalty = penalties[v] + phg.weightOfDisabledEdges(v);
      if (_config.preassign_nodes) {
        _move_from_benefit[v] += phg.weightOfDisabledEdges(v);
      }
      const size_t tag = _config.random_selection ? _g() : penalty;
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        if (i == default_part)
          continue;
        _toPQs[i].insert(v, std::make_pair(penalty, tag));
      }
    }
    initBlockPQ(phg);
    HighResClockTimepoint init_stop = std::chrono::high_resolution_clock::now();
    double init_time =
        std::chrono::duration<double>(init_stop - init_start).count();
    DBG << V(init_time);
  }

  void increaseGain(const HypernodeID v, const HyperedgeWeight w,
                    const PartitionID to) {
    auto key = _toPQs[to].keyOf(v);
    _toPQs[to].increaseKey(v, std::make_pair(key.first - w, key.second));
  }

  void increaseBenefit(const JudiciousPartitionedHypergraph &phg,
                       const HypernodeID v, const HyperedgeID he) {
    ASSERT(_config.preassign_nodes);
    _move_from_benefit[v] += phg.edgeWeight(he);
  }

  bool getNextMove(const JudiciousPartitionedHypergraph &phg, Move &move) {
    if (!updatePQs(phg)) {
      return false;
    }
    return findNextMove(phg, move);
  }

  void updateJudiciousLoad(const JudiciousPartitionedHypergraph &phg,
                           const PartitionID from, const PartitionID to) {
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
    if constexpr (debug) {
      _stats.gain_sequence.push_back(judicious_load_before -
                                     judicious_load_after);
    }
  }

  void disableBlock(const PartitionID p) {
    _disabled_blocks[p] = true;
    _num_disabled_blocks++;
  }

private:
  bool findNextMove(const JudiciousPartitionedHypergraph &phg, Move &m) {
    ASSERT(!_blockPQ.empty());
    const PartitionID to = _blockPQ.top();
    ASSERT(_blockPQ.topKey() == blockGain(phg, to));
    ASSERT(!_toPQs[to].empty() && !_disabled_blocks[to]);
    const HypernodeID u = _toPQs[to].top();
    Gain gain = 0;
    if (_config.preassign_nodes) {
      gain = calculateGainWithPreassignment(phg, to);
    } else {
      gain = _config.use_judicious_increase
                 ? -_blockPQ.topKey().first
                 : std::min(_part_loads.topKey() -
                                (phg.partLoad(to) + _toPQs[to].topKey().first),
                            0);
    }
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    _toPQs[to].deleteTop();
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (i != to && i != m.from) {
        _toPQs[i].remove(u);
      }
    }
    updateOrRemoveToPQFromBlocks(to, phg);
    return true;
  }

  void initBlockPQ(const JudiciousPartitionedHypergraph &phg) {
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

  bool updatePQs(const JudiciousPartitionedHypergraph &phg) {
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      updateOrRemoveToPQFromBlocks(i, phg);
    }
    return !_blockPQ.empty() &&
           _num_disabled_blocks <
               static_cast<size_t>(_context.partition.k + 1 / 2);
  }

  void updateOrRemoveToPQFromBlocks(const PartitionID i,
                                    const JudiciousPartitionedHypergraph &phg) {
    if (!_toPQs[i].empty() && !_disabled_blocks[i]) {
      _blockPQ.adjustKey(i, blockGain(phg, i));
    } else if (_blockPQ.contains(i)) {
      _blockPQ.remove(i);
    }
  }

  std::pair<Gain, size_t> blockGain(const JudiciousPartitionedHypergraph &phg,
                                    const PartitionID p) {
    Gain gain = 0;
    if (_config.use_judicious_increase) {
      if (_config.preassign_nodes) {
        gain = -calculateGainWithPreassignment(phg, p);
      } else {
        gain = std::max(phg.partLoad(p) + _toPQs[p].topKey().first -
                            _part_loads.topKey(),
                        0);
      }
    } else if (_config.use_block_load_only) {
      gain = phg.partLoad(p);
    } else {
      gain = _toPQs[p].topKey().first;
    }
    return std::make_pair(gain, _toPQs[p].topKey().second);
  }

  Gain calculateGainWithPreassignment(const JudiciousPartitionedHypergraph &phg,
                                      const PartitionID p) const {
    ASSERT(_config.preassign_nodes);
    const HyperedgeWeight load_of_first = _part_loads.topKey();
    if (load_of_first == phg.partLoad(p))
      return -_toPQs[p].topKey().first;
    const HypernodeID u = _toPQs[p].top();
    const PartitionID from = phg.partID(u);
    const HyperedgeWeight from_load_after =
        phg.partLoad(from) - _move_from_benefit[u];
    const HyperedgeWeight to_load_after =
        phg.partLoad(p) + _toPQs[p].topKey().first;
    const HyperedgeWeight load_of_second = _part_loads.keyOfSecond();
    if (_part_loads.top() == from) {
      if (from_load_after > to_load_after && from_load_after > load_of_second)
        return _move_from_benefit[u];
      else if (from_load_after == to_load_after &&
               from_load_after >= load_of_second)
        return 0;
      else
        return load_of_first - std::max(to_load_after, load_of_second);
    } else
      return std::min(load_of_first - to_load_after, 0);
  }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>>
      _part_loads;
  vec<bool> _disabled_blocks;
  size_t _num_disabled_blocks = 0;
  vec<HyperedgeWeight> _move_from_benefit;
  std::mt19937 _g;
  GreedyJudiciousInitialPartitionerStats &_stats;
  const GreedyJudiciousInitialPartitionerConfig &_config;
};
} // namespace mt_kahypar
