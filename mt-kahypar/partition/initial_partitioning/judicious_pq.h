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
class JudiciousPQ final {
public:
  using PriorityQueue = ds::ExclusiveHandleHeap<
      ds::Heap<HypernodeWeight, PartitionID, std::greater<>>>;

  explicit JudiciousPQ(const Context &context, HypernodeID num_nodes)
      : _context(context), _toPQs(static_cast<size_t>(context.partition.k),
                                  PriorityQueue(num_nodes)),
        _blockPQ(static_cast<size_t>(context.partition.k)),
        _part_loads(static_cast<size_t>(context.partition.k)),
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

  void initBlockPQ(const PartitionedHypergraph &phg) {
    ASSERT(_blockPQ.empty());
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, phg.partLoad(i));
    }
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _blockPQ.insert(i, blockGain(phg, i));
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

  // Not removing the node from all PQs leads to significantly worse IP results, not exactly sure why. This means we have to delete and reinsert a lot of moves...
  bool getNextMove(PartitionedHypergraph& phg, Move &move, std::mt19937 &g) {
    vec<Move> potential_moves;
    if (!_context.initial_partitioning.random_selection) {
      return findNextMove(phg, move);
    }
    while (findNextMove(phg, move)) {
      if (potential_moves.empty() || move.gain == potential_moves[0].gain) {
        potential_moves.push_back(move);
      } else {
        insert(phg, move.node);
        break;
      }
    }
    if (_context.initial_partitioning.random_selection && potential_moves.size() == 0) {
      return false;
    } else if (potential_moves.size() == 1) {
      move = potential_moves[0];
    } else if (potential_moves.size() > 1) {
      move = chooseRandomMove(potential_moves, g);
      for (const auto &m : potential_moves) {
        if (m.node != move.node) {
          insert(phg, m.node);
        }
      }
    }
    ASSERT(std::all_of(potential_moves.begin(), potential_moves.end(), [&](const auto &m) {
      if (m.node == move.node) {
        return true;
      }
      for (PartitionID i = 0; i < _context.partition.k; ++i) {
        if (!_toPQs[i].contains(m.node)) {
          return false;
        }
      }
      return true;
    }));
    return true;
  }

  /* TODO: move this to getNextMove when doing correct gains there <11-04-22, @noahares> */
  void updateJudiciousLoad(const PartitionedHypergraph& phg, const PartitionID from, const PartitionID to) {
    _part_loads.adjustKey(to, phg.partLoad(to));
    if (from != kInvalidPartition) {
      _part_loads.adjustKey(from, phg.partLoad(from));
    }
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
    const Gain gain = _blockPQ.topKey();
    m.node = u;
    m.from = phg.partID(u);
    m.to = to;
    m.gain = gain;
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (_toPQs[i].contains(u)) {
        _toPQs[i].remove(u);
      }
    }
    return true;
  }

  Move chooseRandomMove(vec<Move> &moves, std::mt19937 &g) {
    std::uniform_int_distribution<> distrib(0, moves.size() - 1);
    return moves[distrib(g)];
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

  // Review Note: if not bipartitioning and code is too slow, calculate gain to all blocks at the same time using connectivity set
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
      return std::max(phg.partLoad(p) + _toPQs[p].topKey() - _part_loads.topKey(), 0);
    } else if (_context.initial_partitioning.use_block_load_only) {
      return phg.partLoad(p);
    } else {
      return _toPQs[p].topKey();
    }
  }

  const Context &_context;
  vec<PriorityQueue> _toPQs;
  PriorityQueue _blockPQ;
  ds::ExclusiveHandleHeap<ds::MaxHeap<HypernodeWeight, PartitionID>> _part_loads;
  vec<bool> _disabled_blocks;
  size_t _num_disabled_blocks = 0;
};
} // namespace mt_kahypar
