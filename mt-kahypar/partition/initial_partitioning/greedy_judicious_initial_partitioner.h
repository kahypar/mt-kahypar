/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2022 Noah Wahl <noah.wahl@student.kit.edu
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/judicious_pq.h"
#include <random>

namespace mt_kahypar {
class GreedyJudiciousInitialPartitioner {

  static constexpr bool debug = false;

public:
  GreedyJudiciousInitialPartitioner(
      PartitionedHypergraph &phg, const Context &context, const size_t seed,
      GreedyJudiciousInitialPartitionerStats &stats)
      : _phg(phg), _context(context),
        _pq(context, phg.initialNumNodes(), seed, stats),
        _preassign_nodes(context.initial_partitioning.preassign_nodes),
        _stats(stats) {
    _default_part = _preassign_nodes ? 0 : -1;
  }

  GreedyJudiciousInitialPartitioner(const GreedyJudiciousInitialPartitioner &) =
      delete;
  GreedyJudiciousInitialPartitioner(GreedyJudiciousInitialPartitioner &&) =
      delete;

  void initialPartition() {
    _phg.resetPartition();
    if (_preassign_nodes) {
      HighResClockTimepoint assign_start = std::chrono::high_resolution_clock::now();
      for (const HypernodeID &hn : _phg.nodes()) {
        _phg.setNodePart(hn, _default_part);
      }
      HighResClockTimepoint assign_stop = std::chrono::high_resolution_clock::now();
      double assign_time = std::chrono::duration<double>(assign_stop - assign_start).count();
      DBG << V(assign_time);
    }
    _pq.init(_phg, _default_part);

    Move move;

    auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following
      // situations.
      if (_preassign_nodes && pin_count_in_from_part_after == 1) {
        for (HypernodeID v : _phg.pins(he)) {
          // being in _default_part means the node is unassigned if
          // _preassign_nodes == false
          if (_phg.partID(v) == _default_part) {
            _pq.increaseBenefit(_phg, v, he);
            if constexpr (debug) {
              _stats.update_hist[_phg.edgeSize(he)]++;
            }
          }
        }
      }
      if (pin_count_in_to_part_after == 1) {
        for (HypernodeID v : _phg.pins(he)) {
          if (_phg.partID(v) == _default_part) {
            _pq.increaseGain(_phg, v, he, move.to);
            if constexpr (debug) {
              _stats.update_hist[_phg.edgeSize(he)]++;
            }
          }
        }
      }
    };

    auto success_func = [&]() {
      if (_preassign_nodes &&
          _phg.partLoad(_default_part) <= _phg.partLoad(move.to)) {
        _pq.disableBlock(move.to);
      }
      if constexpr (debug) {
        _stats.num_moved_nodes++;
      }
    };
    while (_pq.getNextMove(_phg, move)) {
      ASSERT(move.from == _default_part);
      if (_preassign_nodes) {
        _phg.changeNodePart(move.node, move.from, move.to, std::numeric_limits<HyperedgeWeight>::max(), success_func, delta_func);
      } else {
        _phg.setNodePart(move.node, move.to);
        success_func();
        for (const auto he : _phg.incidentEdges(move.node)) {
          delta_func(he, 0, 0, 0, _phg.pinCountInPart(he, move.to));
        }
      }

      _pq.updateJudiciousLoad(_phg, move.from, move.to);
      ASSERT(_stats.gain_sequence.back() == move.gain);
    }
    ASSERT(std::all_of(
        _phg.nodes().begin(), _phg.nodes().end(),
        [&](const auto &hn) { return _phg.partID(hn) != kInvalidPartition; }));
  }

private:
  PartitionedHypergraph &_phg;
  const Context _context;
  JudiciousPQ _pq;
  PartitionID _default_part;
  const bool _preassign_nodes;
  GreedyJudiciousInitialPartitionerStats &_stats;
};
} // namespace mt_kahypar
