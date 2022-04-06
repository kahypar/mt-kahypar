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
  GreedyJudiciousInitialPartitioner(PartitionedHypergraph &phg,
                                    const Context &context)
      : _phg(phg), _context(context), _pq(context, phg.initialNumNodes()),
        _neighbor_deduplicator(phg.initialNumNodes(), 0),
        _preassign_nodes(context.initial_partitioning.preassign_nodes),
        _random_selection(context.initial_partitioning.random_selection),
        _part_loads(static_cast<size_t>(context.partition.k)) {
    _default_part = _preassign_nodes ? 0 : -1;
  }

  GreedyJudiciousInitialPartitioner(const GreedyJudiciousInitialPartitioner &) =
      delete;
  GreedyJudiciousInitialPartitioner(GreedyJudiciousInitialPartitioner &&) =
      delete;

  void initialPartition(const size_t seed) {
    _phg.resetPartition();
    for (const HypernodeID &hn : _phg.nodes()) {
      if (_preassign_nodes) {
        _phg.setNodePart(hn, _default_part);
      }
      _pq.insert(_phg, hn);
    }
    if (_preassign_nodes) {
      _phg.initializePartition();
      _phg.initializeGainCache();
    }
    _pq.initBlockPQ(_phg);

    _part_loads.clear();
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      _part_loads.insert(i, _phg.partLoad(i));
    }

    auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following
      // situations.
      if (!_preassign_nodes && pin_count_in_to_part_after == 1) {
        _edges_with_gain_changes.push_back(he);
      } else if (pin_count_in_from_part_after == 0 ||
                 pin_count_in_from_part_after == 1 ||
                 pin_count_in_to_part_after == 1 ||
                 pin_count_in_to_part_after == 2) {
        _edges_with_gain_changes.push_back(he);
      }
    };
    Move move;
    std::mt19937 g(seed);
    // TODO: other strategies, round robin, first fit?
    do {
      if (!getNextMove(move, g))
        break;
      ASSERT(move.from == _default_part);
      if (_preassign_nodes) {
        _phg.changeNodePartWithGainCacheUpdate(
            move.node, move.from, move.to,
            std::numeric_limits<HypernodeWeight>::max(), [] {}, delta_func);
        _part_loads.adjustKey(move.from, _phg.partLoad(move.from));
      } else {
        _phg.setNodePart(move.node, move.to);
        for (const auto &he : _phg.incidentEdges(move.node)) {
          if (_phg.pinCountInPart(he, move.to) == 1) {
            _edges_with_gain_changes.push_back(he);
          }
        }
      }
      _part_loads.adjustKey(move.to, _phg.partLoad(move.to));
      updateNeighbors(move.to);
    } while (!_preassign_nodes ||
             _phg.partLoad(_default_part) > _phg.partLoad(move.to));
    ASSERT(std::all_of(
        _phg.nodes().begin(), _phg.nodes().end(),
        [&](const auto &hn) { return _phg.partID(hn) != kInvalidPartition; }));
  }

private:
  void updateNeighbors(const PartitionID to) {
    for (HyperedgeID e : _edges_with_gain_changes) {
      for (HypernodeID v : _phg.pins(e)) {
        if (_neighbor_deduplicator[v] != _deduplication_time) {
          if (_phg.partID(v) == _default_part) {
            _pq.updateGain(_phg, v, e, to);
          }
          _neighbor_deduplicator[v] = _deduplication_time;
        }
      }
    }
    _edges_with_gain_changes.clear();
    if (++_deduplication_time == 0) {
      _neighbor_deduplicator.assign(_neighbor_deduplicator.size(), 0);
      _deduplication_time = 1;
    }
  }

  Move chooseRandomMove(vec<Move> &moves, std::mt19937 &g) {
    std::uniform_int_distribution<> distrib(0, moves.size() - 1);
    return moves[distrib(g)];
  }

  bool getNextMove(Move &move, std::mt19937 &g) {
    vec<Move> potential_moves;
    if (!_random_selection) {
      return _pq.findNextMove(_phg, move);
    }
    while (_pq.findNextMove(_phg, move)) {
      if (potential_moves.empty() || move.gain == potential_moves[0].gain) {
        potential_moves.push_back(move);
      } else {
        _pq.insert(_phg, move.node);
        break;
      }
    }
    if (_random_selection && potential_moves.size() == 0)
      return false;
    if (potential_moves.size() == 1)
      move = potential_moves[0];
    if (potential_moves.size() > 1) {
      move = chooseRandomMove(potential_moves, g);
      for (const auto &m : potential_moves) {
        if (m.node != move.node) {
          _pq.insert(_phg, m.node);
        }
      }
    }
    return true;
  }

private:
  PartitionedHypergraph &_phg;
  const Context _context;
  JudiciousPQ _pq;
  PartitionID _default_part;
  vec<HyperedgeID> _edges_with_gain_changes;
  vec<HypernodeID> _neighbor_deduplicator;
  HypernodeID _deduplication_time = 1;
  const bool _preassign_nodes = false;
  const bool _random_selection = false;
  ds::ExclusiveHandleHeap<
      ds::Heap<HyperedgeWeight, PartitionID, std::greater<HyperedgeWeight>>>
      _part_loads;
};
} // namespace mt_kahypar
