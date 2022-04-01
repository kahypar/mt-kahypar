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
#include "mt-kahypar/partition/refinement/judicious/pq_strategy.h"

namespace mt_kahypar {
class GreedyJudiciousInitialPartitioner {

  static constexpr bool debug = false;

 public:
  GreedyJudiciousInitialPartitioner(PartitionedHypergraph& phg,
                                    const Context& context) :
  _phg(phg),
  _context(context),
  _pq(context, phg.initialNumNodes()),
  _default_part(0),
  _neighbor_deduplicator(phg.initialNumNodes(), 0)
  { }

  GreedyJudiciousInitialPartitioner(const GreedyJudiciousInitialPartitioner&) = delete;
  GreedyJudiciousInitialPartitioner(GreedyJudiciousInitialPartitioner&&) = delete;

  void updateNeighbors(const Move& move) {
    for (HyperedgeID e : _edges_with_gain_changes) {
        for (HypernodeID v : _phg.pins(e)) {
          if (_neighbor_deduplicator[v] != _deduplication_time) {
            if (_phg.partID(v) == _default_part) {
              _pq.updateGain(_phg, v, move);
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

  void initialPartition() {
    for (const HypernodeID& hn : _phg.nodes()) {
      _phg.setNodePart(hn, _default_part);
      _pq.insert(_phg, hn);
    }
    _pq.setActivePart(_default_part);
    _pq.initBlockPQ();

    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        _edges_with_gain_changes.push_back(he);
      }
    };
    Move move;
    // TODO: other strategies, round robin, first fit?
    do {
      if (!_pq.findNextMove(_phg, move)) break;
      ASSERT(move.from == _default_part);
      _phg.changeNodePartWithGainCacheUpdate(move.node, move.from, move.to,
                                            std::numeric_limits<HypernodeWeight>::max(),
                                            []{}, delta_func);
      updateNeighbors(move);
    } while (_phg.partLoad(_default_part) > _phg.partLoad(move.to));
  }

 private:
  PartitionedHypergraph& _phg;
  const Context _context;
  JudiciousGainCache _pq;
  const PartitionID _default_part;
  vec<HyperedgeID> _edges_with_gain_changes;
  vec<HypernodeID> _neighbor_deduplicator;
  HypernodeID _deduplication_time = 1;
};
}
