/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <queue>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"

namespace mt_kahypar {
template <template <typename> class GainPolicy>
class Rebalancer {
 private:
  using GainCalculator = GainPolicy<PartitionedHypergraph>;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static constexpr Gain MIN_PQ_GAIN_THRESHOLD = 5;

public:

  struct MoveGainComparator {
    bool operator()(const Move& lhs, const Move& rhs) {
      return lhs.gain > rhs.gain || (lhs.gain == rhs.gain && lhs.node < rhs.node);
    }
  };

  using MovePQ = std::priority_queue<Move, vec<Move>, MoveGainComparator>;

  struct IndexedMovePQ {
    explicit IndexedMovePQ(const size_t idx) :
      idx(idx),
      pq() { }

    size_t idx;
    MovePQ pq;
  };

  explicit Rebalancer(PartitionedHypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _gain(context),
    _part_weights(_context.partition.k) { }

  Rebalancer(const Rebalancer&) = delete;
  Rebalancer(Rebalancer&&) = delete;

  Rebalancer & operator= (const Rebalancer &) = delete;
  Rebalancer & operator= (Rebalancer &&) = delete;

  void rebalance(Metrics& best_metrics);


  vec<Move> repairEmptyBlocks();

private:

  template<typename F>
  bool moveVertex(const HypernodeID hn, const Move& move, const F& objective_delta) {
    ASSERT(_hg.partID(hn) == move.from);
    const PartitionID from = move.from;
    const PartitionID to = move.to;
    const HypernodeWeight node_weight = _hg.nodeWeight(hn);
    if ( from != to ) {
      // Before moving, we ensure that the block we move the vertex to does
      // not become overloaded
      _part_weights[to] += node_weight;
      if ( _part_weights[to] <= _context.partition.max_part_weights[to] ) {
        if ( _hg.changeNodePart(hn, from, to, objective_delta) ) {
          DBG << "Moved vertex" << hn << "from block" << from << "to block" << to
              << "with gain" << move.gain;
          _part_weights[from] -= node_weight;
          return true;
        }
      }
      _part_weights[to] -= node_weight;
    }
    return false;
  }

  PartitionedHypergraph& _hg;
  const Context& _context;
  GainCalculator _gain;
  parallel::scalable_vector<AtomicWeight> _part_weights;
};

using Km1Rebalancer = Rebalancer<Km1Policy>;
using CutRebalancer = Rebalancer<CutPolicy>;
}  // namespace kahypar
