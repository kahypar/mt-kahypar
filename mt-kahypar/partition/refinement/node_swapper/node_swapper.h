/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

class NodeSwapper {
 private:

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

public:
  explicit NodeSwapper(PartitionedHypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _nodes(hypergraph.initialNumNodes()),
    _in_queue(hypergraph.initialNumNodes()) {
    tbb::parallel_for(ID(0), hypergraph.initialNumNodes(),
      [&](const HypernodeID& hn) { _nodes[hn] = hn; });
  }

  NodeSwapper(const NodeSwapper&) = delete;
  NodeSwapper(NodeSwapper&&) = delete;

  NodeSwapper & operator= (const NodeSwapper &) = delete;
  NodeSwapper & operator= (NodeSwapper &&) = delete;

  HyperedgeWeight refine();

private:
  PartitionedHypergraph& _hg;
  const Context& _context;
  vec<HypernodeID> _nodes;
  kahypar::ds::FastResetFlagArray<> _in_queue;
};
}  // namespace kahypar
