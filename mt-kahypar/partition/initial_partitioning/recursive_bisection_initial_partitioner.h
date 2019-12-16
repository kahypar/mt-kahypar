/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <limits>
#include <vector>

#include "tbb/parallel_invoke.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include <libkahypar.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace multilevel {
static inline void partition(Hypergraph& hypergraph, const Context& context, const bool top_level);
}  // namespace multilevel

template <typename TypeTraits>
class RecursiveBisectionInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;
  static constexpr bool kahypar_debug = false;
  static constexpr bool enable_heavy_assert = false;

  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

 public:
  RecursiveBisectionInitialPartitionerT(HyperGraph& hypergraph, const Context& context, const bool top_level) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level) { }

  RecursiveBisectionInitialPartitionerT(const RecursiveBisectionInitialPartitionerT&) = delete;
  RecursiveBisectionInitialPartitionerT(RecursiveBisectionInitialPartitionerT&&) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (const RecursiveBisectionInitialPartitionerT &) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (RecursiveBisectionInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    for ( const HypernodeID& hn : _hg.nodes() ) {
      _hg.setNodePart(hn, 0);
    }

    _hg.updateGlobalPartInfos();
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
};

template <typename TypeTraits>
PartitionID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using RecursiveBisectionInitialPartitioner = RecursiveBisectionInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
