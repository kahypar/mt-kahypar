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

#include <libkahypar.h>

#include "kahypar/partition/context.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
template <typename TypeTraits>
class DirectInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PoolInitialPartitioner = PoolInitialPartitionerT<TypeTraits>;

  static constexpr bool debug = false;

 public:
  DirectInitialPartitionerT(HyperGraph& hypergraph,
                            const Context& context,
                            const bool top_level,
                            const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  DirectInitialPartitionerT(const DirectInitialPartitionerT&) = delete;
  DirectInitialPartitionerT(DirectInitialPartitionerT&&) = delete;
  DirectInitialPartitionerT & operator= (const DirectInitialPartitionerT &) = delete;
  DirectInitialPartitionerT & operator= (DirectInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    if ( _top_level ) {
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    PoolInitialPartitioner& initial_partitioner = *new(tbb::task::allocate_root())
      PoolInitialPartitioner(_hg, _context, _task_group_id);
    tbb::task::spawn_root_and_wait(initial_partitioner);

    if ( _top_level ) {
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

using DirectInitialPartitioner = DirectInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
