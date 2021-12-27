/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
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

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"


namespace mt_kahypar {



// IP algorithm and random seed
using IPTaskList = vec< std::tuple<InitialPartitioningAlgorithm, int, int> >;

class PoolInitialPartitionerContinuation : public tbb::task {
public:
  PoolInitialPartitionerContinuation(PartitionedHypergraph& hypergraph,
                                     const Context& context);

  tbb::task* execute() override ;

  InitialPartitioningDataContainer _ip_data;
  const Context& _context;
  vec<IPTaskList> _ip_task_lists;
};


void spawn_initial_partitioner(PoolInitialPartitionerContinuation& continuation_task );


/*!
 * The pool initial partitioner spawns for each initial partitioning run and algorithm
 * exactly one initial partitioning task. The number of initial partitions computed during
 * an invocation of the pool initial partitioner is exactly ( num IP runs ) * (num IP algos).
 * The best partition is applied to the hypergraph.
 */
class PoolInitialPartitioner : public tbb::task {

  static constexpr bool debug = false;

 public:
  PoolInitialPartitioner(PartitionedHypergraph& hypergraph,
                         const Context& context);

  tbb::task* execute() override ;

 private:
  PartitionedHypergraph& _hg;
  const Context& _context;
};
} // namespace mt_kahypar
