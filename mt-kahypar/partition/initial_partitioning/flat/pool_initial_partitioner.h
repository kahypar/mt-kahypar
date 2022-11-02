/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
