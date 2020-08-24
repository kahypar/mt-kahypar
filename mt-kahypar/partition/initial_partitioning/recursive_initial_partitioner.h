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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"

namespace mt_kahypar {

/*!
 * RECURSIVE INITIAL PARTITIONER
 * For reason of simplicity we assume in the following description of the algorithm that
 * the number of threads p and the number of blocks k is a power of 2 and p = k. The recursive
 * initial partitioner is invoked, if the number of vertices is 2 * c * p (where c is our
 * contraction limit multiplier).
 * The recursive initial partitioner starts by performing parallel coarsening with p threads
 * until c * p vertices are reached. Afterwards, the hypergraph is copied and the hypergraphs
 * are recursively coarsened with p / 2 threads each. Once p = 1 (and the contraction limit is 2 * c)
 * we initially bisect the hypergraph in two blocks. After initial partitioning each thread uncontracts
 * its hypergraph (and performs refinement) until 4 * c hypernodes are rechead. Afterwards, we choose the
 * best partition of both recursions and further bisect each block of the partition to obtain a 4-way
 * partition and continue uncontraction with 2 threads until 8 * c hypernodes. This is repeated until
 * we obtain a k-way partition of the hypergraph.
 * Note, the recursive initial partitioner is written in TBBNumaArena continuation style. The TBBNumaArena continuation
 * style is especially useful for recursive patterns. Each task defines its continuation task. A continuation
 * task defines how computation should continue, if all its child tasks are completed. As a consequence,
 * tasks can be spawned without waiting for their completion, because the continuation task is automatically
 * invoked if all child tasks are terminated. Therefore, no thread will waste CPU time while waiting for
 * their recursive tasks to complete.
 *
 * Implementation Details
 * ----------------------
 * The recursive initial partitioner starts by spawning the root RecursiveTask. The RecursiveTask spawns
 * two RecursiveChildTask. Within such a task the hypergraph is copied and coarsened to the next desired contraction limit.
 * Once that contraction limit is reached the RecursiveChildTask spawns again one RecursiveTask. Once the RecursiveTask
 * of a RecursiveChildTask terminates, the RecursiveChildContinuationTask starts and uncontracts the hypergraph to
 * its original size (and also performs refinement). Once both RecursiveChildTask of a RecursiveTask terminates, the
 * RecursiveContinuationTask starts and chooses the best partition of both recursions and spawns for each block
 * a BisectionTask. The BisectionTask performs a initial partition call to bisect exactly one block of the current
 * partition. Once all BisectionTasks terminates, the BisectionContinuationTask starts and applies all bisections to the
 * current hypergraph.
 */
class RecursiveInitialPartitioner: public IInitialPartitioner {
 private:

  static constexpr bool enable_heavy_assert = false;



 public:
  RecursiveInitialPartitioner(PartitionedHypergraph& hypergraph,
                               const Context& context,
                               const bool top_level,
                               const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  RecursiveInitialPartitioner(const RecursiveInitialPartitioner&) = delete;
  RecursiveInitialPartitioner(RecursiveInitialPartitioner&&) = delete;
  RecursiveInitialPartitioner & operator= (const RecursiveInitialPartitioner &) = delete;
  RecursiveInitialPartitioner & operator= (RecursiveInitialPartitioner &&) = delete;

 private:
  void initialPartitionImpl() final ;

 private:
  PartitionedHypergraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

}  // namespace mt_kahypar
