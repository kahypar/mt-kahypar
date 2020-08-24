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
 * RECURSIVE BISECTION INITIAL PARTITIONER
 * The recursive bisection initial partitioner starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks are reached.
 * Note, the recursive bisection initial partitioner is written in TBBNumaArena continuation style. The TBBNumaArena
 * continuation style is especially useful for recursive patterns. Each task defines its continuation
 * task. A continuation task defines how computation should continue, if all its child tasks are completed.
 * As a consequence, tasks can be spawned without waiting for their completion, because the continuation
 * task is automatically invoked if all child tasks are terminated. Therefore, no thread will waste CPU
 * time while waiting for their recursive tasks to complete.
 *
 * Implementation Details
 * ----------------------
 * The recursive bisection initial partitioner starts by spawning the root RecursiveMultilevelBisectionTask. The RecursiveMultilevelBisectionTask
 * spawns a MultilevelBisectionTask that bisects the hypergraph (multilevel-fashion). Afterwards, the MultilevelBisectionContinuationTask continues
 * and applies the bisection to the hypergraph and spawns two RecursiveBisectionChildTasks. Both are responsible for exactly one block of
 * the partition. The RecursiveBisectionChildTask extracts its corresponding block as unpartitioned hypergraph and spawns
 * recursively a RecursiveMultilevelBisectionTask for that hypergraph. Once that RecursiveMultilevelBisectionTask is completed, a
 * RecursiveBisectionChildContinuationTask is started and the partition of the recursion is applied to the original hypergraph.
 */

class RecursiveBisectionInitialPartitioner : public IInitialPartitioner {
 private:
  static constexpr bool enable_heavy_assert = false;




 public:
  RecursiveBisectionInitialPartitioner(PartitionedHypergraph& hypergraph,
                                        const Context& context,
                                        const bool top_level,
                                        const TaskGroupID task_group_id);

  RecursiveBisectionInitialPartitioner(const RecursiveBisectionInitialPartitioner&) = delete;
  RecursiveBisectionInitialPartitioner(RecursiveBisectionInitialPartitioner&&) = delete;
  RecursiveBisectionInitialPartitioner & operator= (const RecursiveBisectionInitialPartitioner &) = delete;
  RecursiveBisectionInitialPartitioner & operator= (RecursiveBisectionInitialPartitioner &&) = delete;

 private:
  void initialPartitionImpl() final ;

  PartitionedHypergraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};


}  // namespace mt_kahypar
