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

#include "pool_initial_partitioner.h"

#include "mt-kahypar/partition/registries/register_flat_initial_partitioning_algorithms.h"

namespace mt_kahypar {

  class DoNothingContinuation : public tbb::task {
  public:
    tbb::task* execute() override {
      return nullptr;
    }
  };


  class SpawnInitialPartitionerTaskList : public tbb::task {

  public:
    SpawnInitialPartitionerTaskList(InitialPartitioningDataContainer& ip_data,
                                    const Context& context,
                                    IPTaskList ip_tasks) :
            _ip_data(ip_data),
            _context(context),
            _ip_tasks(ip_tasks) { }

    tbb::task* execute() override {
      DoNothingContinuation& task_continuation = *new(allocate_continuation()) DoNothingContinuation();
      task_continuation.set_ref_count(_ip_tasks.size());
      // Spawn Initial Partitioner Tasks
      for ( const auto& [algorithm, seed] : _ip_tasks ) {
        std::unique_ptr<tbb::task> initial_partitioner_ptr =
                FlatInitialPartitionerFactory::getInstance().createObject(
                        algorithm, &task_continuation, algorithm, _ip_data, _context, seed);
        tbb::task* initial_partitioner = initial_partitioner_ptr.release();
        tbb::task::spawn(*initial_partitioner);
      }
      return nullptr;
    }

  private:
    InitialPartitioningDataContainer& _ip_data;
    const Context& _context;
    IPTaskList _ip_tasks;
  };


  PoolInitialPartitionerContinuation::PoolInitialPartitionerContinuation(PartitionedHypergraph& hypergraph,
                                     const Context& context,
                                     const TaskGroupID task_group_id) :
          _ip_data(hypergraph, context, task_group_id),
          _context(context),
          _ip_task_lists(context.shared_memory.num_threads) {

    ASSERT(context.shared_memory.num_threads > 0);
    if ( context.initial_partitioning.enabled_ip_algos.size() <
         static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED) ) {
      ERROR("Size of enabled IP algorithms vector is smaller than number of IP algorithms!");
    }
    // Initial Partitioner tasks are evenly distributed among different task list. For an
    // explanation why we do this, see spawn_initial_partitioner(...)
    size_t task_list_idx = 0;
    std::mt19937 rng(context.partition.seed);
    for ( uint8_t i = 0; i < static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED); ++i ) {
      if ( context.initial_partitioning.enabled_ip_algos[i] ) {
        auto algorithm = static_cast<InitialPartitioningAlgorithm>(i);
        for ( size_t j = 0; j < _context.initial_partitioning.runs; ++j ) {
          _ip_task_lists[task_list_idx].emplace_back(algorithm, rng());
          task_list_idx = (task_list_idx + 1) % _context.shared_memory.num_threads;
        }
      }
    }
  }

  tbb::task* PoolInitialPartitionerContinuation::execute() {
    _ip_data.apply();
    return nullptr;
  }

  void spawn_initial_partitioner(PoolInitialPartitionerContinuation& continuation_task ) {
    // Spawn Initial Partitioner
    const Context& context = continuation_task._context;
    continuation_task.set_ref_count(context.shared_memory.num_threads);
    for ( size_t i = 0; i < context.shared_memory.num_threads; ++i ) {
      // Note, we first spawn exactly num threads tasks that spawns a subset
      // of the initial partitioner tasks. Alternatively, we could also spawn
      // all initial partitioner tasks directly here, but this would insert
      // all tasks into the tbb task queue of one thread from which the other
      // threads have to steal from. This can become a major sequential bottleneck.
      // Therefore, we introduce that indirection such that the initial partitioner
      // tasks are more evenly distributed among the tbb task queues of all threads.
      tbb::task::spawn(*new(continuation_task.allocate_child())
              SpawnInitialPartitionerTaskList(
              continuation_task._ip_data, context, continuation_task._ip_task_lists[i]));
    }
  }

  PoolInitialPartitioner::PoolInitialPartitioner(PartitionedHypergraph& hypergraph,
                         const Context& context,
                         const TaskGroupID task_group_id) :
          _hg(hypergraph),
          _context(context),
          _task_group_id(task_group_id) { }

  tbb::task* PoolInitialPartitioner::execute() {
    PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
            PoolInitialPartitionerContinuation(_hg, _context, _task_group_id);
    spawn_initial_partitioner(ip_continuation);
    return nullptr;
  }

} // namespace mt_kahypar