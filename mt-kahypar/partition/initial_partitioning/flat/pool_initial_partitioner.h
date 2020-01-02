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

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {

template<typename TypeTraits>
class PoolInitialPartitionerContinuationT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;

  public:
  PoolInitialPartitionerContinuationT(HyperGraph& hypergraph,
                                      const Context& context,
                                      const TaskGroupID task_group_id) :
    _ip_data(hypergraph, context, task_group_id) { }

  tbb::task* execute() override {
    _ip_data.apply();
    return nullptr;
  }

  InitialPartitioningDataContainer _ip_data;
};

template<typename TypeTraits>
static void spawn_initial_partitioner(PoolInitialPartitionerContinuationT<TypeTraits>& continuation_task,
                                      const Context& context ) {
  // For each initial partitioner, we create exactly number of initial partitioning runs
  // tasks => each tasks creates exactly one partition.
  tbb::task_list ip_tasks;
  for ( uint8_t i = 0; i < static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED); ++i ) {
    InitialPartitioningAlgorithm algorithm = static_cast<InitialPartitioningAlgorithm>(i);
    for ( size_t j = 0; j < context.initial_partitioning.runs; ++j ) {
      std::unique_ptr<tbb::task> initial_partitioner_ptr =
        FlatInitialPartitionerFactory::getInstance().createObject(
          algorithm, &continuation_task, continuation_task._ip_data, context);
      tbb::task* initial_partitioner = initial_partitioner_ptr.release();
      ip_tasks.push_back(*initial_partitioner);
    }
  }

  // Spawn Initial Partitioner
  continuation_task.set_ref_count(context.initial_partitioning.runs *
    static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED));
  tbb::task::spawn(ip_tasks);
}

/*!
 * The pool initial partitioner spawns for each initial partitioning run and algorithm
 * exactly one initial partitioning task. The number of initial partitions computed during
 * an invocation of the pool initial partitioner is exactly ( num IP runs ) * (num IP algos).
 * The best partition is applied to the hypergraph.
 */
template<typename TypeTraits>
class PoolInitialPartitionerT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PoolInitialPartitionerContinuation = PoolInitialPartitionerContinuationT<TypeTraits>;

  static constexpr bool debug = false;

 public:
  PoolInitialPartitionerT(HyperGraph& hypergraph,
                         const Context& context,
                         const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id) { }

  tbb::task* execute() override {
    PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
      PoolInitialPartitionerContinuation(_hg, _context, _task_group_id);
    spawn_initial_partitioner(ip_continuation, _context);
    return nullptr;
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
};
} // namespace mt_kahypar
