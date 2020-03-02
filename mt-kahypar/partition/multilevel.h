/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <memory>

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"

namespace mt_kahypar {
namespace multilevel {

namespace {

using PartitionedHypergraph = mt_kahypar::PartitionedHypergraph;

class RefinementTask : public tbb::task {

 public:
  RefinementTask(Hypergraph& hypergraph,
                 PartitionedHypergraph& partitioned_hypergraph,
                 const Context& context,
                 const bool top_level,
                 const TaskGroupID task_group_id) :
    _coarsener(nullptr),
    _hg(hypergraph),
    _partitioned_hg(partitioned_hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) {
    // Must be empty, because final partitioned hypergraph
    // is moved into this object
    _coarsener = CoarsenerFactory::getInstance().createObject(
      _context.coarsening.algorithm, _hg, _context, _task_group_id);
  }

  tbb::task* execute() override {
    enableTimerAndStats();
    utils::Timer::instance().stop_timer("initial_partitioning");

    PartitionedHypergraph& coarsest_partitioned_hypergraph =
      _coarsener->coarsestPartitionedHypergraph();
    io::printPartitioningResults(coarsest_partitioned_hypergraph,
      _context, "Initial Partitioning Results:");
    if ( _context.partition.verbose_output ) {
      utils::InitialPartitioningStats::instance().printInitialPartitioningStats();
    }

    // ################## LOCAL SEARCH ##################
    io::printLocalSearchBanner(_context);
    utils::Timer::instance().start_timer("refinement", "Refinement");
    std::unique_ptr<IRefiner> label_propagation =
      LabelPropagationFactory::getInstance().createObject(
        _context.refinement.label_propagation.algorithm,
        coarsest_partitioned_hypergraph, _context, _task_group_id);

    _partitioned_hg = _coarsener->uncoarsen(label_propagation);
    utils::Timer::instance().stop_timer("refinement");

    io::printPartitioningResults(_partitioned_hg, _context, "Local Search Results:");
    return nullptr;
  }

 public:
  std::unique_ptr<ICoarsener> _coarsener;

 private:
  void enableTimerAndStats() {
    if ( _top_level ) {
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

  Hypergraph& _hg;
  PartitionedHypergraph& _partitioned_hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

class CoarseningTask : public tbb::task {
  using PoolInitialPartitionerContinuation = PoolInitialPartitionerContinuationT<GlobalTypeTraits>;

 public:
  CoarseningTask(Hypergraph& hypergraph,
                 const Context& context,
                 ICoarsener& coarsener,
                 const bool top_level,
                 const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _coarsener(coarsener),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  tbb::task* execute() override {
    // ################## COARSENING ##################
    mt_kahypar::io::printCoarseningBanner(_context);
    utils::Timer::instance().start_timer("coarsening", "Coarsening");
    _coarsener.coarsen();
    utils::Timer::instance().stop_timer("coarsening");

    if (_context.partition.verbose_output) {
      mt_kahypar::io::printHypergraphInfo(
        _coarsener.coarsestHypergraph(), "Coarsened Hypergraph");
    }

    // ################## INITIAL PARTITIONING ##################
    PartitionedHypergraph& coarsest_partitioned_hypergraph =
      _coarsener.coarsestPartitionedHypergraph();

    io::printInitialPartitioningBanner(_context);
    utils::Timer::instance().start_timer("initial_partitioning", "Initial Partitioning");
    if ( _context.initial_partitioning.mode == InitialPartitioningMode::direct ) {
      disableTimerAndStats();
      PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
        PoolInitialPartitionerContinuation(
          coarsest_partitioned_hypergraph, _context, _task_group_id);
      spawn_initial_partitioner(ip_continuation);
    } else {
      std::unique_ptr<IInitialPartitioner> initial_partitioner =
        InitialPartitionerFactory::getInstance().createObject(
          _context.initial_partitioning.mode, coarsest_partitioned_hypergraph,
          _context, _top_level, _task_group_id);
      initial_partitioner->initialPartition();
    }
    return nullptr;
  }

 private:
  void disableTimerAndStats() {
    if ( _top_level ) {
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }
  }

  Hypergraph& _hg;
  const Context& _context;
  ICoarsener& _coarsener;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

// ! Helper function that spawns the multilevel partitioner in
// ! TBB continuation style with a given parent task.
static void spawn_multilevel_partitioner(Hypergraph& hypergraph,
                                         PartitionedHypergraph& partitioned_hypergraph,
                                         const Context& context,
                                         const bool top_level,
                                         const TaskGroupID task_group_id,
                                         tbb::task& parent) {
  // The coarsening task is first executed and once it finishes the
  // refinement task continues (without blocking)
  RefinementTask& refinement_task = *new(parent.allocate_continuation())
    RefinementTask(hypergraph, partitioned_hypergraph, context, top_level, task_group_id);
  refinement_task.set_ref_count(1);
  CoarseningTask& coarsening_task = *new(refinement_task.allocate_child()) CoarseningTask(
    hypergraph, context, *refinement_task._coarsener, top_level, task_group_id);
  tbb::task::spawn(coarsening_task);
}

class MultilevelPartitioningTask : public tbb::task {

 public:
  MultilevelPartitioningTask(Hypergraph& hypergraph,
                             PartitionedHypergraph& partitioned_hypergraph,
                             const Context& context,
                             const bool top_level,
                             const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _partitioned_hg(partitioned_hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  tbb::task* execute() override {
    spawn_multilevel_partitioner(
    _hg, _partitioned_hg, _context, _top_level, _task_group_id, *this);
    return nullptr;
  }

 private:
  Hypergraph& _hg;
  PartitionedHypergraph& _partitioned_hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

} // namespace

// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB blocking-style.
static inline PartitionedHypergraph partition(Hypergraph& hypergraph,
                                              const Context& context,
                                              const bool top_level,
                                              const TaskGroupID task_group_id) {
  PartitionedHypergraph partitioned_hypergraph;
  MultilevelPartitioningTask& multilevel_task = *new(tbb::task::allocate_root())
    MultilevelPartitioningTask(
      hypergraph, partitioned_hypergraph,
      context, top_level, task_group_id);
  tbb::task::spawn_root_and_wait(multilevel_task);
  return partitioned_hypergraph;
}

// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB continuation-style.
// ! Note, the final partitioned hypergraph is moved into the
// ! passed partitioned hypergraph object.
static inline void partition_async(Hypergraph& hypergraph,
                                   PartitionedHypergraph& partitioned_hypergraph,
                                   const Context& context,
                                   const bool top_level,
                                   const TaskGroupID task_group_id,
                                   tbb::task* parent) {
  ASSERT(parent);
  spawn_multilevel_partitioner(
    hypergraph, partitioned_hypergraph, context,
    top_level, task_group_id, *parent);
}

}  // namespace multilevel
}  // namespace mt_kahypar
