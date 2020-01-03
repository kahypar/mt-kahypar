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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

/*!
 * RECURSIVE BISECTION INITIAL PARTITIONER
 * The recursive bisection initial partitioner starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks are reached.
 * Note, the recursive bisection initial partitioner is written in TBB continuation style. The TBB
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
template <typename TypeTraits>
class RecursiveBisectionInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  using BlockRange = std::pair<PartitionID, PartitionID>;

  static constexpr bool debug = false;
  static constexpr bool kahypar_debug = false;
  static constexpr bool enable_heavy_assert = false;

  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  /*!
   * A recursive bisection child task extracts a block of the partition
   * and recursively partition it into the desired number of blocks.
   */
  class RecursiveBisectionChildTask : public tbb::task {

   public:
    RecursiveBisectionChildTask(HyperGraph& hypergraph,
                                const Context& context,
                                const PartitionID block,
                                const BlockRange range,
                                const TaskGroupID task_group_id) :
      _hg(hypergraph),
      _context(context),
      _block(block),
      _range(range),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      const PartitionID k = _range.second - _range.first;
      Context rb_context = setupRecursiveBisectionContext(k);

      // Extracts the block of the hypergraph which we recursively want to partition as
      // seperate unpartitioned hypergraph.
      bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
      auto copy_hypergraph = _hg.copy(k, _task_group_id, _block, cut_net_splitting);
      HyperGraph& rb_hypergraph = copy_hypergraph.first;
      auto& mapping = copy_hypergraph.second;

      // Spawns a new recursive bisection task to partition the current block of the hypergraph
      // into the desired number of blocks
      RecursiveBisectionChildContinuationTask& child_continuation = *new(allocate_continuation())
        RecursiveBisectionChildContinuationTask(_hg, std::move(rb_hypergraph),
          std::move(mapping), std::move(rb_context), _block);
      RecursiveMultilevelBisectionTask& recursion = *new(child_continuation.allocate_child()) RecursiveMultilevelBisectionTask(
        child_continuation.recursiveHypergraph(), child_continuation.recursiveContext(), _task_group_id);
      child_continuation.set_ref_count(1);
      tbb::task::spawn(recursion);

      return nullptr;
    }

   private:
    Context setupRecursiveBisectionContext(const PartitionID k) {
      ASSERT(k >= 2);
      Context rb_context(_context);
      rb_context.partition.k = k;

      rb_context.partition.perfect_balance_part_weights.assign(k, 0);
      rb_context.partition.max_part_weights.assign(k, 0);
      for ( PartitionID part_id = _range.first; part_id < _range.second; ++part_id ) {
        rb_context.partition.perfect_balance_part_weights[part_id - _range.first] =
          _context.partition.perfect_balance_part_weights[part_id];
        rb_context.partition.max_part_weights[part_id - _range.first] =
          _context.partition.max_part_weights[part_id];
      }

      return rb_context;
    }

    HyperGraph& _hg;
    const Context _context;
    const PartitionID _block;
    const BlockRange _range;
    const TaskGroupID _task_group_id;
  };

  /*!
   * Continuation task for the recursive bisection child task. Applies the
   * partition obtained by a recursive bisection task to the original
   * hypergraph.
   */
  class RecursiveBisectionChildContinuationTask : public tbb::task {

   public:
    RecursiveBisectionChildContinuationTask(HyperGraph& original_hypergraph,
                                            Hypergraph&& rb_hypergraph,
                                            parallel::scalable_vector<HypernodeID>&& mapping,
                                            Context&& context,
                                            const PartitionID part_id) :
      _original_hg(original_hypergraph),
      _rb_hg(std::move(rb_hypergraph)),
      _mapping(std::move(mapping)),
      _context(std::move(context)),
      _part_id(part_id) { }

    tbb::task* execute() override {
      // Applying partition of the recursively bisected hypergraph (rb_hg) to
      // original hypergraph (original_hg). All hypernodes that belong to block
      // 'part_id' in the original hypergraph are moved to the block defined in
      // rb_hg.
      ASSERT(_original_hg.initialNumNodes() == _mapping.size());
      for ( const HypernodeID& hn : _original_hg.nodes() ) {
        if ( _original_hg.partID(hn) == _part_id ) {
          const HypernodeID original_id = _original_hg.originalNodeID(hn);
          ASSERT(original_id < _mapping.size());
          PartitionID to = _part_id + _rb_hg.partID(_rb_hg.globalNodeID(_mapping[original_id]));
          ASSERT(to != kInvalidPartition && to < _original_hg.k());
          _original_hg.changeNodePart(hn, _part_id, to);
        }
      }
      return nullptr;
    }

    HyperGraph& recursiveHypergraph() {
      return _rb_hg;
    }

    const Context& recursiveContext() const {
      return _context;
    }

   private:
    HyperGraph& _original_hg;
    Hypergraph _rb_hg;
    const parallel::scalable_vector<HypernodeID> _mapping;
    const Context _context;
    const PartitionID _part_id;
  };

  /*!
   * A bisection task bisects the hypergraph into two block. Internally, it
   * calls our multilevel partitioner for k = 2.
   */
  class MultilevelBisectionTask : public tbb::task {
   public:
    MultilevelBisectionTask(HyperGraph& original_hypergraph,
                  HyperGraph& bisection_hypergraph,
                  parallel::scalable_vector<HypernodeID>& bisection_mapping,
                  const Context& bisection_context,
                  const TaskGroupID task_group_id) :
      _original_hg(original_hypergraph),
      _bisection_hg(bisection_hypergraph),
      _bisection_mapping(bisection_mapping),
      _bisection_context(bisection_context),
      _task_group_id(task_group_id) {}

    tbb::task* execute() override {
      auto copy_hypergraph = _original_hg.copy(2, _task_group_id);
      _bisection_hg = std::move(copy_hypergraph.first);
      _bisection_mapping = std::move(copy_hypergraph.second);

      // Bisect hypergraph with parallel multilevel bisection
      multilevel::partition(_bisection_hg, _bisection_context, false, _task_group_id, this);
      return nullptr;
    }

   private:
    Hypergraph& _original_hg;
    HyperGraph& _bisection_hg;
    parallel::scalable_vector<HypernodeID>& _bisection_mapping;
    const Context& _bisection_context;
    const TaskGroupID _task_group_id;
  };

  /*!
   * Continuation task for the bisection task. Applies the bisection
   * obtained by the bisection task and spawns two recursive child
   * bisection task for each block, that further partitions the hypergraph
   * into the desired number of blocks
   */
  class MultilevelBisectionContinuationTask : public tbb::task {

   public:
    MultilevelBisectionContinuationTask(HyperGraph& hypergraph,
                              const Context& context,
                              const TaskGroupID task_group_id) :
      _bisection_hg(),
      _bisection_mapping(),
      _bisection_context(setupBisectionContext(hypergraph, context)),
      _hg(hypergraph),
      _context(context),
      _task_group_id(task_group_id) {}

    tbb::task* execute() override {
      // Apply partition to hypergraph
      const PartitionID block_0 = 0;
      const PartitionID block_1 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      for (const HypernodeID& hn : _hg.nodes()) {
        const HypernodeID original_id = _hg.originalNodeID(hn);
        ASSERT(original_id < _bisection_mapping.size());
        PartitionID part_id = _bisection_hg.partID(_bisection_hg.globalNodeID(_bisection_mapping[original_id]));
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        if ( part_id == 0 ) {
          _hg.setNodePart(hn, block_0);
        } else {
          _hg.setNodePart(hn, block_1);
        }
      }
      _hg.initializeNumCutHyperedges();
      _hg.updateGlobalPartInfos();

      ASSERT(metrics::objective(_bisection_hg, _context.partition.objective) ==
        metrics::objective(_hg, _context.partition.objective));

      ASSERT(_context.partition.k >= 2);
      PartitionID num_blocks_part_0 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      PartitionID num_blocks_part_1 = _context.partition.k / 2;
      BlockRange range_0 = std::make_pair(0, num_blocks_part_0);
      BlockRange range_1 = std::make_pair(num_blocks_part_0, num_blocks_part_0 + num_blocks_part_1);

      if ( num_blocks_part_0 >= 2 && num_blocks_part_1 >= 2 ) {
        // In case we have to partition both blocks from the bisection further into
        // more than one block, we call the recursive bisection initial partitioner
        // recursively in parallel
        DBG << "Current k = " << _context.partition.k << "\n"
            << "Parallel Recursion 0: k =" << num_blocks_part_0 << "\n"
            << "Parallel Recursion 1: k =" << num_blocks_part_1;

        auto tbb_recursion_task_groups = TBB::instance().create_tbb_task_groups_for_recursion();
        RecursiveMultilevelBisectionContinuationTask& recursive_continuation = *new(allocate_continuation()) RecursiveMultilevelBisectionContinuationTask(_hg);
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, 0, range_0, tbb_recursion_task_groups.first);
        RecursiveBisectionChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, num_blocks_part_0, range_1, tbb_recursion_task_groups.second);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else if ( num_blocks_part_0 >= 2 ) {
        ASSERT(num_blocks_part_1 < 2);
        // In case only the first block has to be partitioned into more than one block, we call
        // the recursive bisection initial partitioner recusively on the block 0
        DBG << "Current k = " << _context.partition.k << ","
            << "Recursion 0: k =" << num_blocks_part_0;
        RecursiveMultilevelBisectionContinuationTask& recursive_continuation = *new(allocate_continuation()) RecursiveMultilevelBisectionContinuationTask(_hg);
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, 0, range_0, _task_group_id);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion_0);
      }

      return nullptr;
    }

    Hypergraph _bisection_hg;
    parallel::scalable_vector<HypernodeID> _bisection_mapping;
    Context _bisection_context;

   private:
    Context setupBisectionContext(const HyperGraph& hypergraph, const Context& context) {
      Context bisection_context(context);

      bisection_context.partition.k = 2;
      bisection_context.partition.verbose_output = debug;
      bisection_context.initial_partitioning.mode = InitialPartitioningMode::direct;
      bisection_context.initial_partitioning.technique = kahypar::InitialPartitioningTechnique::flat;

      // Setup Part Weights
      PartitionID num_blocks_part_0 = context.partition.k / 2 + (context.partition.k % 2 != 0 ? 1 : 0);
      ASSERT(num_blocks_part_0 +  context.partition.k / 2 == context.partition.k);
      bisection_context.partition.perfect_balance_part_weights.assign(2, 0);
      bisection_context.partition.max_part_weights.assign(2, 0);
      for ( PartitionID i = 0; i < num_blocks_part_0; ++i ) {
        bisection_context.partition.perfect_balance_part_weights[0] +=
          context.partition.perfect_balance_part_weights[i];
        bisection_context.partition.max_part_weights[0] +=
          context.partition.max_part_weights[i];
      }
      for ( PartitionID i = num_blocks_part_0; i < context.partition.k; ++i ) {
        bisection_context.partition.perfect_balance_part_weights[1] +=
          context.partition.perfect_balance_part_weights[i];
        bisection_context.partition.max_part_weights[1] +=
          context.partition.max_part_weights[i];
      }

      // Special case, if balance constraint will be violated with this bisection
      HypernodeWeight total_weight = hypergraph.totalWeight();
      HypernodeWeight total_max_part_weight = bisection_context.partition.max_part_weights[0] +
        bisection_context.partition.max_part_weights[1];
      if (total_max_part_weight < total_weight) {
        HypernodeWeight delta = total_weight - total_max_part_weight;
        bisection_context.partition.max_part_weights[0] += std::ceil(((double)delta) / 2.0);
        bisection_context.partition.max_part_weights[1] += std::ceil(((double)delta) / 2.0);
      }

      bisection_context.setupContractionLimit(total_weight);

      return bisection_context;
    }

    HyperGraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };

  /*!
   * The recursive bisection task spawns the multilevel bisection
   * task and its continuation task.
   */
  class RecursiveMultilevelBisectionTask : public tbb::task {

   public:
    RecursiveMultilevelBisectionTask(HyperGraph& hypergraph,
                          const Context& context,
                          const TaskGroupID task_group_id) :
      _hg(hypergraph),
      _context(context),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      ASSERT(_context.partition.k >= 2);
      MultilevelBisectionContinuationTask& bisection_continuation_task = *new(allocate_continuation())
        MultilevelBisectionContinuationTask(_hg, _context, _task_group_id);
      bisection_continuation_task.set_ref_count(1);
      tbb::task::spawn(*new(bisection_continuation_task.allocate_child()) MultilevelBisectionTask(
        _hg, bisection_continuation_task._bisection_hg,
        bisection_continuation_task._bisection_mapping,
        bisection_continuation_task._bisection_context,
        _task_group_id));
      return nullptr;
    }

   private:
    HyperGraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };

  /*!
   * Continuation task for the recursive bisection task. Updates
   * global part sizes and weights.
   */
  class RecursiveMultilevelBisectionContinuationTask : public tbb::task {

   public:
    RecursiveMultilevelBisectionContinuationTask(HyperGraph& hypergraph) :
      _hg(hypergraph) { }

    tbb::task* execute() override {
      _hg.updateGlobalPartInfos();
      return nullptr;
    }

   private:
    HyperGraph& _hg;
  };

 public:
  RecursiveBisectionInitialPartitionerT(HyperGraph& hypergraph,
                                        const Context& context,
                                        const bool top_level,
                                        const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  RecursiveBisectionInitialPartitionerT(const RecursiveBisectionInitialPartitionerT&) = delete;
  RecursiveBisectionInitialPartitionerT(RecursiveBisectionInitialPartitionerT&&) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (const RecursiveBisectionInitialPartitionerT &) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (RecursiveBisectionInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    if (_top_level) {
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    RecursiveMultilevelBisectionTask& root_bisection_task = *new(tbb::task::allocate_root()) RecursiveMultilevelBisectionTask(
      _hg, _context, _task_group_id);
    tbb::task::spawn_root_and_wait(root_bisection_task);

    if (_top_level) {
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

template <typename TypeTraits>
PartitionID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using RecursiveBisectionInitialPartitioner = RecursiveBisectionInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
