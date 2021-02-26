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

#include "mt-kahypar/partition/initial_partitioning/recursive_bisection_initial_partitioner.h"


#include <algorithm>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/multilevel.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/metrics.h"

/** RecursiveBisectionInitialPartitioner Implementation Details
  *
  * Note, the recursive bisection initial partitioner is written in TBBNumaArena continuation style. The TBBNumaArena
  * continuation style is especially useful for recursive patterns. Each task defines its continuation
  * task. A continuation task defines how computation should continue, if all its child tasks are completed.
  * As a consequence, tasks can be spawned without waiting for their completion, because the continuation
  * task is automatically invoked if all child tasks are terminated. Therefore, no thread will waste CPU
  * time while waiting for their recursive tasks to complete.
  *
  * ----------------------
  * The recursive bisection initial partitioner starts by spawning the root RecursiveMultilevelBisectionTask. The RecursiveMultilevelBisectionTask
  * spawns a MultilevelBisectionTask that bisects the hypergraph (multilevel-fashion). Afterwards, the MultilevelBisectionContinuationTask continues
  * and applies the bisection to the hypergraph and spawns two RecursiveBisectionChildTasks. Both are responsible for exactly one block of
  * the partition. The RecursiveBisectionChildTask extracts its corresponding block as unpartitioned hypergraph and spawns
  * recursively a RecursiveMultilevelBisectionTask for that hypergraph. Once that RecursiveMultilevelBisectionTask is completed, a
  * RecursiveBisectionChildContinuationTask is started and the partition of the recursion is applied to the original hypergraph.
*/

namespace mt_kahypar {

  class DoNothingContinuation : public tbb::task {
  public:
    tbb::task* execute() override {
      return nullptr;
    }
  };


  using BlockRange = std::pair<PartitionID, PartitionID>;

  struct OriginalHypergraphInfo {

    double computeAdaptiveEpsilon(const HypernodeWeight current_hypergraph_weight,
                                  const PartitionID current_k) const {
      if ( current_hypergraph_weight == 0 ) {
        return 0.0;
      } else {
        double base = ceil(static_cast<double>(original_hypergraph_weight) / original_k)
                      / ceil(static_cast<double>(current_hypergraph_weight) / current_k)
                      * (1.0 + original_epsilon);
        double adaptive_epsilon = std::min(0.99, std::max(std::pow(base, 1.0 /
                                                                         ceil(log2(static_cast<double>(current_k)))) - 1.0,0.0));
        return adaptive_epsilon;
      }
    }

    const HypernodeWeight original_hypergraph_weight;
    const PartitionID original_k;
    const double original_epsilon;
  };


  /*!
   * Continuation task for the recursive bisection child task. Applies the
   * partition obtained by a recursive bisection task to the original
   * hypergraph.
   */
  class RecursiveBisectionChildContinuationTask : public tbb::task {

    using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
#define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  public:
    RecursiveBisectionChildContinuationTask(PartitionedHypergraph& original_hypergraph,
                                            Context&& context,
                                            Hypergraph&& rb_hypergraph,
                                            parallel::scalable_vector<HypernodeID>&& mapping,
                                            const PartitionID k,
                                            const TaskGroupID task_group_id,
                                            const PartitionID part_id) :
            _original_hg(original_hypergraph),
            _context(std::move(context)),
            _rb_hg(std::move(rb_hypergraph)),
            _rb_partitioned_hg(k, task_group_id, _rb_hg),
            _mapping(std::move(mapping)),
            _task_group_id(task_group_id),
            _part_id(part_id) { }

    tbb::task* execute() override {
      // Applying partition of the recursively bisected hypergraph (_rb_partitioned_hg) to
      // original hypergraph (original_hg). All hypernodes that belong to block
      // 'part_id' in the original hypergraph are moved to the block defined in
      // _rb_partitioned_hg.
      ASSERT(_original_hg.initialNumNodes() == _mapping.size());
      _original_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        if ( _original_hg.partID(hn) == _part_id ) {
          ASSERT(hn < _mapping.size());
          PartitionID to = _part_id + _rb_partitioned_hg.partID(_mapping[hn]);
          ASSERT(to != kInvalidPartition && to < _original_hg.k());
          if ( _part_id != to ) {
            _original_hg.changeNodePart(hn, _part_id, to, NOOP_FUNC);
          }
        }
      });
      return nullptr;
    }

    PartitionedHypergraph& recursivePartitionedHypergraph() {
      return _rb_partitioned_hg;
    }

    const Context& recursiveContext() const {
      return _context;
    }

  private:
    PartitionedHypergraph& _original_hg;
    const Context _context;
    Hypergraph _rb_hg;
    PartitionedHypergraph _rb_partitioned_hg;
    const parallel::scalable_vector<HypernodeID> _mapping;
    const TaskGroupID _task_group_id;
    const PartitionID _part_id;
  };

  /*!
   * A bisection task bisects the hypergraph into two block. Internally, it
   * calls our multilevel partitioner for k = 2.
   */
  class MultilevelBisectionTask : public tbb::task {
  public:
    MultilevelBisectionTask(PartitionedHypergraph& original_hypergraph,
                            Hypergraph& bisection_hypergraph,
                            PartitionedHypergraph& bisection_partitioned_hypergraph,
                            const Context& bisection_context,
                            const TaskGroupID task_group_id) :
            _original_hg(original_hypergraph),
            _bisection_hg(bisection_hypergraph),
            _bisection_partitioned_hg(bisection_partitioned_hypergraph),
            _bisection_context(bisection_context),
            _task_group_id(task_group_id) {}

    tbb::task* execute() override {
      // Bisect hypergraph with parallel multilevel bisection
      _bisection_hg = _original_hg.hypergraph().copy(_task_group_id);
      multilevel::partition_async(_bisection_hg, _bisection_partitioned_hg,
                                  _bisection_context, false, _task_group_id, this);
      return nullptr;
    }

  private:
    PartitionedHypergraph& _original_hg;
    Hypergraph& _bisection_hg;
    PartitionedHypergraph& _bisection_partitioned_hg;
    const Context& _bisection_context;
    const TaskGroupID _task_group_id;
  };

  /*!
 * A recursive bisection child task extracts a block of the partition
 * and recursively partition it into the desired number of blocks.
 */
  class RecursiveBisectionChildTask : public tbb::task {

  public:
    RecursiveBisectionChildTask(const OriginalHypergraphInfo original_hypergraph_info,
                                PartitionedHypergraph& hypergraph,
                                const Context& context,
                                const PartitionID block,
                                const BlockRange range,
                                const TaskGroupID task_group_id,
                                const double degree_of_parallelism) :
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _block(block),
            _range(range),
            _task_group_id(task_group_id),
            _degree_of_parallelism(degree_of_parallelism) { }

    tbb::task* execute() override ;

  private:
    Context setupRecursiveBisectionContext(const PartitionID k) {
      ASSERT(k >= 2);
      Context rb_context(_context);
      rb_context.partition.k = k;
      rb_context.type = kahypar::ContextType::initial_partitioning;

      rb_context.partition.perfect_balance_part_weights.assign(k, 0);
      rb_context.partition.max_part_weights.assign(k, 0);
      for ( PartitionID part_id = _range.first; part_id < _range.second; ++part_id ) {
        rb_context.partition.perfect_balance_part_weights[part_id - _range.first] =
                _context.partition.perfect_balance_part_weights[part_id];
        rb_context.partition.max_part_weights[part_id - _range.first] =
                _context.partition.max_part_weights[part_id];
      }

      rb_context.shared_memory.degree_of_parallelism *= _degree_of_parallelism;

      return rb_context;
    }

    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context _context;
    const PartitionID _block;
    const BlockRange _range;
    const TaskGroupID _task_group_id;
    const double _degree_of_parallelism;
  };

  /*!
 * Continuation task for the bisection task. Applies the bisection
 * obtained by the bisection task and spawns two recursive child
 * bisection task for each block, that further partitions the hypergraph
 * into the desired number of blocks
 */
  class MultilevelBisectionContinuationTask : public tbb::task {
    static constexpr bool debug = false;
  public:
    MultilevelBisectionContinuationTask(const OriginalHypergraphInfo original_hypergraph_info,
                                        PartitionedHypergraph& hypergraph,
                                        const Context& context,
                                        const TaskGroupID task_group_id) :
            _bisection_hg(),
            _bisection_partitioned_hg(),
            _bisection_context(),
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _task_group_id(task_group_id) {
      _bisection_context = setupBisectionContext(_hg.hypergraph(), _context);
    }

    tbb::task* execute() override {
      ASSERT(_hg.initialNumNodes() == _bisection_hg.initialNumNodes());
      // Apply partition to hypergraph
      const PartitionID block_0 = 0;
      const PartitionID block_1 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        PartitionID part_id = _bisection_partitioned_hg.partID(hn);
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        if ( part_id == 0 ) {
          _hg.setOnlyNodePart(hn, block_0);
        } else {
          _hg.setOnlyNodePart(hn, block_1);
        }
      });
      _hg.initializePartition(_task_group_id);

      ASSERT(metrics::objective(_bisection_partitioned_hg, _context.partition.objective) ==
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

        auto tbb_recursion_task_groups = TBBNumaArena::instance().create_tbb_task_groups_for_recursion();
        DoNothingContinuation& recursive_continuation = *new(allocate_continuation()) DoNothingContinuation();
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
                _original_hypergraph_info, _hg, _context, 0, range_0, tbb_recursion_task_groups.first, 0.5);
        RecursiveBisectionChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
                _original_hypergraph_info, _hg, _context, num_blocks_part_0, range_1, tbb_recursion_task_groups.second, 0.5);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else if ( num_blocks_part_0 >= 2 ) {
        ASSERT(num_blocks_part_1 < 2);
        // In case only the first block has to be partitioned into more than one block, we call
        // the recursive bisection initial partitioner recusively on the block 0
        DBG << "Current k = " << _context.partition.k << ","
            << "Recursion 0: k =" << num_blocks_part_0;
        DoNothingContinuation& recursive_continuation = *new(allocate_continuation()) DoNothingContinuation();
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
                _original_hypergraph_info, _hg, _context, 0, range_0, _task_group_id, 1.0);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion_0);
      }

      return nullptr;
    }

    Hypergraph _bisection_hg;
    PartitionedHypergraph _bisection_partitioned_hg;
    Context _bisection_context;

  private:
    Context setupBisectionContext(const Hypergraph& hypergraph, const Context& context) {
      Context bisection_context(context);

      bisection_context.partition.k = 2;
      bisection_context.partition.verbose_output = debug;
      bisection_context.initial_partitioning.mode = InitialPartitioningMode::direct;
      bisection_context.type = kahypar::ContextType::initial_partitioning;

      // Setup Part Weights
      const HypernodeWeight total_weight = hypergraph.totalWeight();
      const PartitionID k = context.partition.k;
      const PartitionID k0 = k / 2 + (k % 2 != 0 ? 1 : 0);
      const PartitionID k1 = k / 2;
      ASSERT(k0 + k1 == context.partition.k);
      if ( context.partition.use_individual_part_weights ) {
        const HypernodeWeight max_part_weights_sum = std::accumulate(context.partition.max_part_weights.cbegin(),
                                                                    context.partition.max_part_weights.cend(), 0);
        const double weight_fraction = total_weight / static_cast<double>(max_part_weights_sum);
        ASSERT(weight_fraction <= 1.0);
        bisection_context.partition.perfect_balance_part_weights.clear();
        bisection_context.partition.max_part_weights.clear();
        HypernodeWeight perfect_weight_p0 = 0;
        for ( PartitionID i = 0; i < k0; ++i ) {
          perfect_weight_p0 += ceil(weight_fraction * context.partition.max_part_weights[i]);
        }
        HypernodeWeight perfect_weight_p1 = 0;
        for ( PartitionID i = k0; i < k; ++i ) {
          perfect_weight_p1 += ceil(weight_fraction * context.partition.max_part_weights[i]);
        }
        // In the case of individual part weights, the usual adaptive epsilon formula is not applicable because it
        // assumes equal part weights. However, by observing that ceil(current_weight / current_k) is the current
        // perfect part weight and (1 + epsilon)ceil(original_weight / original_k) is the maximum part weight,
        // we can derive an equivalent formula using the sum of the perfect part weights and the sum of the
        // maximum part weights.
        // Note that the sum of the perfect part weights might be unequal to the hypergraph weight due to rounding.
        // Thus, we need to use the former instead of using the hypergraph weight directly, as otherwise it could
        // happen that (1 + epsilon)perfect_part_weight > max_part_weight because of rounding issues.
        const double base = max_part_weights_sum / static_cast<double>(perfect_weight_p0 + perfect_weight_p1);
        bisection_context.partition.epsilon = total_weight == 0 ? 0 : std::min(0.99, std::max(std::pow(base, 1.0 /
                                                                      ceil(log2(static_cast<double>(k)))) - 1.0,0.0));
        bisection_context.partition.perfect_balance_part_weights.push_back(perfect_weight_p0);
        bisection_context.partition.perfect_balance_part_weights.push_back(perfect_weight_p1);
        bisection_context.partition.max_part_weights.push_back(
                round((1 + bisection_context.partition.epsilon) * perfect_weight_p0));
        bisection_context.partition.max_part_weights.push_back(
                round((1 + bisection_context.partition.epsilon) * perfect_weight_p1));
      } else {
        bisection_context.partition.epsilon = _original_hypergraph_info.computeAdaptiveEpsilon(total_weight, k);

        bisection_context.partition.perfect_balance_part_weights.clear();
        bisection_context.partition.max_part_weights.clear();
        bisection_context.partition.perfect_balance_part_weights.push_back(
                std::ceil(k0 / static_cast<double>(k) * static_cast<double>(total_weight)));
        bisection_context.partition.perfect_balance_part_weights.push_back(
                std::ceil(k1 / static_cast<double>(k) * static_cast<double>(total_weight)));
        bisection_context.partition.max_part_weights.push_back(
                (1 + bisection_context.partition.epsilon) * bisection_context.partition.perfect_balance_part_weights[0]);
        bisection_context.partition.max_part_weights.push_back(
                (1 + bisection_context.partition.epsilon) * bisection_context.partition.perfect_balance_part_weights[1]);
      }
      bisection_context.setupContractionLimit(total_weight);
      bisection_context.setupSparsificationParameters();

      return bisection_context;
    }

    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };


  /*!
  * The recursive bisection task spawns the multilevel bisection
  * task and its continuation task.
  */
  class RecursiveMultilevelBisectionTask : public tbb::task {

  public:
    RecursiveMultilevelBisectionTask(const OriginalHypergraphInfo original_hypergraph_info,
                                     PartitionedHypergraph& hypergraph,
                                     const Context& context,
                                     const TaskGroupID task_group_id) :
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      ASSERT(_context.partition.k >= 2);
      MultilevelBisectionContinuationTask& bisection_continuation_task = *new(allocate_continuation())
              MultilevelBisectionContinuationTask(_original_hypergraph_info, _hg, _context, _task_group_id);
      bisection_continuation_task.set_ref_count(1);
      tbb::task::spawn(*new(bisection_continuation_task.allocate_child())
              MultilevelBisectionTask(
              _hg, bisection_continuation_task._bisection_hg,
              bisection_continuation_task._bisection_partitioned_hg,
              bisection_continuation_task._bisection_context,
              _task_group_id));
      return nullptr;
    }

  private:
    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };


  tbb::task * RecursiveBisectionChildTask::execute() {
    const PartitionID k = _range.second - _range.first;
    Context rb_context = setupRecursiveBisectionContext(k);

    // Extracts the block of the hypergraph which we recursively want to partition as
    // seperate unpartitioned hypergraph.
    bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
    auto copy_hypergraph = _hg.extract(_task_group_id, _block, cut_net_splitting);
    Hypergraph& rb_hypergraph = copy_hypergraph.first;
    auto& mapping = copy_hypergraph.second;

    // Spawns a new recursive bisection task to partition the current block of the hypergraph
    // into the desired number of blocks
    if ( rb_hypergraph.initialNumNodes() > 0 ) {
      RecursiveBisectionChildContinuationTask& child_continuation = *new(allocate_continuation())
              RecursiveBisectionChildContinuationTask(_hg, std::move(rb_context),
                                                      std::move(rb_hypergraph), std::move(mapping), k, _task_group_id, _block);
      RecursiveMultilevelBisectionTask& recursion = *new(child_continuation.allocate_child())
              RecursiveMultilevelBisectionTask(
              _original_hypergraph_info,
              child_continuation.recursivePartitionedHypergraph(),
              child_continuation.recursiveContext(), _task_group_id);
      child_continuation.set_ref_count(1);
      tbb::task::spawn(recursion);
    }

    return nullptr;
  }

  RecursiveBisectionInitialPartitioner::RecursiveBisectionInitialPartitioner(PartitionedHypergraph& hypergraph,
                                                                             const Context& context,
                                                                             const bool top_level,
                                                                             const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  void RecursiveBisectionInitialPartitioner::initialPartitionImpl() {
    if (_top_level) {
      parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    RecursiveMultilevelBisectionTask& root_bisection_task = *new(tbb::task::allocate_root()) RecursiveMultilevelBisectionTask(
            OriginalHypergraphInfo { _hg.totalWeight(), _context.partition.k, _context.partition.epsilon },
            _hg, _context, _task_group_id);
    tbb::task::spawn_root_and_wait(root_bisection_task);

    if (_top_level) {
      parallel::MemoryPool::instance().activate_unused_memory_allocations();
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }
} // namepace mt_kahypar