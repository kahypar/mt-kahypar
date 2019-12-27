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

#include <libkahypar.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

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
      assignPartitionFromRecursionToOriginalHypergraph();
      return nullptr;
    }

    HyperGraph& recursiveHypergraph() {
      return _rb_hg;
    }

    const Context& recursiveContext() const {
      return _context;
    }

   private:
    void assignPartitionFromRecursionToOriginalHypergraph() {
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
    }

    HyperGraph& _original_hg;
    Hypergraph _rb_hg;
    const parallel::scalable_vector<HypernodeID> _mapping;
    const Context _context;
    const PartitionID _part_id;
  };

  class RecursiveBisectionChildTask : public tbb::task {

   public:
    RecursiveBisectionChildTask(HyperGraph& hypergraph,
                                const Context& context,
                                const PartitionID block,
                                const BlockRange range,
                                const bool top_level,
                                const TaskGroupID task_group_id) :
      _hg(hypergraph),
      _context(context),
      _block(block),
      _range(range),
      _top_level(top_level),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      const PartitionID k = _range.second - _range.first;
      Context rb_context = setupRecursiveBisectionContext(k);

      // Extracts the block, which we want to recursively partition, as new hypergraph
      // and calls the recursive bisection initial partitioner recursively.
      bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
      auto copy_hypergraph = _hg.copy(k, _task_group_id, _block, cut_net_splitting);
      HyperGraph& rb_hypergraph = copy_hypergraph.first;
      auto& mapping = copy_hypergraph.second;

      RecursiveBisectionChildContinuationTask& child_continuation = *new(allocate_continuation())
        RecursiveBisectionChildContinuationTask(_hg, std::move(rb_hypergraph),
          std::move(mapping), std::move(rb_context), _block);
      RecursiveBisectionTask& recursion = *new(child_continuation.allocate_child()) RecursiveBisectionTask(
        child_continuation.recursiveHypergraph(), child_continuation.recursiveContext(), false, _task_group_id);
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
    const bool _top_level;
    const TaskGroupID _task_group_id;
  };

  class RecursiveBisectionContinuationTask : public tbb::task {

   public:
    RecursiveBisectionContinuationTask(HyperGraph& hypergraph) :
      _hg(hypergraph) { }

    tbb::task* execute() override {
      _hg.updateGlobalPartInfos();
      return nullptr;
    }

   private:
    HyperGraph& _hg;
  };

  class RecursiveBisectionTask : public tbb::task {

   public:
    RecursiveBisectionTask(HyperGraph& hypergraph,
                          const Context& context,
                          const bool top_level,
                          const TaskGroupID task_group_id) :
      _hg(hypergraph),
      _context(context),
      _top_level(top_level),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      ASSERT(_context.partition.k >= 2);

      PartitionID num_blocks_part_0 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      PartitionID num_blocks_part_1 = _context.partition.k / 2;
      BlockRange range_0 = std::make_pair(0, num_blocks_part_0);
      BlockRange range_1 = std::make_pair(num_blocks_part_0, num_blocks_part_0 + num_blocks_part_1);

      // Bisecting the hypergraph into two blocks
      utils::Timer::instance().start_timer("top_level_bisection", "Top Level Bisection", false, _top_level);
      bisect(0, num_blocks_part_0);
      utils::Timer::instance().stop_timer("top_level_bisection", _top_level);

      if ( num_blocks_part_0 >= 2 && num_blocks_part_1 >= 2 ) {
        // In case we have to partition both blocks from the bisection further into
        // more than one block, we call the recursive bisection initial partitioner
        // recursively in parallel
        DBG << "Current k = " << _context.partition.k << "\n"
            << "Parallel Recursion 0: k =" << num_blocks_part_0 << "\n"
            << "Parallel Recursion 1: k =" << num_blocks_part_1;

        auto tbb_recursion_task_groups = TBB::instance().create_tbb_task_groups_for_recursion();
        RecursiveBisectionContinuationTask& recursive_continuation = *new(allocate_continuation()) RecursiveBisectionContinuationTask(_hg);
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, 0, range_0, _top_level, tbb_recursion_task_groups.first);
        RecursiveBisectionChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, num_blocks_part_0, range_1, _top_level, tbb_recursion_task_groups.second);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else if ( num_blocks_part_0 >= 2 ) {
        ASSERT(num_blocks_part_1 < 2);
        // In case only the first block has to be partitioned into more than one block, we call
        // the recursive bisection initial partitioner recusively on the block 0
        DBG << "Current k = " << _context.partition.k << ","
            << "Sequential Recursion 0: k =" << num_blocks_part_0;
        RecursiveBisectionContinuationTask& recursive_continuation = *new(allocate_continuation()) RecursiveBisectionContinuationTask(_hg);
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _hg, _context, 0, range_0, _top_level, _task_group_id);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion_0);
      }

      return nullptr;
    }

   private:

    void bisect(const PartitionID block_0, const PartitionID block_1) {
      ASSERT(block_0 < _context.partition.k);
      ASSERT(block_1 < _context.partition.k);
      Context bisection_context = setupBisectionContext();

      utils::Timer::instance().start_timer("top_level_copy", "Top Level Copy", false, _top_level);
      auto copy_hypergraph = _hg.copy(2, _task_group_id);
      HyperGraph& tmp_hg = copy_hypergraph.first;
      auto& mapping = copy_hypergraph.second;
      utils::Timer::instance().stop_timer("top_level_copy", _top_level);

      utils::Timer::instance().start_timer("top_level_parallel_bisection", "Top Level Parallel Bisection", false, _top_level);
      // Bisect hypergraph with parallel multilevel bisection
      multilevel::partition(tmp_hg, bisection_context, false, _task_group_id);
      utils::Timer::instance().stop_timer("top_level_parallel_bisection", _top_level);

      // Apply partition to hypergraph
      utils::Timer::instance().start_timer("top_level_apply_bisection", "Top Level Apply Bisection", false, _top_level);
      for (const HypernodeID& hn : _hg.nodes()) {
        const HypernodeID original_id = _hg.originalNodeID(hn);
        ASSERT(original_id < mapping.size());
        PartitionID part_id = tmp_hg.partID(tmp_hg.globalNodeID(mapping[original_id]));
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        if ( part_id == 0 ) {
          _hg.setNodePart(hn, block_0);
        } else {
          _hg.setNodePart(hn, block_1);
        }
      }
      _hg.initializeNumCutHyperedges();
      _hg.updateGlobalPartInfos();

      utils::Timer::instance().stop_timer("top_level_apply_bisection", _top_level);

      ASSERT(metrics::objective(tmp_hg, _context.partition.objective) ==
        metrics::objective(_hg, _context.partition.objective));
    }

    Context setupBisectionContext() {
      Context bisection_context(_context);

      bisection_context.partition.k = 2;
      bisection_context.partition.verbose_output = debug;
      bisection_context.initial_partitioning.mode = InitialPartitioningMode::direct;
      bisection_context.initial_partitioning.technique = kahypar::InitialPartitioningTechnique::flat;

      // Setup Part Weights
      PartitionID num_blocks_part_0 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      ASSERT(num_blocks_part_0 +  _context.partition.k / 2 == _context.partition.k);
      bisection_context.partition.perfect_balance_part_weights.assign(2, 0);
      bisection_context.partition.max_part_weights.assign(2, 0);
      for ( PartitionID i = 0; i < num_blocks_part_0; ++i ) {
        bisection_context.partition.perfect_balance_part_weights[0] +=
          _context.partition.perfect_balance_part_weights[i];
        bisection_context.partition.max_part_weights[0] +=
          _context.partition.max_part_weights[i];
      }
      for ( PartitionID i = num_blocks_part_0; i < _context.partition.k; ++i ) {
        bisection_context.partition.perfect_balance_part_weights[1] +=
          _context.partition.perfect_balance_part_weights[i];
        bisection_context.partition.max_part_weights[1] +=
          _context.partition.max_part_weights[i];
      }

      // Special case, if balance constraint will be violated with this bisection
      // => causes KaHyPar to exit with failure
      HypernodeWeight total_weight = _hg.totalWeight();
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
    const bool _top_level;
    const TaskGroupID _task_group_id;
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

    RecursiveBisectionTask& root_bisection_task = *new(tbb::task::allocate_root()) RecursiveBisectionTask(
      _hg, _context, _top_level, _task_group_id);
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
