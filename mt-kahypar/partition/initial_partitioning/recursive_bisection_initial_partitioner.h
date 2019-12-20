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

 public:
  RecursiveBisectionInitialPartitionerT(HyperGraph& hypergraph,
                                        const Context& context,
                                        const bool top_level,
                                        TBB& tbb_arena) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _tbb_arena(tbb_arena) { }

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

    recursiveBisection(_hg, _context, _top_level, _tbb_arena);

    if (_top_level) {
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

  void recursiveBisection(HyperGraph& hypergraph, const Context& context, const bool top_level, TBB& tbb_arena) {
    ASSERT(_context.partition.k >= 2);

    PartitionID num_blocks_part_0 = context.partition.k / 2 + (context.partition.k % 2 != 0 ? 1 : 0);
    PartitionID num_blocks_part_1 = context.partition.k / 2;
    BlockRange range_0 = std::make_pair(0, num_blocks_part_0);
    BlockRange range_1 = std::make_pair(num_blocks_part_0, num_blocks_part_0 + num_blocks_part_1);

    // Bisecting the hypergraph into two blocks
    utils::Timer::instance().start_timer("top_level_bisection", "Top Level Bisection", false, top_level);
    bisect(hypergraph, context, 0, num_blocks_part_0, top_level, tbb_arena);
    utils::Timer::instance().stop_timer("top_level_bisection", top_level);

    bool do_parallel_recursion = num_blocks_part_0 >= 2 && num_blocks_part_1 >= 2;
    if ( do_parallel_recursion ) {

      if ( context.shared_memory.num_threads > 1 ) {
        // In case we have to partition both blocks from the bisection further into
        // more than one block, we call the recursive bisection initial partitioner
        // recursively in parallel with half of the number of threads.
        size_t num_threads_0 = std::max(context.shared_memory.num_threads / 2 +
          (context.shared_memory.num_threads % 2 == 1 ? 1 : 0), 1UL);
        size_t num_threads_1 = std::max(context.shared_memory.num_threads / 2, 1UL);

        DBG << "Current k = " << context.partition.k << ","
            << "Parallel Recursion 0: k =" << num_blocks_part_0 << ", num_threads =" << num_threads_0
            << "Parallel Recursion 1: k =" << num_blocks_part_1 << ", num_threads =" << num_threads_1;

        auto tbb_splitted_arena = tbb_arena.split_tbb_numa_arena(num_threads_0, num_threads_1);
        tbb::parallel_invoke([&] {
          utils::Timer::instance().start_timer("top_level_recursion_0", "Top Level Recursion 0", true, top_level);
          recursivelyBisectBlock(hypergraph, context, 0, num_threads_0, range_0, *tbb_splitted_arena.first);
          utils::Timer::instance().stop_timer("top_level_recursion_0", top_level);
        }, [&] {
          utils::Timer::instance().start_timer("top_level_recursion_1", "Top Level Recursion 1", true, top_level);
          recursivelyBisectBlock(hypergraph, context, num_blocks_part_0, num_threads_1, range_1, *tbb_splitted_arena.second);
          utils::Timer::instance().stop_timer("top_level_recursion_1", top_level);
        });
        tbb_splitted_arena.first->terminate();
        tbb_splitted_arena.second->terminate();
      } else {
        DBG << "Current k = " << context.partition.k << ","
            << "Sequential Recursion 0: k =" << num_blocks_part_0 << ", num_threads =" << context.shared_memory.num_threads
            << "Sequential Recursion 1: k =" << num_blocks_part_1 << ", num_threads =" << context.shared_memory.num_threads;

        utils::Timer::instance().start_timer("top_level_recursion_0", "Top Level Recursion 0", true, top_level);
        recursivelyBisectBlock(hypergraph, context, 0, 1UL, range_0, tbb_arena);
        utils::Timer::instance().stop_timer("top_level_recursion_0", top_level);

        utils::Timer::instance().start_timer("top_level_recursion_1", "Top Level Recursion 1", true, top_level);
        recursivelyBisectBlock(hypergraph, context, num_blocks_part_0, 1UL, range_1, tbb_arena);
        utils::Timer::instance().stop_timer("top_level_recursion_1", top_level);
      }
    } else if ( num_blocks_part_0 >= 2 ) {
      // In case only the first block has to be partitioned into more than one block, we call
      // the recursive bisection initial partitioner recusively on the block 0
      DBG << "Current k = " << context.partition.k << ","
          << "Sequential Recursion 0: k =" << num_blocks_part_0 << ", num_threads =" << context.shared_memory.num_threads;

      utils::Timer::instance().start_timer("top_level_recursion_0", "Top Level Recursion 0", true, top_level);
      recursivelyBisectBlock(hypergraph, context, 0, context.shared_memory.num_threads, range_0, tbb_arena);
      utils::Timer::instance().stop_timer("top_level_recursion_0", top_level);
    }

    hypergraph.updateGlobalPartInfos();
  }

  void recursivelyBisectBlock(HyperGraph& hypergraph,
                              const Context& context,
                              const PartitionID block,
                              const size_t num_threads,
                              const BlockRange& range,
                              TBB& tbb_arena) {
    const PartitionID k = range.second - range.first;
    ASSERT(k >= 2);
    Context rb_context = setupRecursiveBisectionContext(context, k, num_threads, range);

    // Extracts the block, which we want to recursively partition, as new hypergraph
    // and calls the recursive bisection initial partitioner recursively.
    bool cut_net_splitting = context.partition.objective == kahypar::Objective::km1;
    auto copy_hypergraph = hypergraph.copy(k, tbb_arena, block, cut_net_splitting);
    HyperGraph& rb_hypergraph = copy_hypergraph.first;
    auto& mapping = copy_hypergraph.second;
    recursiveBisection(rb_hypergraph, rb_context, false, tbb_arena);

    // Assigns partition to top level hypergraph
    assignPartitionFromRecursionToOriginalHypergraph(hypergraph, rb_hypergraph, mapping, block);
  }

  void bisect(HyperGraph& hypergraph, const Context& context,
              const PartitionID block_0, const PartitionID block_1,
              bool top_level, TBB& tbb_arena) {
    ASSERT(block_0 < context.partition.k);
    ASSERT(block_1 < context.partition.k);
    Context bisection_context = setupBisectionContext(hypergraph.totalWeight(), context);

    utils::Timer::instance().start_timer("top_level_copy", "Top Level Copy", false, top_level);
    auto copy_hypergraph = hypergraph.copy(2, tbb_arena);
    HyperGraph& tmp_hg = copy_hypergraph.first;
    auto& mapping = copy_hypergraph.second;
    utils::Timer::instance().stop_timer("top_level_copy", top_level);

    utils::Timer::instance().start_timer("top_level_parallel_bisection", "Top Level Parallel Bisection", false, top_level);
    // Bisect hypergraph with parallel multilevel bisection
    multilevel::partition(tmp_hg, bisection_context, false, tbb_arena);
    utils::Timer::instance().stop_timer("top_level_parallel_bisection", top_level);

    // Apply partition to hypergraph
    utils::Timer::instance().start_timer("top_level_apply_bisection", "Top Level Apply Bisection", false, top_level);
    for (const HypernodeID& hn : hypergraph.nodes()) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < mapping.size());
      PartitionID part_id = tmp_hg.partID(tmp_hg.globalNodeID(mapping[original_id]));
      ASSERT(part_id != kInvalidPartition && part_id < hypergraph.k());
      if ( part_id == 0 ) {
        hypergraph.setNodePart(hn, block_0);
      } else {
        hypergraph.setNodePart(hn, block_1);
      }
    }
    hypergraph.initializeNumCutHyperedges();
    hypergraph.updateGlobalPartInfos();

    utils::Timer::instance().stop_timer("top_level_apply_bisection", top_level);

    ASSERT(metrics::objective(tmp_hg, context.partition.objective) ==
      metrics::objective(hypergraph, context.partition.objective));
  }

  void assignPartitionFromRecursionToOriginalHypergraph(HyperGraph& original_hg,
                                                        const HyperGraph& rb_hg,
                                                        const parallel::scalable_vector<HypernodeID>& mapping,
                                                        const PartitionID part_id) {
    ASSERT(original_hg.initialNumNodes() == mapping.size());
    for ( const HypernodeID& hn : original_hg.nodes() ) {
      if ( original_hg.partID(hn) == part_id ) {
        const HypernodeID original_id = original_hg.originalNodeID(hn);
        ASSERT(original_id < mapping.size());
        PartitionID to = part_id + rb_hg.partID(rb_hg.globalNodeID(mapping[original_id]));
        ASSERT(to != kInvalidPartition && to < original_hg.k());
        original_hg.changeNodePart(hn, part_id, to);
      }
    }
  }

  Context setupRecursiveBisectionContext(const Context& context,
                                         const PartitionID k,
                                         const size_t num_threads,
                                         const BlockRange& range) {
    Context rb_context(context);
    rb_context.partition.k = k;
    rb_context.shared_memory.num_threads = num_threads;

    ASSERT(range.first < range.second);
    ASSERT(range.second - range.first == k);
    rb_context.partition.perfect_balance_part_weights.assign(k, 0);
    rb_context.partition.max_part_weights.assign(k, 0);
    for ( PartitionID part_id = range.first; part_id < range.second; ++part_id ) {
      rb_context.partition.perfect_balance_part_weights[part_id - range.first] =
        context.partition.perfect_balance_part_weights[part_id];
      rb_context.partition.max_part_weights[part_id - range.first] =
        context.partition.max_part_weights[part_id];
    }

    return rb_context;
  }

  Context setupBisectionContext(const HypernodeWeight total_weight, const Context& context) {
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
    // => causes KaHyPar to exit with failure
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

 private:
  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
  TBB& _tbb_arena;
};

template <typename TypeTraits>
PartitionID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using RecursiveBisectionInitialPartitioner = RecursiveBisectionInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
