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
class RecursiveInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;
  static constexpr bool kahypar_debug = false;
  static constexpr bool enable_heavy_assert = false;

  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct RecursivePartitionResult {
    RecursivePartitionResult() :
      hypergraph(),
      context(),
      mapping(),
      objective(std::numeric_limits<HyperedgeWeight>::max()),
      imbalance(1.0) { }

    explicit RecursivePartitionResult(Context&& c) :
      hypergraph(),
      context(c),
      mapping(),
      objective(std::numeric_limits<HyperedgeWeight>::max()),
      imbalance(1.0) { }

    HyperGraph hypergraph;
    Context context;
    parallel::scalable_vector<HypernodeID> mapping;
    HyperedgeWeight objective;
    double imbalance;
  };

  class RecursiveChildContinuationTask : public tbb::task {

   public:
    RecursiveChildContinuationTask(std::unique_ptr<ICoarsener>&& coarsener,
                                   RecursivePartitionResult& result,
                                   const TaskGroupID task_group_id) :
      _coarsener(std::move(coarsener)),
      _result(result),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      std::unique_ptr<IRefiner> label_propagation =
        LabelPropagationFactory::getInstance().createObject(
          _result.context.refinement.label_propagation.algorithm, _result.hypergraph, _result.context, _task_group_id);
      _coarsener->uncoarsen(label_propagation);

      _result.objective = metrics::objective(_result.hypergraph, _result.context.partition.objective);
      _result.imbalance = metrics::imbalance(_result.hypergraph, _result.context);
      return nullptr;
    }

   private:
    std::unique_ptr<ICoarsener> _coarsener;
    RecursivePartitionResult& _result;
    const TaskGroupID _task_group_id;
  };

  class RecursiveChildTask : public tbb::task {

   public:
    RecursiveChildTask(HyperGraph& hypergraph,
                       const Context& context,
                       RecursivePartitionResult& result,
                       const bool top_level,
                       const TaskGroupID task_group_id,
                       const size_t num_threads,
                       const size_t recursion_number) :
      _hg(hypergraph),
      _context(context),
      _result(result),
      _top_level(top_level),
      _task_group_id(task_group_id),
      _num_threads(num_threads),
      _recursion_number(recursion_number) { }

    tbb::task* execute() override {
      _result = RecursivePartitionResult(setupRecursiveContext());

      // Copy hypergraph
      utils::Timer::instance().start_timer("top_level_hypergraph_copy_" + std::to_string(_recursion_number),
                                          "Top Level Hypergraph Copy " + std::to_string(_recursion_number), true, _top_level);
      auto copy = _hg.copy(_result.context.partition.k, _task_group_id);
      _result.hypergraph = std::move(copy.first);
      _result.mapping = std::move(copy.second);
      utils::Timer::instance().stop_timer("top_level_hypergraph_copy_" + std::to_string(_recursion_number), _top_level);

      // Call multilevel partitioner recursively
      DBG << "Perform recursive multilevel partitioner call with"
          << "k =" << _result.context.partition.k << ","
          << "p =" << _result.context.shared_memory.num_threads << ","
          << "c =" << _result.context.coarsening.contraction_limit << "and"
          << "rep =" << _result.context.initial_partitioning.runs;

      std::unique_ptr<ICoarsener> coarsener =
        CoarsenerFactory::getInstance().createObject(
          _result.context.coarsening.algorithm, _result.hypergraph, _result.context, _task_group_id);
      coarsener->coarsen();

      RecursiveChildContinuationTask& child_continuation = *new(allocate_continuation())
        RecursiveChildContinuationTask(std::move(coarsener), _result, _task_group_id);
      RecursiveTask& recursive_task = *new(child_continuation.allocate_child()) RecursiveTask(
        _result.hypergraph, _result.context, false, _task_group_id);
      child_continuation.set_ref_count(1);
      tbb::task::spawn(recursive_task);
      return nullptr;
    }

   private:
    Context setupRecursiveContext() {
      ASSERT(_num_threads >= 1);
      Context context(_context);

      // Shared Memory Parameters
      context.shared_memory.num_threads = _num_threads;

      // Partitioning Parameters
      bool reduce_k = !_top_level && _context.shared_memory.num_threads < (size_t)_context.partition.k && _context.partition.k > 2;
      if (reduce_k) {
        context.partition.k = std::ceil(((double)context.partition.k) / 2.0);
        context.partition.perfect_balance_part_weights.assign(context.partition.k, 0);
        context.partition.max_part_weights.assign(context.partition.k, 0);
        for (PartitionID part = 0; part < _context.partition.k; ++part) {
          context.partition.perfect_balance_part_weights[part / 2] +=
            _context.partition.perfect_balance_part_weights[part];
          context.partition.max_part_weights[part / 2] +=
            _context.partition.max_part_weights[part];
        }
      }
      context.partition.verbose_output = debug;
      context.type = kahypar::ContextType::initial_partitioning;

      // Coarsening Parameters
      context.coarsening.contraction_limit = std::max(context.partition.k * context.coarsening.contraction_limit_multiplier,
                                                      2 * context.shared_memory.num_threads * context.coarsening.contraction_limit_multiplier);
      context.setupMaximumAllowedNodeWeight(_hg.totalWeight());

      // Initial Partitioning Parameters
      bool is_parallel_recursion = _context.shared_memory.num_threads != context.shared_memory.num_threads;
      context.initial_partitioning.runs = std::max(context.initial_partitioning.runs / (is_parallel_recursion ? 2 : 1), 1UL);

      return context;
    }

    HyperGraph& _hg;
    const Context& _context;
    RecursivePartitionResult& _result;
    const bool _top_level;
    const TaskGroupID _task_group_id;
    const size_t _num_threads;
    const size_t _recursion_number;
  };

  class BisectionContinuationTask : public tbb::task {

   public:
    BisectionContinuationTask(HyperGraph& hypergraph,
                              const Context& context,
                              const HyperedgeWeight current_objective,
                              const PartitionID num_bisections) :
      _hg(hypergraph),
      _context(context),
      _current_objective(current_objective),
      _node_mapping(_hg.initialNumNodes(), kInvalidHypernode),
      _results() {
      _results.reserve(num_bisections);
      for ( PartitionID block = 0; block < num_bisections; ++block ) {
        ASSERT(block < _hg.k());
        _results.emplace_back(_hg.partSize(block));
      }
    }

    tbb::task* execute() override {
      // Apply all bisections to current hypergraph
      PartitionID unbisected_block = (_context.partition.k % 2 == 1 ? (PartitionID) _results.size() : kInvalidPartition);
      for (const HypernodeID& hn : _hg.nodes()) {
        PartitionID from = _hg.partID(hn);
        PartitionID to = kInvalidPartition;
        if (from != unbisected_block) {
          ASSERT(from != kInvalidPartition && from < (PartitionID)_results.size());
          ASSERT(_node_mapping[_hg.originalNodeID(hn)] < _results[from].partition.size());
          to = _results[from].partition[_node_mapping[_hg.originalNodeID(hn)]] == 0 ? 2 * from : 2 * from + 1;
        } else {
          to = _context.partition.k - 1;
        }

        ASSERT(to != kInvalidPartition && to < _hg.k());
        if (from != to) {
          _hg.changeNodePart(hn, from, to);
        }
      }

      HEAVY_INITIAL_PARTITIONING_ASSERT([&] {
          HyperedgeWeight expected_objective = _current_objective;
          HyperedgeWeight actual_objective = metrics::objective(_hg, _context.partition.objective);
          for (size_t i = 0; i < _results.size(); ++i) {
            expected_objective += _results[i].objective;
          }

          if (expected_objective != actual_objective) {
            LOG << V(expected_objective) << V(actual_objective);
            return false;
          }
          return true;
        } ());

      _hg.updateGlobalPartInfos();
      return nullptr;
    }

   private:
    HyperGraph& _hg;
    const Context& _context;
    const HyperedgeWeight _current_objective;

   public:
    parallel::scalable_vector<HypernodeID> _node_mapping;
    parallel::scalable_vector<KaHyParPartitioningResult> _results;
  };

  class BisectionTask : public tbb::task {

   public:
    BisectionTask(HyperGraph& hypergraph,
                  const Context& context,
                  const PartitionID block,
                  KaHyParPartitioningResult& result,
                  parallel::scalable_vector<HypernodeID>& node_mapping) :
      _hg(hypergraph),
      _context(context),
      _block(block),
      _result(result),
      _node_mapping(node_mapping) { }

    tbb::task* execute() override {
      bisectPartition();
      return nullptr;
    }

   private:
    void bisectPartition() {
      kahypar_context_t* kahypar_context = setupKaHyParBisectionContext();
      std::vector<HypernodeWeight> max_part_weights;
      max_part_weights.emplace_back(_context.partition.max_part_weights[2 * _block]);
      max_part_weights.emplace_back(_context.partition.max_part_weights[2 * _block + 1]);
      // Special case, if balance constraint will be violated with this bisection
      // => would cause KaHyPar to exit with failure
      if (max_part_weights[0] + max_part_weights[1] < _hg.partWeight(_block)) {
        HypernodeWeight delta = _hg.partWeight(_block) - (max_part_weights[0] + max_part_weights[1]);
        max_part_weights[0] += std::ceil(((double)delta) / 2.0);
        max_part_weights[1] += std::ceil(((double)delta) / 2.0);
      }
      ASSERT(max_part_weights[0] + max_part_weights[1] >= _hg.partWeight(_block));
      kahypar_set_custom_target_block_weights(2, max_part_weights.data(), kahypar_context);

      // Build KaHyPar Hypergraph Data Structure
      bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
      KaHyParHypergraph kahypar_hypergraph = extractBlockAsKaHyParHypergraph(_hg, _block, _node_mapping, cut_net_splitting);

      // Bisect KaHyPar Hypergraph
      kahypar_partition(kahypar_hypergraph, *reinterpret_cast<kahypar::Context*>(kahypar_context), _result);

      kahypar_context_free(kahypar_context);
    }

    kahypar_context_t* setupKaHyParBisectionContext() {
      kahypar_context_t* kahypar_c = mt_kahypar::setupContext(_context, false, kahypar_debug);
      kahypar::Context& kahypar_context = *reinterpret_cast<kahypar::Context*>(kahypar_c);
      kahypar_context.partition.k = 2;
      kahypar_context.initial_partitioning.technique = kahypar::InitialPartitioningTechnique::flat;
      kahypar_context.initial_partitioning.coarsening.algorithm = kahypar::CoarseningAlgorithm::do_nothing;
      kahypar_context.initial_partitioning.local_search.algorithm = kahypar::RefinementAlgorithm::do_nothing;
      kahypar_context.coarsening.algorithm = kahypar::CoarseningAlgorithm::do_nothing;
      kahypar_context.local_search.algorithm = kahypar::RefinementAlgorithm::do_nothing;
      kahypar_context.initial_partitioning.nruns = _context.initial_partitioning.runs;
      kahypar_context.partition.seed = utils::Randomize::instance().getRandomInt(0, 2000, sched_getcpu());
      return kahypar_c;
    }

    HyperGraph& _hg;
    const Context& _context;
    const PartitionID _block;
    KaHyParPartitioningResult& _result;
    parallel::scalable_vector<HypernodeID>& _node_mapping;
  };

  class RecursiveContinuationTask : public tbb::task {

   public:
    RecursiveContinuationTask(HyperGraph& hypergraph,
                              const Context& context,
                              const bool top_level,
                              const TaskGroupID task_group_id,
                              const bool was_recursion) :
      _hg(hypergraph),
      _context(context),
      _top_level(top_level),
      _task_group_id(task_group_id),
      _was_recursion(was_recursion) { }

    RecursivePartitionResult r1;
    RecursivePartitionResult r2;

    tbb::task* execute() override {
      ASSERT(r1.objective < std::numeric_limits<HyperedgeWeight>::max());

      RecursivePartitionResult best;
      // Choose best partition of both parallel recursion
      bool r1_has_better_quality = r1.objective < r2.objective;
      bool r1_is_balanced = r1.imbalance < r1.context.partition.epsilon;
      bool r2_is_balanced = r2.imbalance < r2.context.partition.epsilon;
      if (!_was_recursion ||
          (r1_has_better_quality && r1_is_balanced) ||
          (r1_is_balanced && !r2_is_balanced) ||
          (r1_has_better_quality && !r1_is_balanced && !r2_is_balanced)) {
        best = std::move(r1);
      } else {
        best = std::move(r2);
      }

      HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective == metrics::objective(best.hypergraph, _context.partition.objective));

      // Apply best partition to hypergraph
      for (const HypernodeID& hn : _hg.nodes()) {
        const HypernodeID original_id = _hg.originalNodeID(hn);
        ASSERT(original_id < best.mapping.size());
        // The partID function of the hypergraph takes a global node id and
        // returns the partition id of the vertex.
        // The mapping function of RecursivePartitionResult object (best.mapping) stores a mapping from
        // the original node id of hypergraph _hg to the original node ids of the copied hypergraph
        // (best.hypergraph).
        // Note, original node ids are the node ids of the input hypergraph and the global node ids are
        // the internal node ids of the hypergraph.
        PartitionID part_id = best.hypergraph.partID(best.hypergraph.globalNodeID(best.mapping[original_id]));
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        _hg.setNodePart(hn, part_id);
      }
      _hg.initializeNumCutHyperedges();
      _hg.updateGlobalPartInfos();

      // The hypergraph is now partitioned into the number of blocks of the recursive context (best.context.partition.k).
      // Based on wheter we reduced k in recursion, we have to bisect the blocks of the partition
      // in the desired number of blocks of the current context (_context.partition.k).

      HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective == metrics::objective(_hg, _context.partition.objective));

      // Bisect all blocks of best partition, if we are not on the top level of recursive initial partitioning
      // and the number of threads is small than k
      bool perform_bisections = !_top_level && _context.shared_memory.num_threads < (size_t)_context.partition.k;
      if (perform_bisections) {
        BisectionContinuationTask& bisection_continuation = *new(allocate_continuation())
          BisectionContinuationTask(_hg, _context, best.objective, _context.partition.k / 2);
        bisection_continuation.set_ref_count(_context.partition.k / 2 );
        for (PartitionID block = 0; block < _context.partition.k / 2; ++block) {
          tbb::task::spawn(*new(bisection_continuation.allocate_child()) BisectionTask(
            _hg, _context, block, bisection_continuation._results[block], bisection_continuation._node_mapping));
        }
      }
      return nullptr;
    }

   private:
    HyperGraph& _hg;
    const Context& _context;
    const bool _top_level;
    const TaskGroupID _task_group_id;
    const bool _was_recursion;
  };

  class RecursiveTask : public tbb::task {

   public:
    RecursiveTask(HyperGraph& hypergraph,
                  const Context& context,
                  const bool top_level,
                  const TaskGroupID task_group_id) :
      _hg(hypergraph),
      _context(context),
      _top_level(top_level),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      if (_context.shared_memory.num_threads == 1 &&
          _context.coarsening.contraction_limit == 2 * _context.coarsening.contraction_limit_multiplier) {
        // Base Case -> Bisect Hypergraph
        initialBisection();
        return nullptr;
      }

      // We do parallel recursion, if the contract limit is equal to 2 * p * t
      // ( where p is the number of threads and t the contract limit multiplier )
      bool do_parallel_recursion = _context.coarsening.contraction_limit /
                                  (2 * _context.coarsening.contraction_limit_multiplier) ==
                                  _context.shared_memory.num_threads;
      if (do_parallel_recursion) {
        // Perform parallel recursion
        size_t num_threads_1 = std::ceil(((double) std::max(_context.shared_memory.num_threads, 2UL)) / 2.0);
        size_t num_threads_2 = std::floor(((double) std::max(_context.shared_memory.num_threads, 2UL)) / 2.0);
        auto tbb_recursion_task_groups = TBB::instance().create_tbb_task_groups_for_recursion();

        RecursiveContinuationTask& recursive_continuation = *new(allocate_continuation())
          RecursiveContinuationTask(_hg, _context, _top_level, _task_group_id, true);
        RecursiveChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveChildTask(
          _hg, _context, recursive_continuation.r1, _top_level, tbb_recursion_task_groups.first, num_threads_1, 0);
        RecursiveChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) RecursiveChildTask(
          _hg, _context, recursive_continuation.r2, _top_level, tbb_recursion_task_groups.second, num_threads_2, 0);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else {
        RecursiveContinuationTask& recursive_continuation = *new(allocate_continuation())
          RecursiveContinuationTask(_hg, _context, _top_level, _task_group_id, false);
        RecursiveChildTask& recursion = *new(recursive_continuation.allocate_child()) RecursiveChildTask(
          _hg, _context, recursive_continuation.r1, _top_level, _task_group_id, _context.shared_memory.num_threads, 0);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion);
      }

      return nullptr;
    }

   private:
    void initialBisection() {
      ASSERT(_context.partition.k == 2);
      ASSERT(_context.partition.max_part_weights.size() == 2);
      ASSERT(_context.shared_memory.num_threads == 1);

      kahypar_context_t* kahypar_context = setupKaHyParBisectionContext();
      kahypar_set_custom_target_block_weights(2, _context.partition.max_part_weights.data(), kahypar_context);

      // Compute node mapping that maps all hypernodes of the coarsened hypergraph
      // to a consecutive range of node ids
      HypernodeID num_vertices = 0;
      std::vector<HypernodeID> node_mapping(_hg.initialNumNodes(), kInvalidHypernode);
      for (const HypernodeID& hn : _hg.nodes()) {
        ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
        node_mapping[_hg.originalNodeID(hn)] = num_vertices++;
      }

      // Build KaHyPar Hypergraph Data Structure
      const KaHyParHypergraph kahypar_hypergraph = convertToKaHyParHypergraph(_hg, node_mapping);
      KaHyParPartitioningResult result(num_vertices);

      // Bisect KaHyPar Hypergraph
      kahypar_partition(kahypar_hypergraph, *reinterpret_cast<kahypar::Context*>(kahypar_context), result);

      // Apply partition to MT-KaHyPar Hypergraph
      for (HypernodeID u = 0; u < result.partition.size(); ++u) {
        ASSERT(u < kahypar_hypergraph.reverse_mapping.size());
        HypernodeID hn = kahypar_hypergraph.reverse_mapping[u];
        ASSERT(node_mapping[_hg.originalNodeID(hn)] != kInvalidHypernode);
        PartitionID part_id = result.partition[u];
        ASSERT(part_id >= 0 && part_id < _context.partition.k);
        _hg.setNodePart(hn, part_id);
      }
      _hg.updateGlobalPartInfos();
      _hg.initializeNumCutHyperedges();

      kahypar_context_free(kahypar_context);
    }

    kahypar_context_t* setupKaHyParBisectionContext() {
      kahypar_context_t* kahypar_c = mt_kahypar::setupContext(_context, false, kahypar_debug);
      kahypar::Context& kahypar_context = *reinterpret_cast<kahypar::Context*>(kahypar_c);
      kahypar_context.partition.k = 2;
      kahypar_context.initial_partitioning.technique = kahypar::InitialPartitioningTechnique::flat;
      kahypar_context.initial_partitioning.coarsening.algorithm = kahypar::CoarseningAlgorithm::do_nothing;
      kahypar_context.initial_partitioning.local_search.algorithm = kahypar::RefinementAlgorithm::do_nothing;
      kahypar_context.coarsening.algorithm = kahypar::CoarseningAlgorithm::do_nothing;
      kahypar_context.local_search.algorithm = kahypar::RefinementAlgorithm::do_nothing;
      kahypar_context.initial_partitioning.nruns = _context.initial_partitioning.runs;
      kahypar_context.partition.seed = utils::Randomize::instance().getRandomInt(0, 2000, sched_getcpu());
      return kahypar_c;
    }

    HyperGraph& _hg;
    const Context& _context;
    const bool _top_level;
    const TaskGroupID _task_group_id;
  };

 public:
  RecursiveInitialPartitionerT(HyperGraph& hypergraph,
                               const Context& context,
                               const bool top_level,
                               const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  RecursiveInitialPartitionerT(const RecursiveInitialPartitionerT&) = delete;
  RecursiveInitialPartitionerT(RecursiveInitialPartitionerT&&) = delete;
  RecursiveInitialPartitionerT & operator= (const RecursiveInitialPartitionerT &) = delete;
  RecursiveInitialPartitionerT & operator= (RecursiveInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    if (_top_level) {
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    RecursiveTask& root_recursive_task = *new(tbb::task::allocate_root()) RecursiveTask(
      _hg, _context, _top_level, _task_group_id);
    tbb::task::spawn_root_and_wait(root_recursive_task);

    if (_top_level) {
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

template <typename TypeTraits>
PartitionID RecursiveInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID RecursiveInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using RecursiveInitialPartitioner = RecursiveInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
