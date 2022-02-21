/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/deep_multilevel.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "tbb/parallel_invoke.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"

#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"


/*!
 * DEEP MULTILEVEL
 * For reason of simplicity we assume in the following description of the algorithm that
 * the number of threads p and the number of blocks k is a power of 2 and p = k. The deep
 * partitioning algorithm is invoked, if the number of vertices is 2 * c * p (where c is our
 * contraction limit multiplier).
 * The deep multilevel algorithm starts by performing parallel coarsening with p threads
 * until c * p vertices are reached. Afterwards, the hypergraph is copied and the hypergraphs
 * are recursively coarsened with p / 2 threads each. Once p = 1 (and the contraction limit is 2 * c)
 * we initially bisect the hypergraph in two blocks. After initial partitioning each thread uncontracts
 * its hypergraph (and performs refinement) until 4 * c hypernodes are rechead. Afterwards, we choose the
 * best partition of both recursions and further bisect each block of the partition to obtain a 4-way
 * partition and continue uncontraction with 2 threads until 8 * c hypernodes. This is repeated until
 * we obtain a k-way partition of the hypergraph.
 * Note, the deep multilevel algorithm is written in TBBInitializer continuation style. The TBBInitializer continuation
 * style is especially useful for recursive patterns. Each task defines its continuation task. A continuation
 * task defines how computation should continue, if all its child tasks are completed. As a consequence,
 * tasks can be spawned without waiting for their completion, because the continuation task is automatically
 * invoked if all child tasks are terminated. Therefore, no thread will waste CPU time while waiting for
 * their recursive tasks to complete.
 *
 * Implementation Details
 * ----------------------
 * The deep multilevel algorithm starts by spawning the root DeepPartitionTask. The DeepPartitionTask spawns
 * two DeepPartitionChildTask. Within such a task the hypergraph is copied and coarsened to the next desired contraction limit.
 * Once that contraction limit is reached the DeepPartitionChildTask spawns again one DeepPartitionTask. Once the DeepPartitionTask
 * of a DeepPartitionChildTask terminates, the DeepChildContinuationTask starts and uncontracts the hypergraph to
 * its original size (and also performs refinement). Once both DeepPartitionChildTask of a DeepPartitionTask terminates, the
 * DeepPartitionContinuationTask starts and chooses the best partition of both recursions and spawns for each block
 * a DeepBisectionTask. The DeepBisectionTask performs a initial partition call to bisect exactly one block of the current
 * partition. Once all DeepBisectionTasks terminates, the DeepBisectionContinuationTask starts and applies all bisections to the
 * current hypergraph.
 */

namespace mt_kahypar {

  struct DeepPartitionResult {
    DeepPartitionResult() :
            hypergraph(),
            partitioned_hypergraph(),
            mapping(),
            context(),
            objective(std::numeric_limits<HyperedgeWeight>::max()),
            imbalance(1.0) { }

    explicit DeepPartitionResult(Context&& c) :
            hypergraph(),
            partitioned_hypergraph(),
            mapping(),
            context(c),
            objective(std::numeric_limits<HyperedgeWeight>::max()),
            imbalance(1.0) { }

    explicit DeepPartitionResult(const Context& c) :
            hypergraph(),
            partitioned_hypergraph(),
            mapping(),
            context(c),
            objective(std::numeric_limits<HyperedgeWeight>::max()),
            imbalance(1.0) { }

    Hypergraph hypergraph;
    PartitionedHypergraph partitioned_hypergraph;
    parallel::scalable_vector<HypernodeID> mapping;
    Context context;
    HyperedgeWeight objective;
    double imbalance;
  };

  struct OriginalHypergraphInfo {

    double computeAdaptiveEpsilon(const PartitionID current_k) const {
      return std::min(0.99, std::max(std::pow(1.0 + original_epsilon, 1.0 /
        log2(ceil(static_cast<double>(original_k) / static_cast<double>(current_k)) + 1.0)) - 1.0,0.0));
    }

    const PartitionID original_k;
    const double original_epsilon;
  };

  /*!
   * Continuation task for the deep child task. It is automatically called
   * after the deep child task terminates and responsible for uncontracting
   * the hypergraph.
   */
  class DeepChildContinuationTask : public tbb::task {

  public:
    DeepChildContinuationTask(DeepPartitionResult& result) :
            _coarsener(nullptr),
            _sparsifier(nullptr),
            _result(result) {
      bool nlevel = _result.context.coarsening.algorithm == CoarseningAlgorithm::nlevel_coarsener;
      _uncoarseningData = std::make_shared<UncoarseningData>(nlevel, _result.hypergraph, _result.context);
      _coarsener = CoarsenerFactory::getInstance().createObject(
              _result.context.coarsening.algorithm, _result.hypergraph, _result.context, *_uncoarseningData);
      _sparsifier = HypergraphSparsifierFactory::getInstance().createObject(
              _result.context.sparsification.similiar_net_combiner_strategy, _result.context);

    }

    tbb::task* execute() override {
      if ( _sparsifier->isSparsified() ) {
        // In that case, the sparsified hypergraph generated by the
        // heavy hyperedge remover was used for initial partitioning.
        // => Partition has to mapped from sparsified hypergraph to
        // coarsest partitioned hypergraph.
        _sparsifier->undoSparsification(_coarsener->coarsestPartitionedHypergraph());
      }

      _coarsener.reset();

      // Uncontraction
      std::unique_ptr<IRefiner> label_propagation =
              LabelPropagationFactory::getInstance().createObject(
                      _result.context.refinement.label_propagation.algorithm, _result.hypergraph,
                      _result.context);
      std::unique_ptr<IRefiner> fm =
              FMFactory::getInstance().createObject(
                      _result.context.refinement.fm.algorithm, _result.hypergraph,
                      _result.context);

      if (_uncoarseningData->nlevel) {
        _uncoarsener = std::make_unique<NLevelUncoarsener>(_result.hypergraph, _result.context, *_uncoarseningData);
      } else {
        _uncoarsener = std::make_unique<MultilevelUncoarsener>(_result.hypergraph, _result.context, *_uncoarseningData);
      }
      _result.partitioned_hypergraph = _uncoarsener->uncoarsen(label_propagation, fm);

      // Compute metrics
      _result.objective = metrics::objective(_result.partitioned_hypergraph, _result.context.partition.objective);
      _result.imbalance = metrics::imbalance(_result.partitioned_hypergraph, _result.context);
      return nullptr;
    }

  public:
    std::unique_ptr<ICoarsener> _coarsener;
    std::unique_ptr<IHypergraphSparsifier> _sparsifier;
    std::unique_ptr<IUncoarsener> _uncoarsener;
    std::shared_ptr<UncoarseningData> _uncoarseningData;

  private:
    DeepPartitionResult& _result;
  };

  /*!
 * Continuation task for the deep bisection task. The deep bisection continuation task
 * is called after all deep bisection tasks terminated and is responsible for applying
 * all bisections done by the deep bisection tasks to the current hypergraph.
 */
  class DeepBisectionContinuationTask : public tbb::task {
    static constexpr bool enable_heavy_assert = false;
  public:
    DeepBisectionContinuationTask(PartitionedHypergraph& hypergraph,
                                  const Context& context,
                                  const HyperedgeWeight current_objective,
                                  const PartitionID num_bisections) :
            _hg(hypergraph),
            _context(context),
            _current_objective(current_objective),
            _results() {
      _results.reserve(num_bisections);
      for ( PartitionID block = 0; block < num_bisections; ++block ) {
        _results.emplace_back(_context);
      }
    }

    tbb::task* execute() override {
      // Apply all bisections to current hypergraph
      PartitionID unbisected_block = (_context.partition.k % 2 == 1 ? (PartitionID) _results.size() : kInvalidPartition);
      _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID from = _hg.partID(hn);
        PartitionID to = kInvalidPartition;
        if ( from != unbisected_block ) {
          ASSERT(from != kInvalidPartition && static_cast<size_t>(from) < _results.size());
          ASSERT(hn < _results[from].mapping.size());
          const PartitionedHypergraph& from_hg = _results[from].partitioned_hypergraph;
          to = from_hg.partID(_results[from].mapping[hn]) == 0 ? 2 * from : 2 * from + 1;
        } else {
          to = _context.partition.k - 1;
        }

        ASSERT(to != kInvalidPartition && to < _hg.k());
        if (from != to) {
          _hg.changeNodePart(hn, from, to);
        }
      });

      HEAVY_INITIAL_PARTITIONING_ASSERT([&] {
        HyperedgeWeight expected_objective = _current_objective;
        HyperedgeWeight actual_objective = metrics::objective(_hg, _context.partition.objective);
        for (size_t i = 0; i < _results.size(); ++i) {
          expected_objective += metrics::objective(
                  _results[i].partitioned_hypergraph, _context.partition.objective);
        }

        if (expected_objective != actual_objective) {
          LOG << V(expected_objective) << V(actual_objective);
          return false;
        }
        return true;
      } ());

      return nullptr;
    }

  private:
    PartitionedHypergraph& _hg;
    const Context& _context;
    const HyperedgeWeight _current_objective;

  public:
    parallel::scalable_vector<DeepPartitionResult> _results;
  };


  /*!
   * A deep bisection task is started after we return from the recursion. It is
   * responsible for bisecting one block of the current k'-way partition (k' < k).
   */
  class DeepBisectionTask : public tbb::task {

  public:
    DeepBisectionTask(PartitionedHypergraph& hypergraph,
                      const PartitionID block,
                      DeepPartitionResult& result) :
            _hg(hypergraph),
            _stable_construction_of_incident_edges(
                result.context.preprocessing.stable_construction_of_incident_edges),
            _block(block),
            _result(result) { }

    tbb::task* execute() override {
      // Setup Initial Partitioning Context
      std::vector<HypernodeWeight> perfect_balance_part_weights;
      std::vector<HypernodeWeight> max_part_weights;
      perfect_balance_part_weights.emplace_back(_result.context.partition.perfect_balance_part_weights[2 * _block]);
      perfect_balance_part_weights.emplace_back(_result.context.partition.perfect_balance_part_weights[2 * _block + 1]);
      max_part_weights.emplace_back(_result.context.partition.max_part_weights[2 * _block]);
      max_part_weights.emplace_back(_result.context.partition.max_part_weights[2 * _block + 1]);
      _result.context.partition.perfect_balance_part_weights = std::move(perfect_balance_part_weights);
      _result.context.partition.max_part_weights = std::move(max_part_weights);
      _result.context.partition.k = 2;

      // Extract Block of Hypergraph
      bool cut_net_splitting = _result.context.partition.objective == kahypar::Objective::km1;
      auto tmp_hypergraph = _hg.extract(_block, cut_net_splitting, _stable_construction_of_incident_edges);
      _result.hypergraph = std::move(tmp_hypergraph.first);
      _result.mapping = std::move(tmp_hypergraph.second);
      _result.partitioned_hypergraph = PartitionedHypergraph(
              2, _result.hypergraph, parallel_tag_t());

      if ( _result.hypergraph.initialNumNodes() > 0 ) {
        // Spawn Initial Partitioner
        PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
                PoolInitialPartitionerContinuation(
                _result.partitioned_hypergraph, _result.context);
        spawn_initial_partitioner(ip_continuation);
      }
      return nullptr;
    }

  private:
    PartitionedHypergraph& _hg;
    bool _stable_construction_of_incident_edges;
    const PartitionID _block;
    DeepPartitionResult& _result;
  };



  /*!
 * Continuation task for the deep partition task. The continuation task
 * is called after all child tasks of the deep partition task terminated
 * and is responsible for choosing the best partition of the child tasks
 * and spawn deep bisection tasks to further transform the k'-way partition into
 * a 2*k'-way partition.
 */
  class DeepPartitionContinuationTask : public tbb::task {
    static constexpr bool enable_heavy_assert = false;
  public:
    DeepPartitionContinuationTask(const OriginalHypergraphInfo original_hypergraph_info,
                                  PartitionedHypergraph& hypergraph,
                                  const Context& context,
                                  const bool was_recursion,
                                  const bool is_top_level) :
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _was_recursion(was_recursion),
            _is_top_level(is_top_level) { }

    DeepPartitionResult r1;
    DeepPartitionResult r2;

    tbb::task* execute() override {
      ASSERT(r1.objective < std::numeric_limits<HyperedgeWeight>::max());

      DeepPartitionResult best;
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
      // Note, we move r1 or r2 into best, both contain the the
      // hypergraph and the partitioned hypergraph, whereas the
      // partitioned hypergraph contains a pointer to the hypergraph.
      // Moving r1 or r2 invalidates the pointer to the original
      // hypergraph. Therefore, we explicitly set it here.
      best.partitioned_hypergraph.setHypergraph(best.hypergraph);

      HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective ==
                                        metrics::objective(best.partitioned_hypergraph, _context.partition.objective));

      // Apply best partition to hypergraph
      _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        PartitionID part_id = best.partitioned_hypergraph.partID(hn);
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        _hg.setOnlyNodePart(hn, part_id);
      });
      _hg.initializePartition();

      // The hypergraph is now partitioned into the number of blocks of the recursive context (best.context.partition.k).
      // Based on whether we reduced k in recursion, we have to bisect the blocks of the partition
      // in the desired number of blocks of the current context (_context.partition.k).

      HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective == metrics::objective(_hg, _context.partition.objective));

      // Bisect all blocks of best partition, if we are not on the top level of recursive initial partitioning
      // and the number of threads is small than k
      bool perform_bisections = !_is_top_level &&
        _context.shared_memory.num_threads < (size_t)_context.partition.k;
      if (perform_bisections) {
        DeepBisectionContinuationTask& bisection_continuation = *new(allocate_continuation())
                DeepBisectionContinuationTask(_hg, _context,
                                              best.objective, _context.partition.k / 2);
        bisection_continuation.set_ref_count(_context.partition.k / 2 );
        for (PartitionID block = 0; block < _context.partition.k / 2; ++block) {
          tbb::task::spawn(*new(bisection_continuation.allocate_child()) DeepBisectionTask(
                  _hg, block, bisection_continuation._results[block]));
        }
      }
      return nullptr;
    }

  private:
    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context& _context;
    const bool _was_recursion;
    const bool _is_top_level;
  };



  /*!
   * The deep partition task contains the base case for initial bisecting the hypergraph
   * (if p = 1) and performing recursion by calling the deep partition child tasks.
   */
  class DeepPartitionTask : public tbb::task {

  public:
    DeepPartitionTask(const OriginalHypergraphInfo original_hypergraph_info,
                      PartitionedHypergraph& hypergraph,
                      const Context& context,
                      bool is_top_level) :
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _is_top_level(is_top_level) { }

    tbb::task* execute() override ;

  private:
    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context& _context;
    bool _is_top_level;
  };


  /*!
   * The recursive child task is responsible for copying the hypergraph
   * and coarsen the hypergraph until the next contraction limit is reached.
   */
  class DeepPartitionChildTask : public tbb::task {
    static constexpr bool debug = false;
  public:
    DeepPartitionChildTask(const OriginalHypergraphInfo original_hypergraph_info,
                           PartitionedHypergraph& hypergraph,
                           const Context& context,
                           DeepPartitionResult& result,
                           const size_t num_threads,
                           const size_t recursion_number,
                           const double degree_of_parallelism,
                           bool is_top_level) :
            _original_hypergraph_info(original_hypergraph_info),
            _hg(hypergraph),
            _context(context),
            _result(result),
            _num_threads(num_threads),
            _recursion_number(recursion_number),
            _degree_of_parallelism(degree_of_parallelism),
            _is_top_level(is_top_level) { }

    tbb::task* execute() override {
      // Copy hypergraph
      _result = DeepPartitionResult(setupRecursiveContext(_is_top_level));
      _result.hypergraph = _hg.hypergraph().copy(parallel_tag_t());

      DBG << "Perform recursive multilevel partitioner call with"
          << "k =" << _result.context.partition.k << ","
          << "p =" << _result.context.shared_memory.num_threads << ","
          << "c =" << _result.context.coarsening.contraction_limit << "and"
          << "rep =" << _result.context.initial_partitioning.runs;

      DeepChildContinuationTask& child_continuation = *new(allocate_continuation())
              DeepChildContinuationTask(_result);

      // Coarsening
      child_continuation._coarsener->coarsen();

      // Call deep multilevel algorithm
      if ( _context.useSparsification() ) {
        // Sparsify Hypergraph, if heavy hyperedge removal is enabled
        child_continuation._sparsifier->sparsify(child_continuation._coarsener->coarsestHypergraph());
      }

      if ( child_continuation._sparsifier->isSparsified() ) {
        deepPartition(child_continuation._sparsifier->sparsifiedPartitionedHypergraph(), child_continuation);
      } else {
        deepPartition(child_continuation._coarsener->coarsestPartitionedHypergraph(), child_continuation);
      }

      return nullptr;
    }

  private:
    void deepPartition(PartitionedHypergraph& partitioned_hypergraph,
                       DeepChildContinuationTask& child_continuation) {
      DeepPartitionTask& recursive_task = *new(child_continuation.allocate_child()) DeepPartitionTask(
              _original_hypergraph_info, partitioned_hypergraph, _result.context, false);
      child_continuation.set_ref_count(1);
      tbb::task::spawn(recursive_task);
    }

    Context setupRecursiveContext(bool is_top_level) {
      ASSERT(_num_threads >= 1);
      Context context(_context);

      if (!is_top_level) {
        context.type = kahypar::ContextType::initial_partitioning;
      }
      context.partition.verbose_output = debug;

      // Shared Memory Parameters
      context.shared_memory.num_threads = _num_threads;
      context.shared_memory.degree_of_parallelism *= _degree_of_parallelism;

      // Partitioning Parameters
      bool reduce_k = !is_top_level &&
        _context.shared_memory.num_threads < (size_t)_context.partition.k && _context.partition.k > 2;
      if (reduce_k) {
        context.partition.k = std::ceil(((double)context.partition.k) / 2.0);
        context.partition.perfect_balance_part_weights.assign(context.partition.k, 0);
        context.partition.max_part_weights.assign(context.partition.k, 0);
        for (PartitionID part = 0; part < _context.partition.k; ++part) {
          context.partition.perfect_balance_part_weights[part / 2] +=
                  _context.partition.perfect_balance_part_weights[part];
        }

        context.partition.epsilon = _original_hypergraph_info.computeAdaptiveEpsilon(context.partition.k);
        for (PartitionID part = 0; part < context.partition.k; ++part) {
          context.partition.max_part_weights[part] = std::ceil(( 1.0 + context.partition.epsilon ) *
                                                                context.partition.perfect_balance_part_weights[part]);
        }
      }

      // Coarsening Parameters
      context.coarsening.contraction_limit = std::max(
              context.partition.k * context.coarsening.contraction_limit_multiplier,
              2 * ID(context.shared_memory.num_threads) *
              context.coarsening.contraction_limit_multiplier);
      context.setupMaximumAllowedNodeWeight(_hg.totalWeight());
      context.setupSparsificationParameters();
      context.setupThreadsPerFlowSearch();

      // Initial Partitioning Parameters
      bool is_parallel_recursion = _context.shared_memory.num_threads != context.shared_memory.num_threads;
      context.initial_partitioning.runs = std::max(context.initial_partitioning.runs / (is_parallel_recursion ? 2 : 1), 1UL);

      return context;
    }

    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHypergraph& _hg;
    const Context& _context;
    DeepPartitionResult& _result;
    const size_t _num_threads;
    const size_t _recursion_number;
    const double _degree_of_parallelism;
    const bool _is_top_level;
  };


  tbb::task* DeepPartitionTask::execute() {
    if (_context.shared_memory.num_threads == 1 &&
        _context.coarsening.contraction_limit == 2 * _context.coarsening.contraction_limit_multiplier) {
      // Base Case -> Bisect Hypergraph
      ASSERT(_context.partition.k == 2);
      ASSERT(_context.partition.max_part_weights.size() == 2);
      PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
              PoolInitialPartitionerContinuation(_hg, _context);
      spawn_initial_partitioner(ip_continuation);
    } else {
      // We do parallel recursion, if the contract limit is equal to 2 * p * t
      // ( where p is the number of threads and t the contract limit multiplier )
      bool do_parallel_recursion = _context.coarsening.contraction_limit /
                                   (2 * _context.coarsening.contraction_limit_multiplier) ==
                                   _context.shared_memory.num_threads;
      if (do_parallel_recursion) {
        // Perform parallel recursion
        size_t num_threads_1 = std::ceil(((double) std::max(_context.shared_memory.num_threads, 2UL)) / 2.0);
        size_t num_threads_2 = std::floor(((double) std::max(_context.shared_memory.num_threads, 2UL)) / 2.0);

        DeepPartitionContinuationTask& recursive_continuation = *new(allocate_continuation())
                DeepPartitionContinuationTask(_original_hypergraph_info, _hg, _context, true, _is_top_level);
        DeepPartitionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) DeepPartitionChildTask(
                _original_hypergraph_info, _hg, _context, recursive_continuation.r1,
                num_threads_1, 0, 0.5, _is_top_level);
        DeepPartitionChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) DeepPartitionChildTask(
                _original_hypergraph_info, _hg, _context, recursive_continuation.r2,
                num_threads_2, 0, 0.5, _is_top_level);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else {
        DeepPartitionContinuationTask& recursive_continuation = *new(allocate_continuation())
                DeepPartitionContinuationTask(_original_hypergraph_info, _hg, _context, false, _is_top_level);
        DeepPartitionChildTask& recursion = *new(recursive_continuation.allocate_child()) DeepPartitionChildTask(
                _original_hypergraph_info, _hg, _context, recursive_continuation.r1,
                _context.shared_memory.num_threads, 0, 1.0, _is_top_level);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion);
      }
    }
    return nullptr;
  }


namespace deep_multilevel {
  PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context) {
    PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph, parallel_tag_t());
    partition(partitioned_hypergraph, context);
    return partitioned_hypergraph;
  }

  void partition(PartitionedHypergraph& hypergraph, const Context& context) {
    if (context.partition.mode == Mode::deep_multilevel) {
      utils::Timer::instance().start_timer("deep", "Deep Multilevel");
    }
    if (context.type == kahypar::ContextType::main) {
      parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    DeepPartitionTask& root_recursive_task = *new(tbb::task::allocate_root()) DeepPartitionTask(
            OriginalHypergraphInfo { context.partition.k, context.partition.epsilon },
            hypergraph, context, context.partition.mode == Mode::deep_multilevel);
    tbb::task::spawn_root_and_wait(root_recursive_task);

    if (context.partition.num_vcycles > 0 && context.type == kahypar::ContextType::main) {
      multilevel::partitionVCycle(hypergraph.hypergraph(), hypergraph, context);
    }

    if (context.type == kahypar::ContextType::main) {
      parallel::MemoryPool::instance().activate_unused_memory_allocations();
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
    if (context.partition.mode == Mode::deep_multilevel) {
      utils::Timer::instance().stop_timer("deep");
    }
  }
} // namespace deep_multilevel
} // namepace mt_kahypar
