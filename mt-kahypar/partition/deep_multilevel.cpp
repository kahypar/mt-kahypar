/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/io/partitioning_output.h"

namespace mt_kahypar {

namespace deep_multilevel {

struct OriginalHypergraphInfo {

  double computeAdaptiveEpsilon(const PartitionID current_k) const {
    return std::min(0.99, std::max(std::pow(1.0 + original_epsilon, 1.0 /
      log2(ceil(static_cast<double>(original_k) / static_cast<double>(current_k)) + 1.0)) - 1.0,0.0));
  }

  const PartitionID original_k;
  const double original_epsilon;
};

namespace tmp {

static constexpr bool enable_heavy_assert = false;
static constexpr bool debug = false;

struct DeepPartitionResult {
  explicit DeepPartitionResult(Context&& c) :
          hypergraph(),
          partitioned_hypergraph(),
          context(c),
          objective(std::numeric_limits<HyperedgeWeight>::max()),
          imbalance(1.0) { }

  explicit DeepPartitionResult(const Context& c) :
          hypergraph(),
          partitioned_hypergraph(),
          context(c),
          objective(std::numeric_limits<HyperedgeWeight>::max()),
          imbalance(1.0) { }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  HyperedgeWeight objective;
  double imbalance;
};

class DeepMultilevelOutput {

  struct PartitioningResult {
    HypernodeID num_nodes;
    HyperedgeID num_hes;
    HypernodeID num_pins;
    HyperedgeWeight objective;
    double imbalance;
  };

  struct Results {
    Results(const PartitionID num_blocks,
            const double eps,
            const PartitionID limit) :
      k(num_blocks),
      epsilon(eps),
      c_limit(limit),
      results(),
      remaining(0),
      copies(0) { }

    void increment_calls() {
      ++remaining;
      ++copies;
    }

    void decrement_remaining() {
      --remaining;
    }

    const PartitionID k;
    const double epsilon;
    const PartitionID c_limit;
    vec<PartitioningResult> results;
    size_t remaining;
    size_t copies;
  };

 public:
  explicit DeepMultilevelOutput(const Context& context) :
    _enabled(context.partition.verbose_output),
    _m(),
    _num_threads(context.shared_memory.num_threads),
    _remaining_bipartitions(context.shared_memory.num_threads),
    _cur_level(0),
    _levels() { }

  void bipartitionDone() {
    if ( _enabled ) {
      --_remaining_bipartitions;
    }
  }

  void addPartitioningCall(const PartitionID k,
                           const double epsilon,
                           const HypernodeID limit) {
    if ( _enabled ) {
      std::lock_guard<std::mutex> lock(_m);
      size_t level = std::numeric_limits<size_t>::max();
      if ( _levels.count(limit) == 0 ) {
        level = _cur_level++;
        _levels[limit] = level;
        _results.emplace_back(k, epsilon, limit);
      } else {
        level = _levels[limit];
      }
      _results[level].increment_calls();
    }
  }

  void addPartitioningResult(const HypernodeID limit,
                             const PartitionedHypergraph& partitioned_hg,
                             const HyperedgeWeight objective,
                             const double imbalance) {
    if ( _enabled ) {
      std::lock_guard<std::mutex> lock(_m);
      ASSERT(_levels.count(limit) == 1);
      const size_t level = _levels[limit];
      _results[level].results.push_back(PartitioningResult {
        partitioned_hg.initialNumNodes(), partitioned_hg.initialNumEdges(),
        partitioned_hg.initialNumPins(), objective, imbalance });
      _results[level].decrement_remaining();

      if ( is_done(level) ) {
        printResult(_results[level]);
      }
    }
  }

 private:

  void printResult(const Results& result) {
    LOG << "";
    LOG << "############# k =" << result.k
        << "- Eps =" << result.epsilon
        << "- Contraction Limit =" << result.c_limit
        << "#############";
    const uint8_t width = 12;
    const size_t best_idx = result.results.size() > 1 ? best_result(result) : 2;
    LOG << BOLD
        << std::left << std::setw(width) << "Thread(s)"
        << std::left << std::setw(width) << "Num Nodes"
        << std::left << std::setw(width) << "Num HEs"
        << std::left << std::setw(width) << "Num Pins"
        << std::left << std::setw(width) << "Objective"
        << std::left << std::setw(width) << "Imbalance" << END;
    for ( size_t i = 0; i < result.results.size(); ++i ) {
      LOG << (i == best_idx ? BOLD : "")
          << std::left << std::setw(width) << thread_num(i, result.results.size(), _num_threads)
          << std::left << std::setw(width) << result.results[i].num_nodes
          << std::left << std::setw(width) << result.results[i].num_hes
          << std::left << std::setw(width) << result.results[i].num_pins
          << std::left << std::setw(width) << result.results[i].objective
          << std::left << std::setw(width)<< result.results[i].imbalance
          << (i == best_idx ? END : "") ;
    }
    LOG << "";
  }

  size_t best_result(const Results& result) {
    size_t best_idx = 0;
    for ( size_t i = 1; i < result.results.size(); ++i ) {
      if ( result.results[best_idx].objective > result.results[i].objective ) {
        best_idx = i;
      }
    }
    return best_idx;
  }

  std::string thread_num(const size_t i, const size_t num_copies, const size_t num_threads) {
    std::stringstream ss;
    if ( num_copies == num_threads ) {
      ss << (i + 1);
    } else {
      const size_t step = num_threads / num_copies + (num_threads % num_copies > 0);
      ss << (i * step + 1) << " - " << std::min((i + 1) * step, num_threads);
    }
    return ss.str();
  }

  bool is_done(const size_t level) const {
    ASSERT(level < _results.size());
    return _remaining_bipartitions == 0 && _results[level].remaining == 0;
  }

  const bool _enabled;
  std::mutex _m;
  const size_t _num_threads;
  std::atomic<size_t> _remaining_bipartitions;
  std::atomic<size_t> _cur_level;
  std::unordered_map<HypernodeID, size_t> _levels;
  vec<Results> _results;

};

Context setupMultilevelContext(const Hypergraph& hypergraph,
                               const Context& context,
                               const OriginalHypergraphInfo& info,
                               const size_t num_threads,
                               const double degree_of_parallelism,
                               const bool is_top_level) {
  ASSERT(num_threads >= 1);
  Context d_context(context);

  if (!is_top_level) {
    d_context.type = ContextType::initial_partitioning;
  }
  d_context.partition.verbose_output =
    context.partition.verbose_output && is_top_level;

  // Shared Memory Parameters
  d_context.shared_memory.num_threads = num_threads;
  d_context.shared_memory.degree_of_parallelism *= degree_of_parallelism;

  // Partitioning Parameters
  bool reduce_k = !is_top_level && context.partition.k > 2 &&
    context.shared_memory.num_threads < static_cast<size_t>(context.partition.k);
  if (reduce_k) {
    d_context.partition.k = std::ceil(((double)d_context.partition.k) / 2.0);
    d_context.partition.perfect_balance_part_weights.assign(d_context.partition.k, 0);
    d_context.partition.max_part_weights.assign(d_context.partition.k, 0);
    for (PartitionID part = 0; part < context.partition.k; ++part) {
      d_context.partition.perfect_balance_part_weights[part / 2] +=
        context.partition.perfect_balance_part_weights[part];
    }

    d_context.partition.epsilon = info.computeAdaptiveEpsilon(d_context.partition.k);
    for (PartitionID part = 0; part < d_context.partition.k; ++part) {
      d_context.partition.max_part_weights[part] =
        std::ceil(( 1.0 + d_context.partition.epsilon ) * d_context.partition.perfect_balance_part_weights[part]);
    }
  }

  // Coarsening Parameters
  d_context.coarsening.contraction_limit =
    std::max(static_cast<size_t>(d_context.partition.k), 2 * d_context.shared_memory.num_threads) *
      d_context.coarsening.contraction_limit_multiplier;
  d_context.setupMaximumAllowedNodeWeight(hypergraph.totalWeight());
  d_context.setupThreadsPerFlowSearch();

  // Initial Partitioning Parameters
  const bool is_parallel_recursion = context.shared_memory.num_threads != d_context.shared_memory.num_threads;
  d_context.initial_partitioning.runs = std::max(
    d_context.initial_partitioning.runs / (is_parallel_recursion ? 2 : 1), 1UL);

  return d_context;
}

void recursively_perform_multilevel_partitioning(PartitionedHypergraph& partitioned_hg,
                                                 DeepPartitionResult& result,
                                                 const OriginalHypergraphInfo& info,
                                                 DeepMultilevelOutput& output,
                                                 const size_t num_threads,
                                                 const bool is_parallel_recursion,
                                                 const bool is_top_level);

struct BipartitioningResult {
  Hypergraph hg;
  PartitionedHypergraph partitioned_hg;
  vec<HypernodeID> mapping;
};

BipartitioningResult bipartion_block(PartitionedHypergraph& partitioned_hg,
                                     const Context& context,
                                     const PartitionID block) {
  // Setup context for bipartitioning call
  Context b_context(context);
  b_context.partition.perfect_balance_part_weights.clear();
  b_context.partition.max_part_weights.clear();
  b_context.partition.perfect_balance_part_weights.emplace_back(context.partition.perfect_balance_part_weights[2 * block]);
  b_context.partition.perfect_balance_part_weights.emplace_back(context.partition.perfect_balance_part_weights[2 * block + 1]);
  b_context.partition.max_part_weights.emplace_back(context.partition.max_part_weights[2 * block]);
  b_context.partition.max_part_weights.emplace_back(context.partition.max_part_weights[2 * block + 1]);
  b_context.partition.k = 2;

  // Extract block of partition
  const bool cut_net_splitting = context.partition.objective == Objective::km1;
  auto tmp_hypergraph = partitioned_hg.extract(block, cut_net_splitting,
    context.preprocessing.stable_construction_of_incident_edges);
  Hypergraph& hg  = tmp_hypergraph.first;
  PartitionedHypergraph bipartitioned_hg(2, hg, parallel_tag_t());
  auto& mapping = tmp_hypergraph.second;

  if ( hg.initialNumNodes() > 0 ) {
    // Bipartition block
    pool::bipartition(bipartitioned_hg, b_context);
  }

  return BipartitioningResult { std::move(hg), std::move(bipartitioned_hg), std::move(mapping) };
}

void deep_multilevel_partitioning(PartitionedHypergraph& partitioned_hg,
                                  const Context& context,
                                  const OriginalHypergraphInfo& info,
                                  DeepMultilevelOutput& output,
                                  const HypernodeID current_num_nodes,
                                  const bool is_top_level) {
  const HypernodeID contraction_limit_for_bipartitioning = 2 * context.coarsening.contraction_limit_multiplier;
  if ( context.shared_memory.num_threads == 1 &&
       context.coarsening.contraction_limit == contraction_limit_for_bipartitioning ) {
    // When only one tread remains and we reach the contraction limit for bipartitioning,
    // we perform initial bipartitioning.
    ASSERT(context.partition.k == 2);
    ASSERT(context.partition.max_part_weights.size() == 2);
    pool::bipartition(partitioned_hg, context);
    output.bipartitionDone();
  } else {
    // The deep multilevel scheme maintains the invariant that t threads process a
    // a hypergraph with at least t * C nodes (C = contraction_limit_for_bipartitioning).
    // If we reach the contraction limit where this invariant is violated, we copy the
    // hypergraph and continue the deep multilevel partitioning on both copies recursively.
    const bool do_parallel_recursion = current_num_nodes ==
      context.shared_memory.num_threads * contraction_limit_for_bipartitioning;
    tbb::task_group tg;
    DeepPartitionResult r1(context);
    DeepPartitionResult r2(context);
    if ( do_parallel_recursion ) {
      size_t num_threads_1 = std::ceil(((double) std::max(context.shared_memory.num_threads, 2UL)) / 2.0);
      size_t num_threads_2 = std::floor(((double) std::max(context.shared_memory.num_threads, 2UL)) / 2.0);
      tg.run([&] { recursively_perform_multilevel_partitioning(
          partitioned_hg, r1, info, output, num_threads_1, true, is_top_level); });
      tg.run([&] { recursively_perform_multilevel_partitioning(
        partitioned_hg, r2, info, output, num_threads_2, true, is_top_level); });
      tg.wait();
    } else {
      recursively_perform_multilevel_partitioning(
        partitioned_hg, r1, info, output, context.shared_memory.num_threads, false, is_top_level);
    }

    DBG << "\nRecursion 1:"
        << "k =" << r1.context.partition.k << "Epsilon =" << r1.context.partition.epsilon
        << "Objective =" << r1.objective << "Imbalance =" << r1.imbalance
        << "\nRecursion 2:"
        << "k =" << r2.context.partition.k << "Epsilon =" << r2.context.partition.epsilon
        << "Objective =" << r2.objective << "Imbalance =" << r2.imbalance;


    // After returning from the recursion, we continue deep multilevel partitioning
    // with the better partition from the two recursive calls;
    DeepPartitionResult best(context);
    bool r1_has_better_quality = r1.objective < r2.objective;
    bool r1_is_balanced = r1.imbalance < r1.context.partition.epsilon;
    bool r2_is_balanced = r2.imbalance < r2.context.partition.epsilon;
    if ((r1_has_better_quality && r1_is_balanced) ||
        (r1_is_balanced && !r2_is_balanced) ||
        (r1_has_better_quality && !r1_is_balanced && !r2_is_balanced)) {
      best = std::move(r1);
    } else {
      best = std::move(r2);
    }
    best.partitioned_hypergraph.setHypergraph(best.hypergraph);

    HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective ==
      metrics::objective(best.partitioned_hypergraph, context.partition.objective));

    // Apply best partition to hypergraph
    partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      PartitionID part_id = best.partitioned_hypergraph.partID(hn);
      ASSERT(part_id != kInvalidPartition && part_id < partitioned_hg.k());
      partitioned_hg.setOnlyNodePart(hn, part_id);
    });
    partitioned_hg.initializePartition();

    HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective ==
      metrics::objective(partitioned_hg, context.partition.objective));

    // As long as the number of threads are smaller than the number of blocks,
    // we bipartition each block until the number of blocks equals the desired
    // number of blocks.
    const bool bipartition_all_blocks = !is_top_level &&
      context.shared_memory.num_threads < (size_t) context.partition.k;
    if (bipartition_all_blocks) {
      vec<BipartitioningResult> result(context.partition.k / 2);
      for (PartitionID block = 0; block < context.partition.k / 2; ++block) {
        tg.run([&, block] { result[block] = bipartion_block(partitioned_hg, context, block); });
      }
      tg.wait();
      for (PartitionID block = 0; block < context.partition.k / 2; ++block) {
        result[block].partitioned_hg.setHypergraph(result[block].hg);
      }

      // We now assign the nodes of block b to the blocks 2*b and 2*b + 1 based on
      // the bipartition of block b. Note that when k is odd, there exists one block
      // that we do not bipartition, which we then assign to block k - 1
      const PartitionID non_bipartitioned_block = context.partition.k % 2 == 1 ?
        context.partition.k / 2 : kInvalidPartition;
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID from = partitioned_hg.partID(hn);
        PartitionID to = context.partition.k - 1;
        if ( from != non_bipartitioned_block ) {
          to = result[from].partitioned_hg.partID(result[from].mapping[hn]) == 0 ? 2 * from : 2 * from + 1;
        }
        if ( from != to ) {
          partitioned_hg.changeNodePart(hn, from, to);
        }
      });
      DBG << "Increased number of blocks from" << best.partitioned_hypergraph.k()
          << "to" << partitioned_hg.k() << "and objective from" << best.objective
          << "to" << metrics::objective(partitioned_hg, context.partition.objective);
    }
  }
}
}

void tmp::recursively_perform_multilevel_partitioning(PartitionedHypergraph& partitioned_hg,
                                                      DeepPartitionResult& result,
                                                      const OriginalHypergraphInfo& info,
                                                      DeepMultilevelOutput& output,
                                                      const size_t num_threads,
                                                      const bool is_parallel_recursion,
                                                      const bool is_top_level) {
  result.context = setupMultilevelContext(partitioned_hg.hypergraph(),
    result.context, info, num_threads, is_parallel_recursion ? 0.5 : 1.0, is_top_level);
  result.hypergraph = partitioned_hg.hypergraph().copy(parallel_tag_t());
  const Context& context = result.context;
  Hypergraph& hypergraph = result.hypergraph;

  DBG << "Perform deep multilevel partitioning call with"
      << "k =" << result.context.partition.k << ","
      << "p =" << result.context.shared_memory.num_threads << ","
      << "c =" << result.context.coarsening.contraction_limit << "and"
      << "rep =" << result.context.initial_partitioning.runs;

  output.addPartitioningCall(context.partition.k,
    context.partition.epsilon, context.coarsening.contraction_limit);

  bool nlevel = context.coarsening.algorithm == CoarseningAlgorithm::nlevel_coarsener;
  UncoarseningData uncoarseningData(nlevel, hypergraph, context);

  // ################## COARSENING ##################
  {
    mt_kahypar::io::printCoarseningBanner(context);
    std::unique_ptr<ICoarsener> coarsener = CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context, uncoarseningData);
    coarsener->coarsen();
  }

  // ################## DEEP MULTILEVEL PARTITIONING ##################
  mt_kahypar::io::printDeepMultilevelBanner(context);
  deep_multilevel_partitioning(uncoarseningData.coarsestPartitionedHypergraph(),
    context, info, output, context.coarsening.contraction_limit, false);

  // ################## UNCOARSENING ##################
  io::printLocalSearchBanner(context);
  std::unique_ptr<IUncoarsener> uncoarsener;
  if ( uncoarseningData.nlevel ) {
    uncoarsener = std::make_unique<NLevelUncoarsener>(hypergraph, context, uncoarseningData);
  } else {
    uncoarsener = std::make_unique<MultilevelUncoarsener>(hypergraph, context, uncoarseningData);
  }
  result.partitioned_hypergraph = uncoarsener->uncoarsen();

  // Compute metrics
  result.objective = metrics::objective(result.partitioned_hypergraph, context.partition.objective);
  result.imbalance = metrics::imbalance(result.partitioned_hypergraph, context);

  output.addPartitioningResult(context.coarsening.contraction_limit,
    result.partitioned_hypergraph, result.objective, result.imbalance);
}

PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context) {
  PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph, parallel_tag_t());
  partition(partitioned_hypergraph, context);
  return partitioned_hypergraph;
}

void partition(PartitionedHypergraph& hypergraph, const Context& context) {
  utils::Utilities& utils = utils::Utilities::instance();
  if (context.partition.mode == Mode::deep_multilevel) {
    utils.getTimer(context.utility_id).start_timer("deep", "Deep Multilevel");
  }
  if (context.type == ContextType::main) {
    parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
    utils.getTimer(context.utility_id).disable();
    utils.getStats(context.utility_id).disable();
  }

  // ################## DEEP MULTILEVEL PARTITIONING ##################
  const bool is_top_level = context.partition.mode == Mode::deep_multilevel;
  const HypernodeID current_num_nodes = !is_top_level ?
      context.coarsening.contraction_limit : std::numeric_limits<HypernodeID>::max();
  tmp::DeepMultilevelOutput output(context);
  tmp::deep_multilevel_partitioning(hypergraph, context,
    OriginalHypergraphInfo { context.partition.k, context.partition.epsilon },
    output, current_num_nodes, is_top_level);

  // ################## V-CYCLES ##################
  if (context.partition.num_vcycles > 0 && context.type == ContextType::main) {
    multilevel::partitionVCycle(hypergraph.hypergraph(), hypergraph, context);
  }

  if (context.type == ContextType::main) {
    parallel::MemoryPool::instance().activate_unused_memory_allocations();
    utils.getTimer(context.utility_id).enable();
    utils.getStats(context.utility_id).enable();
  }
  if (context.partition.mode == Mode::deep_multilevel) {
    utils.getTimer(context.utility_id).stop_timer("deep");
  }
}

} // namespace deep_multilevel
} // namepace mt_kahypar
