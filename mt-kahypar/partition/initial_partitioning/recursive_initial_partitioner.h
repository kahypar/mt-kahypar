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

#include <future>

#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <libkahypar.h>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"

namespace mt_kahypar {

namespace multilevel {
  static inline void partition(Hypergraph& hypergraph, const Context& context, const bool top_level);
} // namespace multilevel


template< typename TypeTraits >
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

 public:
  RecursiveInitialPartitionerT(HyperGraph& hypergraph, const Context& context, const bool top_level) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level) { }

  RecursiveInitialPartitionerT(const RecursiveInitialPartitionerT&) = delete;
  RecursiveInitialPartitionerT(RecursiveInitialPartitionerT&&) = delete;
  RecursiveInitialPartitionerT& operator= (const RecursiveInitialPartitionerT&) = delete;
  RecursiveInitialPartitionerT& operator= (RecursiveInitialPartitionerT&&) = delete;

 private:
  void initialPartitionImpl() override final {
    if ( _context.shared_memory.num_threads == 1 &&
         _context.coarsening.contraction_limit == 2 * _context.coarsening.contraction_limit_multiplier ) {
      // Base Case -> Bisect Hypergraph
      initialBisection();
      return;
    }

    bool do_parallel_recursion = _context.coarsening.contraction_limit /
                                 ( 2 * _context.coarsening.contraction_limit_multiplier ) ==
                                 _context.shared_memory.num_threads;

    RecursivePartitionResult r1;
    RecursivePartitionResult r2;
    if ( do_parallel_recursion ) {
      // Perform parallel recursion
      size_t num_threads_1 = std::ceil( ( (double) _context.shared_memory.num_threads ) / 2.0 );
      size_t num_threads_2 = std::floor( ( (double) _context.shared_memory.num_threads ) / 2.0 );
      std::future<RecursivePartitionResult> f1 = std::async(std::launch::async, [&] { return recursivePartition(num_threads_1); } );
      std::future<RecursivePartitionResult> f2 = std::async(std::launch::async, [&] { return recursivePartition(num_threads_2); } );
      r1 = f1.get();
      r2 = f2.get();
      /*r1 = recursivePartition(num_threads_1);
      r2 = recursivePartition(num_threads_2);*/
    } else {
      r1 = recursivePartition(_context.shared_memory.num_threads);
      r2.objective = std::numeric_limits<HyperedgeWeight>::max();
      r2.imbalance = std::numeric_limits<double>::max();
    }

    // Choose best partition of both parallel recursion
    RecursivePartitionResult best;
    bool r1_has_better_quality = r1.objective < r2.objective;
    bool r1_is_balanced = r1.imbalance < r1.context.partition.epsilon;
    bool r2_is_balanced = r2.imbalance < r2.context.partition.epsilon;
    if ( ( r1_has_better_quality && r1_is_balanced ) ||
         ( r1_is_balanced && !r2_is_balanced ) ||
         ( r1_has_better_quality && !r1_is_balanced && !r2_is_balanced ) ) {
      best = std::move(r1);
    } else {
      best = std::move(r2);
    }

    HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective == metrics::objective(best.hypergraph, _context.partition.objective));

    // Apply best partition to hypergraph
    for ( const HypernodeID& hn : _hg.nodes() ) {
      const HypernodeID original_id = _hg.originalNodeID(hn);
      ASSERT(original_id < best.mapping.size());
      PartitionID part_id = best.hypergraph.partID(best.hypergraph.globalNodeID(best.mapping[original_id]));
      ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
      _hg.setNodePart(hn, part_id);
    }
    _hg.initializeNumCutHyperedges();
    _hg.updateGlobalPartInfos();


    HEAVY_INITIAL_PARTITIONING_ASSERT(best.objective == metrics::objective(_hg, _context.partition.objective));

    // Bisect all blocks of best partition
    bool perform_bisections = !_top_level && _context.shared_memory.num_threads < (size_t) _context.partition.k;
    if ( perform_bisections ) {
      std::vector<HypernodeID> node_mapping(_hg.initialNumNodes(), kInvalidHypernode);
      std::vector<KaHyParParitioningResult> results;
      results.reserve(_context.partition.k / 2);

      // Bisect all blocks of the current partition in parallel
      tbb::task_arena arena(_context.shared_memory.num_threads, 0);
      tbb::task_group group;
      for ( PartitionID k = 0; k < _context.partition.k / 2; ++k ) {
        results.emplace_back(_hg.partSize(k));
        arena.execute([&, k] {
          group.run([&, k] {
            bisectPartition(k, results[k], node_mapping);
          });
        });
      }
      arena.execute([&] { group.wait(); });

      // Apply all bisections to current hypergraph
      PartitionID unbisected_block = ( _context.partition.k % 2 == 1 ? (PartitionID) results.size() : kInvalidPartition );
      for ( const HypernodeID& hn : _hg.nodes() ) {
        PartitionID from = _hg.partID(hn);
        PartitionID to = kInvalidPartition;
        if ( from != unbisected_block ) {
          ASSERT(from != kInvalidPartition && from < (PartitionID) results.size());
          ASSERT(node_mapping[_hg.originalNodeID(hn)] < results[from].partition.size());
          to = results[from].partition[node_mapping[_hg.originalNodeID(hn)]] == 0 ? 2 * from : 2 * from + 1;
        } else {
          to = _context.partition.k - 1;
        }

        ASSERT(to != kInvalidPartition && to < _hg.k());
        if ( from != to ) {
          _hg.changeNodePart(hn, from, to);
        }
      }

      HEAVY_INITIAL_PARTITIONING_ASSERT([&] {
        HyperedgeWeight expected_objective = best.objective;
        HyperedgeWeight actual_objective = metrics::objective(_hg, _context.partition.objective);
        for ( size_t i = 0; i < results.size(); ++i ) {
          expected_objective += results[i].objective;
        }

        if ( expected_objective != actual_objective ) {
          LOG << V(expected_objective) << V(actual_objective);
          return false;
        }
        return true;
      }());
    }

    _hg.updateGlobalPartInfos();
  }

  RecursivePartitionResult recursivePartition(const size_t num_threads) {
    RecursivePartitionResult result(setupRecursiveContext(num_threads));

    // Copy hypergraph
    auto copy = _hg.copy(result.context.partition.k);
    result.hypergraph = std::move(copy.first);
    result.mapping = std::move(copy.second);

    // Call multilevel partitioner recursively
    DBG << "Perform recursive multilevel partitioner call with"
        << "k =" << result.context.partition.k << ","
        << "p =" << result.context.shared_memory.num_threads << ","
        << "c =" << result.context.coarsening.contraction_limit << "and"
        << "rep =" << result.context.initial_partitioning.runs;
    multilevel::partition(result.hypergraph, result.context, false);

    result.objective = metrics::objective(result.hypergraph, result.context.partition.objective);
    result.imbalance = metrics::imbalance(result.hypergraph, result.context);
    return result;
  }



  void initialBisection() {
    ASSERT(_context.partition.k == 2);
    ASSERT(_context.partition.max_part_weights.size() == 2);
    ASSERT(_context.shared_memory.num_threads == 1);

    kahypar_context_t* kahypar_context = setupKaHyParContext();
    kahypar_set_custom_target_block_weights(2, _context.partition.max_part_weights.data(), kahypar_context);

    // Compute node mapping that maps all hypernodes of the coarsened hypergraph
    // to a consecutive range of node ids
    HypernodeID num_vertices = 0;
    std::vector<HypernodeID> node_mapping(_hg.initialNumNodes(), kInvalidHypernode);
    for ( const HypernodeID& hn : _hg.nodes() ) {
      ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
      node_mapping[_hg.originalNodeID(hn)] = num_vertices++;
    }

    // Build KaHyPar Hypergraph Data Structure
    const KaHyParHypergraph kahypar_hypergraph = convertToKaHyParHypergraph(_hg, node_mapping);
    KaHyParParitioningResult result(num_vertices);

    // Bisect KaHyPar Hypergraph
    kahypar_partition(kahypar_hypergraph, *reinterpret_cast<kahypar::Context*>(kahypar_context), result);

    // Apply partition to MT-KaHyPar Hypergraph
    for ( HypernodeID u = 0; u < result.partition.size(); ++u ) {
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

  void bisectPartition(const PartitionID k,
                       KaHyParParitioningResult& result,
                       std::vector<HypernodeID>& node_mapping) {

    kahypar_context_t* kahypar_context = setupKaHyParContext();
    std::vector<HypernodeWeight> max_part_weights;
    max_part_weights.emplace_back(_context.partition.max_part_weights[2 * k]);
    max_part_weights.emplace_back(_context.partition.max_part_weights[2 * k + 1]);
    kahypar_set_custom_target_block_weights(2, max_part_weights.data(), kahypar_context);

    // Build KaHyPar Hypergraph Data Structure
    bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
    KaHyParHypergraph kahypar_hypergraph = extractBlockAsKaHyParHypergraph(_hg, k, node_mapping, cut_net_splitting);

    // Bisect KaHyPar Hypergraph
    kahypar_partition(kahypar_hypergraph, *reinterpret_cast<kahypar::Context*>(kahypar_context), result);

    kahypar_context_free(kahypar_context);
  }

  Context setupRecursiveContext(const size_t num_threads) {
    ASSERT(num_threads >= 1);
    Context context(_context);
    context.shared_memory.num_threads = num_threads;

    bool reduce_k = !_top_level && _context.shared_memory.num_threads < (size_t) _context.partition.k && _context.partition.k > 2;
    if ( reduce_k ) {
      context.partition.k = std::ceil( ( (double) context.partition.k ) / 2.0 );

      context.partition.perfect_balance_part_weights.clear();
      context.partition.max_part_weights.clear();
      for ( PartitionID part = 0; part < _context.partition.k / 2; ++part ) {
        context.partition.perfect_balance_part_weights.emplace_back(
          _context.partition.perfect_balance_part_weights[2 * part] +
          _context.partition.perfect_balance_part_weights[2 * part + 1]);
        context.partition.max_part_weights.emplace_back(
          _context.partition.max_part_weights[2 * part] +
          _context.partition.max_part_weights[2 * part + 1]);
      }
      if ( _context.partition.k % 2 == 1 ) {
        context.partition.perfect_balance_part_weights.emplace_back(
          _context.partition.perfect_balance_part_weights.back());
        context.partition.max_part_weights.emplace_back(
          _context.partition.max_part_weights.back());
      }
      ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t) context.partition.k);
      ASSERT(context.partition.max_part_weights.size() == (size_t) context.partition.k);
    }
    context.partition.verbose_output = debug;

    context.coarsening.contraction_limit = std::max( context.partition.k * context.coarsening.contraction_limit_multiplier,
                                                     2 * context.shared_memory.num_threads * context.coarsening.contraction_limit_multiplier );
    context.coarsening.hypernode_weight_fraction = context.coarsening.max_allowed_weight_multiplier /
                                                   context.coarsening.contraction_limit;
    context.coarsening.max_allowed_node_weight = ceil(context.coarsening.hypernode_weight_fraction * _hg.totalWeight());

    bool is_parallel_recursion = _context.shared_memory.num_threads != context.shared_memory.num_threads;
    context.initial_partitioning.runs = std::max( context.initial_partitioning.runs / ( is_parallel_recursion ? 2 : 1 ), 1UL );

    return context;
  }

  kahypar_context_t* setupKaHyParContext() {
    kahypar_context_t* kahypar_c = mt_kahypar::setupContext(_context, kahypar_debug);
    kahypar::Context& kahypar_context = *reinterpret_cast<kahypar::Context*>(kahypar_c);
    kahypar_context.partition.k = 2;
    kahypar_context.coarsening.algorithm = kahypar::CoarseningAlgorithm::do_nothing;
    kahypar_context.local_search.algorithm = kahypar::RefinementAlgorithm::do_nothing;
    kahypar_context.initial_partitioning.nruns = _context.initial_partitioning.runs;
    kahypar_context.partition.seed = utils::Randomize::instance().getRandomInt(0, 2000, sched_getcpu());
    return kahypar_c;
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
};

template< typename TypeTraits >
PartitionID RecursiveInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template< typename TypeTraits >
HypernodeID RecursiveInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();


using RecursiveInitialPartitioner = RecursiveInitialPartitionerT<GlobalTypeTraits>;

}  // namespace mt_kahypar
