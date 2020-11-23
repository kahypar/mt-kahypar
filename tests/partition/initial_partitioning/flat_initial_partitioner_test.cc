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

#include "gmock/gmock.h"

#include <atomic>

#include "tbb/parallel_invoke.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/initial_partitioning/flat/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/greedy_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/label_propagation_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/gain_computation_policy.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/pq_selection_policy.h"

using ::testing::Test;

namespace mt_kahypar {

template<class InitialPartitioner,
         InitialPartitioningAlgorithm algorithm,
         PartitionID k, size_t runs>
struct TestConfig {
  using InitialPartitionerTask = InitialPartitioner;
  static constexpr InitialPartitioningAlgorithm ALGORITHM = algorithm;
  static constexpr PartitionID K = k;
  static constexpr size_t RUNS = runs;
};

template<class InitialPartitionerTask>
class InitialPartitionerRootTaskT : public tbb::task {

 public:
  InitialPartitionerRootTaskT(PartitionedHypergraph& hypergraph,
                              const Context& context,
                              const InitialPartitioningAlgorithm algorithm,
                              const size_t runs) :
    _ip_data(hypergraph, context, TBBNumaArena::GLOBAL_TASK_GROUP),
    _context(context),
    _algorithm(algorithm),
    _runs(runs) {}

  tbb::task* execute() override {
    tbb::task::set_ref_count(_runs + 1);
    const int seed = 420;
    for ( size_t i = 0; i < _runs; ++i ) {
      tbb::task::spawn(*new(tbb::task::allocate_child())
        InitialPartitionerTask(_algorithm, _ip_data, _context, seed + i));
    }
    tbb::task::wait_for_all();
    _ip_data.apply();
    return nullptr;
  }

 private:
  InitialPartitioningDataContainer _ip_data;
  const Context& _context;
  const InitialPartitioningAlgorithm _algorithm;
  const size_t _runs;
};

template<typename Config>
class AFlatInitialPartitionerTest : public Test {

 public:
  using InitialPartitionerRootTask = InitialPartitionerRootTaskT<typename Config::InitialPartitionerTask>;

  AFlatInitialPartitionerTest() :
    hypergraph(),
    partitioned_hypergraph(),
    context() {
    context.partition.k = Config::K;
    context.partition.epsilon = 0.2;
    context.partition.objective = kahypar::Objective::km1;
    context.initial_partitioning.lp_initial_block_size = 5;
    context.initial_partitioning.lp_maximum_iterations = 100;
    hypergraph = io::readHypergraphFile(
      "../tests/instances/test_instance.hgr", TBBNumaArena::GLOBAL_TASK_GROUP);
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, TBBNumaArena::GLOBAL_TASK_GROUP, hypergraph);
    context.setupPartWeights(hypergraph.totalWeight());
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    context.initial_partitioning.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    utils::Timer::instance().disable();
  }

  void execute() {
    InitialPartitionerRootTask& root_ip_task = *new(tbb::task::allocate_root())
      InitialPartitionerRootTask(partitioned_hypergraph, context, Config::ALGORITHM, Config::RUNS);
    tbb::task::spawn_root_and_wait(root_ip_task);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
};

using GreedyRoundRobinFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, SequentialPQSelectionPolicy>;
using GreedyRoundRobinMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, SequentialPQSelectionPolicy>;

typedef ::testing::Types<TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 2, 1>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 2, 2>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 2, 5>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 3, 1>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 3, 2>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 3, 5>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 4, 1>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 4, 2>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 4, 5>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 5, 1>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 5, 2>,
                         TestConfig<RandomInitialPartitioner, InitialPartitioningAlgorithm::random, 5, 5>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 2, 1>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 2, 2>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 2, 5>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 3, 1>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 3, 2>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 3, 5>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 4, 1>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 4, 2>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 4, 5>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 5, 1>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 5, 2>,
                         TestConfig<BFSInitialPartitioner, InitialPartitioningAlgorithm::bfs, 5, 5>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 2, 1>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 2, 2>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 2, 5>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 3, 1>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 3, 2>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 3, 5>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 4, 1>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 4, 2>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 4, 5>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 5, 1>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 5, 2>,
                         TestConfig<GreedyRoundRobinFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_fm, 5, 5>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 2, 1>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 2, 2>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 2, 5>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 3, 1>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 3, 2>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 3, 5>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 4, 1>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 4, 2>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 4, 5>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 5, 1>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 5, 2>,
                         TestConfig<GreedyGlobalFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_fm, 5, 5>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 2, 1>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 2, 2>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 2, 5>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 3, 1>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 3, 2>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 3, 5>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 4, 1>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 4, 2>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 4, 5>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 5, 1>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 5, 2>,
                         TestConfig<GreedySequentialFMInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_fm, 5, 5>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 2, 1>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 2, 2>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 2, 5>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 3, 1>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 3, 2>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 3, 5>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 4, 1>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 4, 2>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 4, 5>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 5, 1>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 5, 2>,
                         TestConfig<GreedyRoundRobinMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_round_robin_max_net, 5, 5>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 2, 1>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 2, 2>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 2, 5>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 3, 1>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 3, 2>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 3, 5>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 4, 1>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 4, 2>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 4, 5>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 5, 1>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 5, 2>,
                         TestConfig<GreedyGlobalMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_global_max_net, 5, 5>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 2, 1>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 2, 2>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 2, 5>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 3, 1>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 3, 2>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 3, 5>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 4, 1>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 4, 2>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 4, 5>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 5, 1>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 5, 2>,
                         TestConfig<GreedySequentialMaxNetInitialPartitioner, InitialPartitioningAlgorithm::greedy_sequential_max_net, 5, 5>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 2, 1>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 2, 2>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 2, 5>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 3, 1>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 3, 2>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 3, 5>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 4, 1>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 4, 2>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 4, 5>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 5, 1>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 5, 2>,
                         TestConfig<LabelPropagationInitialPartitioner, InitialPartitioningAlgorithm::label_propagation, 5, 5> > TestConfigs;

TYPED_TEST_CASE(AFlatInitialPartitionerTest, TestConfigs);

TYPED_TEST(AFlatInitialPartitionerTest, HasValidImbalance) {
  this->execute();

  ASSERT_LE(metrics::imbalance(this->partitioned_hypergraph, this->context),
            this->context.partition.epsilon);
}

TYPED_TEST(AFlatInitialPartitionerTest, AssginsEachHypernode) {
  this->execute();

  for ( const HypernodeID& hn : this->hypergraph.nodes() ) {
    ASSERT_NE(this->partitioned_hypergraph.partID(hn), -1);
  }
}

TYPED_TEST(AFlatInitialPartitionerTest, HasNoSignificantLowPartitionWeights) {
  this->execute();

  // Each block should have a weight greater or equal than 20% of the average
  // block weight.
  for ( PartitionID block = 0; block < this->context.partition.k; ++block ) {
    ASSERT_GE(this->partitioned_hypergraph.partWeight(block),
              this->context.partition.perfect_balance_part_weights[block] / 5);
  }
}


}  // namespace mt_kahypar
