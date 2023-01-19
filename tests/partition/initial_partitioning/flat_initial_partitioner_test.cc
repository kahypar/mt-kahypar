/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#include "gmock/gmock.h"

#include <atomic>

#include "tbb/parallel_invoke.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/greedy_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/label_propagation_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/policies/gain_computation_policy.h"
#include "mt-kahypar/partition/initial_partitioning/policies/pq_selection_policy.h"

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

template<typename Config>
class AFlatInitialPartitionerTest : public Test {

 public:
  using InitialPartitioner = typename Config::InitialPartitionerTask;

  AFlatInitialPartitionerTest() :
    hypergraph(),
    partitioned_hypergraph(),
    context(),
    ip_data(nullptr) {
    context.partition.k = Config::K;
    context.partition.epsilon = 0.2;
    context.partition.objective = Objective::km1;
    context.initial_partitioning.lp_initial_block_size = 5;
    context.initial_partitioning.lp_maximum_iterations = 100;
    hypergraph = io::readHypergraphFile(
      "../tests/instances/test_instance.hgr");
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    context.initial_partitioning.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    utils::Utilities::instance().getTimer(context.utility_id).disable();
    ip_data = std::make_unique<InitialPartitioningDataContainer>(partitioned_hypergraph, context);
  }

  void execute() {
    tbb::task_group tg;
    const int seed = 420;
    for ( size_t i = 0; i < Config::RUNS; ++i ) {
      tg.run([&, i] {
        InitialPartitioner ip(Config::ALGORITHM, *ip_data.get(), context, seed + i, i);
        ip.partition();
      });
    }
    tg.wait();
    ip_data->apply();
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<InitialPartitioningDataContainer> ip_data;
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
