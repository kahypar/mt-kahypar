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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/registries/register_initial_partitioning_algorithms.h"
#include "mt-kahypar/partition/registries/register_refinement_algorithms.cpp"
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"

using ::testing::Test;

namespace mt_kahypar {

template<PartitionID k, size_t runs>
struct TestConfig {
  static constexpr PartitionID K = k;
  static constexpr size_t RUNS = runs;
};

template<typename Config>
class APoolInitialPartitionerTest : public Test {
 public:
  APoolInitialPartitionerTest() :
    hypergraph(),
    partitioned_hypergraph(),
    context() {
    context.partition.k = Config::K;
    context.partition.epsilon = 0.2;
    context.partition.objective = Objective::km1;
    context.initial_partitioning.runs = Config::RUNS;
    context.refinement.label_propagation.algorithm =
      LabelPropagationAlgorithm::label_propagation_km1;
    context.initial_partitioning.refinement.label_propagation.algorithm =
      LabelPropagationAlgorithm::label_propagation_km1;
    hypergraph = io::readHypergraphFile(
      "../tests/instances/test_instance.hgr");
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    utils::Utilities::instance().getTimer(context.utility_id).disable();
  }

  static void SetUpTestSuite() {
    TBBInitializer::instance(HardwareTopology::instance().num_cpus());
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
};

typedef ::testing::Types<TestConfig<2, 1>,
                         TestConfig<2, 2>,
                         TestConfig<2, 5>,
                         TestConfig<3, 1>,
                         TestConfig<3, 2>,
                         TestConfig<3, 5>,
                         TestConfig<4, 1>,
                         TestConfig<4, 2>,
                         TestConfig<4, 5>,
                         TestConfig<5, 1>,
                         TestConfig<5, 2>,
                         TestConfig<5, 5>> TestConfigs;

TYPED_TEST_CASE(APoolInitialPartitionerTest, TestConfigs);

TYPED_TEST(APoolInitialPartitionerTest, HasValidImbalance) {
  pool::bipartition(this->partitioned_hypergraph, this->context);

  ASSERT_LE(metrics::imbalance(this->partitioned_hypergraph, this->context),
            this->context.partition.epsilon);
}

TYPED_TEST(APoolInitialPartitionerTest, AssginsEachHypernode) {
  pool::bipartition(this->partitioned_hypergraph, this->context);

  for ( const HypernodeID& hn : this->partitioned_hypergraph.nodes() ) {
    ASSERT_NE(this->partitioned_hypergraph.partID(hn), -1);
  }
}

TYPED_TEST(APoolInitialPartitionerTest, HasNoSignificantLowPartitionWeights) {
  pool::bipartition(this->partitioned_hypergraph, this->context);

  // Each block should have a weight greater or equal than 20% of the average
  // block weight.
  for ( PartitionID block = 0; block < this->context.partition.k; ++block ) {
    ASSERT_GE(this->partitioned_hypergraph.partWeight(block),
              this->context.partition.perfect_balance_part_weights[block] / 5);
  }
}


}  // namespace mt_kahypar
