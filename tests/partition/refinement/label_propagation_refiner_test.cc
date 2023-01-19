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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/registries/register_refinement_algorithms.cpp"
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {
template <PartitionID k, Objective objective>
struct TestConfig { };

template <PartitionID k>
struct TestConfig<k, Objective::km1> {
  using Refiner = LabelPropagationKm1Refiner;
  static constexpr PartitionID K = k;
  static constexpr Objective OBJECTIVE = Objective::km1;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_km1;
};

template <PartitionID k>
struct TestConfig<k, Objective::cut> {
  using Refiner = LabelPropagationCutRefiner;
  static constexpr PartitionID K = k;
  static constexpr Objective OBJECTIVE = Objective::cut;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_cut;
};

template <typename Config>
class ALabelPropagationRefiner : public Test {
  using Refiner = typename Config::Refiner;

  static size_t num_threads;

 public:
  ALabelPropagationRefiner() :
    hypergraph(),
    partitioned_hypergraph(),
    context(),
    refiner(nullptr),
    metrics() {
    context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = Mode::direct;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.epsilon = 0.25;
    context.partition.k = Config::K;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Initial Partitioning
    context.initial_partitioning.mode = Mode::deep_multilevel;
    context.initial_partitioning.runs = 1;

    // Label Propagation
    context.refinement.label_propagation.algorithm = Config::LP_ALGO;
    context.initial_partitioning.refinement.label_propagation.algorithm = Config::LP_ALGO;

    // Read hypergraph
    hypergraph = io::readHypergraphFile(
      "../tests/instances/contracted_unweighted_ibm01.hgr");
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    initialPartition();

    refiner = std::make_unique<Refiner>(hypergraph, context);
    refiner->initialize(partitioned_hypergraph);
  }

  void initialPartition() {
    Context ip_context(context);
    ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    InitialPartitioningDataContainer ip_data(partitioned_hypergraph, ip_context);
    BFSInitialPartitioner initial_partitioner(InitialPartitioningAlgorithm::bfs, ip_data, ip_context, 420, 0);
    initial_partitioner.partition();
    ip_data.apply();
    metrics.km1 = metrics::km1(partitioned_hypergraph);
    metrics.cut = metrics::hyperedgeCut(partitioned_hypergraph);
    metrics.imbalance = metrics::imbalance(partitioned_hypergraph, context);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<Refiner> refiner;
  Metrics metrics;
};

template <typename Config>
size_t ALabelPropagationRefiner<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<2, Objective::cut>,
                         TestConfig<4, Objective::cut>,
                         TestConfig<8, Objective::cut>,
                         TestConfig<2, Objective::km1>,
                         TestConfig<4, Objective::km1>,
                         TestConfig<8, Objective::km1> > TestConfigs;

TYPED_TEST_CASE(ALabelPropagationRefiner, TestConfigs);

TYPED_TEST(ALabelPropagationRefiner, UpdatesImbalanceCorrectly) {
  this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
}

TYPED_TEST(ALabelPropagationRefiner, DoesNotViolateBalanceConstraint) {
  this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon + EPS);
}

TYPED_TEST(ALabelPropagationRefiner, UpdatesMetricsCorrectly) {
  this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_EQ(metrics::objective(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.getMetric(Mode::direct, this->context.partition.objective));
}

TYPED_TEST(ALabelPropagationRefiner, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::objective(this->partitioned_hypergraph, this->context.partition.objective);
  this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.getMetric(Mode::direct, this->context.partition.objective), objective_before);
}
}  // namespace mt_kahypar
