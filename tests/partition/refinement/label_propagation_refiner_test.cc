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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/initial_partitioning/direct_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/policies/execution_policy.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/partition/refinement/label_propagation_refiner.h"

#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {

template< PartitionID k, kahypar::Objective objective >
struct TestConfig { };

template< PartitionID k >
struct TestConfig<k, kahypar::Objective::km1> {
  using TypeTraits = ds::TestTypeTraits<2>;
  using Refiner = LabelPropagationRefinerT<TypeTraits, ExponentialExecutionPolicy, Km1Policy>;
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = kahypar::Objective::km1;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_km1;
};

template< PartitionID k >
struct TestConfig<k, kahypar::Objective::cut> {
  using TypeTraits = ds::TestTypeTraits<2>;
  using Refiner = LabelPropagationRefinerT<TypeTraits, ExponentialExecutionPolicy, CutPolicy>;
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = kahypar::Objective::cut;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_cut;
};

template< typename Config >
class ALabelPropagationRefiner : public Test {

using Refiner = typename Config::Refiner;
using TypeTraits = typename Config::TypeTraits;
using HyperGraph = typename TypeTraits::HyperGraph;
using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
using TBB = typename TypeTraits::TBB;
using HwTopology = typename TypeTraits::HwTopology;

static size_t num_threads;

 public:
  ALabelPropagationRefiner() :
    hypergraph(),
    context(),
    refiner(nullptr),
    metrics() {
    context.partition.graph_filename = "test_instances/ibm01.hgr";
    context.partition.graph_community_filename = "test_instances/ibm01.hgr.community";
    context.partition.mode = Mode::direct_kway;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.epsilon = 0.25;
    context.partition.k = Config::K;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;
    context.shared_memory.use_community_redistribution = true;
    context.shared_memory.assignment_strategy = CommunityAssignmentStrategy::bin_packing;
    context.shared_memory.assignment_objective = CommunityAssignmentObjective::pin_objective;

    // Initial Partitioning
    context.initial_partitioning.mode = InitialPartitioningMode::recursive;
    context.initial_partitioning.runs = 1;
    context.initial_partitioning.context_file = "../test_instances/fast_initial_partitioning.ini";

    // Label Propagation
    context.refinement.label_propagation.algorithm = Config::LP_ALGO;
    context.refinement.label_propagation.execution_policy = ExecutionType::exponential;
    context.refinement.label_propagation.part_weight_update_frequency = 5;

    // Read hypergraph
    hypergraph = io::readHypergraphFile<HyperGraph, StreamingHyperGraph, TBB, HwTopology>(
      "../test_instances/unweighted_ibm01.hgr", context.partition.k, InitialHyperedgeDistribution::equally);
    context.setupPartWeights(hypergraph.totalWeight());
    initialPartition();

    refiner = std::make_unique<Refiner>(hypergraph, context);
  }

  static void SetUpTestSuite() {
    TBB::instance(num_threads);
  }

  void initialPartition() {
    DirectInitialPartitionerT<TypeTraits> initial_partitioner(hypergraph, context, true);
    initial_partitioner.initialPartition();
    metrics.km1 = metrics::km1(hypergraph);
    metrics.cut = metrics::hyperedgeCut(hypergraph);
    metrics.imbalance = metrics::imbalance(hypergraph, context);
  }

  HyperGraph hypergraph;
  Context context;
  std::unique_ptr<Refiner> refiner;
  kahypar::Metrics metrics;
};

template< typename Config >
size_t ALabelPropagationRefiner<Config>::num_threads = std::thread::hardware_concurrency();

static constexpr double EPS = 10e-6;

typedef ::testing::Types<TestConfig<2,   kahypar::Objective::cut>,
                         TestConfig<4,   kahypar::Objective::cut>,
                         TestConfig<8,   kahypar::Objective::cut>,
                         TestConfig<2,   kahypar::Objective::km1>,
                         TestConfig<4,   kahypar::Objective::km1>,
                         TestConfig<8,   kahypar::Objective::km1> > TestConfigs;

TYPED_TEST_CASE(ALabelPropagationRefiner, TestConfigs);


TYPED_TEST(ALabelPropagationRefiner, UpdatesImbalanceCorrectly) {
  this->refiner->refine({}, this->metrics);
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->hypergraph, this->context), this->metrics.imbalance);
}

TYPED_TEST(ALabelPropagationRefiner, DoesNotViolateBalanceConstraint) {
  this->refiner->refine({}, this->metrics);
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon + EPS);
}

TYPED_TEST(ALabelPropagationRefiner, UpdatesMetricsCorrectly) {
  this->refiner->refine({}, this->metrics);
  ASSERT_EQ(metrics::objective(this->hypergraph, this->context.partition.objective),
            this->metrics.getMetric(kahypar::Mode::direct_kway, this->context.partition.objective));
}

TYPED_TEST(ALabelPropagationRefiner, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::objective(this->hypergraph, this->context.partition.objective);
  this->refiner->refine({}, this->metrics);
  ASSERT_LE(this->metrics.getMetric(kahypar::Mode::direct_kway, this->context.partition.objective), objective_before);
}



} // namespace mt_kahypar
