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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/registries/register_refinement_algorithms.cpp"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/label_propagation/async_lp_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/utils/randomize.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

using ::testing::Test;
using ::testing::Return;
using ::testing::_;

namespace mt_kahypar {
template <PartitionID k, kahypar::Objective objective>
struct TestConfig { };

template <PartitionID k>
struct TestConfig<k, kahypar::Objective::km1> {
  using Refiner = AsyncLPKm1Refiner;
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = kahypar::Objective::km1;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_km1;
};

template <PartitionID k>
struct TestConfig<k, kahypar::Objective::cut> {
  using Refiner = AsyncLPCutRefiner;
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = kahypar::Objective::cut;
  static constexpr LabelPropagationAlgorithm LP_ALGO = LabelPropagationAlgorithm::label_propagation_cut;
};

template <typename Config>
class AAsynchLPRefiner : public Test {
  using Refiner = typename Config::Refiner;

  static size_t num_threads;

 public:
    AAsynchLPRefiner() :
    hypergraph(),
    partitioned_hypergraph(),
    context(),
    refiner(nullptr),
    metrics(),
    refinement_nodes(),
    lock_manager(),
    contraction_group_id(0) {
    context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = kahypar::Mode::direct_kway;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.epsilon = 0.25;
    context.partition.k = Config::K;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Initial Partitioning
    context.initial_partitioning.mode = InitialPartitioningMode::recursive;
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

    lock_manager = std::make_unique<ds::GroupLockManager>(partitioned_hypergraph.hypergraph().initialNumNodes(),ds::invalidGroupID);
    initialPartition();

    node_anti_duplicator = std::make_unique<ds::ThreadSafeFlagArray<HypernodeID>>(hypergraph.initialNumNodes());
    edge_anti_duplicator = std::make_unique<ds::ThreadSafeFlagArray<HyperedgeID>>(hypergraph.initialNumEdges());
    fm_node_tracker_dummy = std::make_unique<AsyncNodeTracker>();
    fm_node_tracker_dummy->resize(hypergraph.initialNumNodes());

    refiner = AsyncLPRefinerFactory::getInstance().createObject(
            context.refinement.label_propagation.algorithm,
            partitioned_hypergraph.hypergraph(),
            context,
            lock_manager.get(),
            node_anti_duplicator.get(),
            edge_anti_duplicator.get(),
            fm_node_tracker_dummy.get()
            );
  }

  void initialPartition() {
    Context ip_context(context);
    ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    InitialPartitioningDataContainer ip_data(partitioned_hypergraph, ip_context);
    BFSInitialPartitioner& initial_partitioner = *new(tbb::task::allocate_root())
      BFSInitialPartitioner(InitialPartitioningAlgorithm::bfs, ip_data, ip_context, 420, 0);
    tbb::task::spawn_root_and_wait(initial_partitioner);
    ip_data.apply();
    metrics.update_km1_strong(metrics::km1(partitioned_hypergraph));
    metrics.update_cut_strong(metrics::hyperedgeCut(partitioned_hypergraph));
    metrics.update_imbalance_strong(metrics::imbalance(partitioned_hypergraph, context));

    // Set refinement nodes to all border nodes (as the AsyncLPRefiner needs input seed nodes)
//    auto add_border_node_to_refinement_nodes = [&](const HypernodeID& u) {
//        if (partitioned_hypergraph.isBorderNode(u)) {
//            refinement_nodes.push_back(u);
//        }
//    };
//    partitioned_hypergraph.doParallelForAllNodes(add_border_node_to_refinement_nodes);

    for (auto u : partitioned_hypergraph.nodes()) {
        if (partitioned_hypergraph.isBorderNode(u)) {
            refinement_nodes.push_back(u);
        }
    }

  }

  ~AAsynchLPRefiner() {
    refiner.reset();
    node_anti_duplicator.reset();
    edge_anti_duplicator.reset();
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<IAsyncRefiner> refiner;
  metrics::ThreadSafeMetrics metrics;
  parallel::scalable_vector<HypernodeID> refinement_nodes;

  std::unique_ptr<ds::GroupLockManager> lock_manager;
  ds::ContractionGroupID contraction_group_id;
  std::unique_ptr<ds::ThreadSafeFlagArray<HypernodeID>> node_anti_duplicator;
  std::unique_ptr<ds::ThreadSafeFlagArray<HyperedgeID>> edge_anti_duplicator;
  std::unique_ptr<AsyncNodeTracker> fm_node_tracker_dummy;
};

template <typename Config>
size_t AAsynchLPRefiner<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<2, kahypar::Objective::cut>,
                         TestConfig<4, kahypar::Objective::cut>,
                         TestConfig<8, kahypar::Objective::cut>,
                         TestConfig<2, kahypar::Objective::km1>,
                         TestConfig<4, kahypar::Objective::km1>,
                         TestConfig<8, kahypar::Objective::km1> > TestConfigs;

TYPED_TEST_CASE(AAsynchLPRefiner, TestConfigs);

TYPED_TEST(AAsynchLPRefiner, UpdatesImbalanceCorrectly) {
  this->refiner->refine(this->partitioned_hypergraph, this->refinement_nodes, this->metrics, std::numeric_limits<double>::max(),this->contraction_group_id);
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.loadImbalance());
}

TYPED_TEST(AAsynchLPRefiner, DoesNotViolateBalanceConstraint) {
  this->refiner->refine(this->partitioned_hypergraph, this->refinement_nodes, this->metrics, std::numeric_limits<double>::max(),this->contraction_group_id);
  ASSERT_LE(this->metrics.loadImbalance(), this->context.partition.epsilon + EPS);
}

TYPED_TEST(AAsynchLPRefiner, UpdatesMetricsCorrectly) {
  this->refiner->refine(this->partitioned_hypergraph, this->refinement_nodes, this->metrics, std::numeric_limits<double>::max(),this->contraction_group_id);
  ASSERT_EQ(metrics::objective(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.getMetric(kahypar::Mode::direct_kway, this->context.partition.objective));
}

TYPED_TEST(AAsynchLPRefiner, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::objective(this->partitioned_hypergraph, this->context.partition.objective);
  this->refiner->refine(this->partitioned_hypergraph, this->refinement_nodes, this->metrics, std::numeric_limits<double>::max(),this->contraction_group_id);
  ASSERT_LE(this->metrics.getMetric(kahypar::Mode::direct_kway, this->context.partition.objective), objective_before);
}

TYPED_TEST(AAsynchLPRefiner, RefineWithEmptyRefinementNodesDeathTest) {
    testing::FLAGS_gtest_death_test_style="threadsafe";
    ASSERT_DEATH(this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max(),this->contraction_group_id), "");
}

}  // namespace mt_kahypar
