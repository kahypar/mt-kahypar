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

#include "mt-kahypar/application/command_line_options.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/mt_kahypar.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/multilevel.h"

using ::testing::Test;

namespace mt_kahypar {


template <template<typename> class InitialPartitioner, InitialPartitioningMode mode, PartitionID k>
struct TestConfig {
  using TypeTraits = GlobalTypeTraits;
  using Partitioner = InitialPartitioner<TypeTraits>;
  static constexpr InitialPartitioningMode MODE = mode;
  static constexpr PartitionID K = k;
};

template <typename Config>
class AInitialPartitionerTest : public Test {
  using TypeTraits = typename Config::TypeTraits;
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using InitialPartitioner = typename Config::Partitioner;

  static size_t num_threads;

 public:
  AInitialPartitionerTest() :
    hypergraph(),
    context() {

    parseIniToContext(context, "../../../../config/fast_preset.ini");

    context.partition.graph_filename = "../test_instances/unweighted_ibm01.hgr";
    context.partition.graph_community_filename = "../test_instances/ibm01.hgr.community";
    context.partition.mode = kahypar::Mode::direct_kway;
    context.partition.objective = kahypar::Objective::km1;
    context.partition.epsilon = 0.2;
    context.partition.k = Config::K;
    context.partition.verbose_output = true;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Community Assignment Strategy
    context.preprocessing.community_redistribution.assignment_objective =
      CommunityAssignmentObjective::pin_objective;
    context.preprocessing.community_redistribution.assignment_strategy =
      CommunityAssignmentStrategy::bin_packing;

    // Initial Partitioning
    context.initial_partitioning.runs = 1;
    context.initial_partitioning.use_adaptive_epsilon = false;
    context.initial_partitioning.mode = Config::MODE;

    // Label Propagation
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;

    // Read hypergraph
    hypergraph = io::readHypergraphFile<HyperGraph, HyperGraphFactory>(
      "../test_instances/unweighted_ibm01.hgr", TBB::GLOBAL_TASK_GROUP);
    partitioned_hypergraph = PartitionedHyperGraph(
      context.partition.k, TBB::GLOBAL_TASK_GROUP, hypergraph);
    context.setupPartWeights(hypergraph.totalWeight());
    context.setupContractionLimit(hypergraph.totalWeight());
    assignCommunities();

    initial_partitioner = std::make_unique<InitialPartitioner>(
      partitioned_hypergraph, context, true, TBB::GLOBAL_TASK_GROUP);
  }

  void assignCommunities() {
    std::vector<PartitionID> communities;
    io::readPartitionFile(context.partition.graph_community_filename, communities);

    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, communities[hypergraph.originalNodeID(hn)]);
    }
    hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

    std::unique_ptr<preprocessing::ICommunityAssignment> community_assignment =
      RedistributionFactory::getInstance().createObject(
        context.preprocessing.community_redistribution.assignment_strategy, hypergraph, context);
    parallel::scalable_vector<PartitionID> community_node_mapping =
      community_assignment->computeAssignment();
    hypergraph.setCommunityNodeMapping(std::move(community_node_mapping));
  }

  static void SetUpTestSuite() {
    TBB::instance(num_threads);
  }

  HyperGraph hypergraph;
  PartitionedHyperGraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<InitialPartitioner> initial_partitioner;
};

template <typename Config>
size_t AInitialPartitionerTest<Config>::num_threads = HwTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<RecursiveInitialPartitionerT, InitialPartitioningMode::recursive, 2>,
                         TestConfig<RecursiveInitialPartitionerT, InitialPartitioningMode::recursive, 3>,
                         TestConfig<RecursiveInitialPartitionerT, InitialPartitioningMode::recursive, 4>,
                         TestConfig<RecursiveInitialPartitionerT, InitialPartitioningMode::recursive, 5>,
                         TestConfig<RecursiveBisectionInitialPartitionerT, InitialPartitioningMode::recursive_bisection, 2>,
                         TestConfig<RecursiveBisectionInitialPartitionerT, InitialPartitioningMode::recursive_bisection, 3>,
                         TestConfig<RecursiveBisectionInitialPartitionerT, InitialPartitioningMode::recursive_bisection, 4>,
                         TestConfig<RecursiveBisectionInitialPartitionerT, InitialPartitioningMode::recursive_bisection, 5> > TestConfigs;

TYPED_TEST_CASE(AInitialPartitionerTest, TestConfigs);

TYPED_TEST(AInitialPartitionerTest, VerifiesThatAllPartsAreNonEmpty) {
  this->initial_partitioner->initialPartition();

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_GT(this->partitioned_hypergraph.partSize(part_id), 0);
  }
}

TYPED_TEST(AInitialPartitionerTest, VerifiesThatPartSizesAndWeightsAreCorrect) {
  this->initial_partitioner->initialPartition();

  std::vector<HypernodeID> part_size(this->context.partition.k, 0);
  std::vector<HypernodeWeight> part_weight(this->context.partition.k, 0);
  for ( const HypernodeID& hn : this->hypergraph.nodes() ) {
    PartitionID part_id = this->partitioned_hypergraph.partID(hn);
    ++part_size[part_id];
    part_weight[part_id] += this->partitioned_hypergraph.nodeWeight(hn);
  }

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_EQ(this->partitioned_hypergraph.partSize(part_id), part_size[part_id]);
    ASSERT_EQ(this->partitioned_hypergraph.partWeight(part_id), part_weight[part_id]);
  }
}

TYPED_TEST(AInitialPartitionerTest, VerifiesThatAllPartWeightsAreSmallerThanMaxPartWeight) {
  this->initial_partitioner->initialPartition();

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_LE(this->partitioned_hypergraph.partWeight(part_id),
              this->context.partition.max_part_weights[part_id]);
  }
}

}  // namespace mt_kahypar
