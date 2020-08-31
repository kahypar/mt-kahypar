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

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/initial_partitioning/recursive_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_bisection_initial_partitioner.h"

using ::testing::Test;

namespace mt_kahypar {


template <class InitialPartitioner, InitialPartitioningMode mode, PartitionID k>
struct TestConfig {
  using Partitioner = InitialPartitioner;
  static constexpr InitialPartitioningMode MODE = mode;
  static constexpr PartitionID K = k;
};

template <typename Config>
class AInitialPartitionerTest : public Test {
  using InitialPartitioner = typename Config::Partitioner;

  static size_t num_threads;

 public:
  AInitialPartitionerTest() :
    hypergraph(),
    context() {

    parseIniToContext(context, "../config/speed_preset.ini");

    context.partition.graph_filename = "../tests/instances/contracted_unweighted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = kahypar::Mode::direct_kway;
    context.partition.objective = kahypar::Objective::km1;
    context.partition.epsilon = 0.2;
    context.partition.k = Config::K;
    context.partition.verbose_output = true;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Coarsening
    if ( context.partition.paradigm == Paradigm::multilevel ) {
      context.coarsening.algorithm = CoarseningAlgorithm::multilevel_coarsener;
    } else {
      context.coarsening.algorithm = CoarseningAlgorithm::nlevel_coarsener;
      context.refinement.max_batch_size = 100;
    }

    // Initial Partitioning
    context.initial_partitioning.runs = 1;
    context.sparsification.use_degree_zero_contractions = false;
    context.sparsification.use_heavy_net_removal = false;
    context.sparsification.use_similiar_net_removal = false;
    context.sparsification.use_degree_zero_contractions = false;
    context.initial_partitioning.use_adaptive_epsilon = false;
    context.initial_partitioning.mode = Config::MODE;
    context.initial_partitioning.remove_degree_zero_hns_before_ip = false;

    // Label Propagation
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    context.initial_partitioning.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;

    // Read hypergraph
    hypergraph = io::readHypergraphFile(
      "../tests/instances/contracted_unweighted_ibm01.hgr", TBBNumaArena::GLOBAL_TASK_GROUP);
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, TBBNumaArena::GLOBAL_TASK_GROUP, hypergraph);
    context.setupPartWeights(hypergraph.totalWeight());
    context.setupContractionLimit(hypergraph.totalWeight());
    assignCommunities();

    initial_partitioner = std::make_unique<InitialPartitioner>(
      partitioned_hypergraph, context, true, TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  void assignCommunities() {
    std::vector<PartitionID> communities;
    io::readPartitionFile(context.partition.graph_community_filename, communities);

    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, communities[hn]);
    }
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<InitialPartitioner> initial_partitioner;
};

template <typename Config>
size_t AInitialPartitionerTest<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<RecursiveInitialPartitioner, InitialPartitioningMode::recursive, 2>,
                         TestConfig<RecursiveInitialPartitioner, InitialPartitioningMode::recursive, 3>,
                         TestConfig<RecursiveInitialPartitioner, InitialPartitioningMode::recursive, 4>,
                         TestConfig<RecursiveInitialPartitioner, InitialPartitioningMode::recursive, 5>,
                         TestConfig<RecursiveBisectionInitialPartitioner, InitialPartitioningMode::recursive_bisection, 2>,
                         TestConfig<RecursiveBisectionInitialPartitioner, InitialPartitioningMode::recursive_bisection, 3>,
                         TestConfig<RecursiveBisectionInitialPartitioner, InitialPartitioningMode::recursive_bisection, 4>,
                         TestConfig<RecursiveBisectionInitialPartitioner, InitialPartitioningMode::recursive_bisection, 5> > TestConfigs;

TYPED_TEST_CASE(AInitialPartitionerTest, TestConfigs);

TYPED_TEST(AInitialPartitionerTest, VerifiesThatAllPartsAreNonEmpty) {
  this->initial_partitioner->initialPartition();

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_GT(this->partitioned_hypergraph.partWeight(part_id), 0);
  }
}

TYPED_TEST(AInitialPartitionerTest, VerifiesThatPartSizesAndWeightsAreCorrect) {
  this->initial_partitioner->initialPartition();

  std::vector<HypernodeWeight> part_weight(this->context.partition.k, 0);
  for ( const HypernodeID& hn : this->hypergraph.nodes() ) {
    PartitionID part_id = this->partitioned_hypergraph.partID(hn);
    part_weight[part_id] += this->partitioned_hypergraph.nodeWeight(hn);
  }

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
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
