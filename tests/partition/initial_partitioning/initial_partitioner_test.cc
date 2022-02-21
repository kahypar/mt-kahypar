/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/initial_partitioning/deep_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_bipartitioning_initial_partitioner.h"

using ::testing::Test;

namespace mt_kahypar {


template <class InitialPartitioner, Mode mode, PartitionID k>
struct TestConfig {
  using Partitioner = InitialPartitioner;
  static constexpr Mode MODE = mode;
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

    if ( context.partition.paradigm == Paradigm::multilevel ) {
      parseIniToContext(context, "../config/default_preset.ini");
    } else {
      parseIniToContext(context, "../config/quality_preset.ini");
    }

    context.partition.graph_filename = "../tests/instances/contracted_unweighted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = Mode::direct;
    context.partition.objective = kahypar::Objective::km1;
    context.partition.epsilon = 0.2;
    context.partition.k = Config::K;
    context.partition.verbose_output = true;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Coarsening
    if ( context.partition.paradigm == Paradigm::nlevel ) {
      context.refinement.max_batch_size = 100;
    }

    // Initial Partitioning
    context.initial_partitioning.runs = 1;
    context.sparsification.use_degree_zero_contractions = false;
    context.sparsification.use_heavy_net_removal = false;
    context.sparsification.use_similiar_net_removal = false;
    context.sparsification.use_degree_zero_contractions = false;
    context.initial_partitioning.mode = Config::MODE;
    context.initial_partitioning.remove_degree_zero_hns_before_ip = false;

    // Label Propagation
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    context.initial_partitioning.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;

    // FM
    context.refinement.fm.algorithm = FMAlgorithm::do_nothing;
    context.initial_partitioning.refinement.fm.algorithm = FMAlgorithm::do_nothing;

    // Flows
    context.refinement.flows.algorithm = FlowAlgorithm::do_nothing;
    context.initial_partitioning.refinement.flows.algorithm = FlowAlgorithm::do_nothing;

    // Read hypergraph
    hypergraph = io::readHypergraphFile(
      "../tests/instances/contracted_unweighted_ibm01.hgr");
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    context.setupContractionLimit(hypergraph.totalWeight());
    assignCommunities();

    initial_partitioner = std::make_unique<InitialPartitioner>(partitioned_hypergraph, context);
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

typedef ::testing::Types<TestConfig<DeepInitialPartitioner, Mode::deep_multilevel, 2>,
                         TestConfig<DeepInitialPartitioner, Mode::deep_multilevel, 3>,
                         TestConfig<DeepInitialPartitioner, Mode::deep_multilevel, 4>,
                         TestConfig<RecursiveBipartitioningInitialPartitioner, Mode::recursive_bipartitioning, 2>,
                         TestConfig<RecursiveBipartitioningInitialPartitioner, Mode::recursive_bipartitioning, 3>,
                         TestConfig<RecursiveBipartitioningInitialPartitioner, Mode::recursive_bipartitioning, 4> > TestConfigs;

TYPED_TEST_CASE(AInitialPartitionerTest, TestConfigs);

TYPED_TEST(AInitialPartitionerTest, VerifiesComputedPartition) {
  this->initial_partitioner->initialPartition();

  // Check that each vertex is assigned to a block
  for ( const HypernodeID& hn : this->partitioned_hypergraph.nodes() ) {
    ASSERT_NE(this->partitioned_hypergraph.partID(hn), kInvalidPartition)
      << "Hypernode " << hn << " is unassigned!";
  }

  // Check that non of the blocks is empty
  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_GT(this->partitioned_hypergraph.partWeight(part_id), 0)
      << "Block " << part_id << " is empty!";
  }

  // Check that part weights are correct
  std::vector<HypernodeWeight> part_weight(this->context.partition.k, 0);
  for ( const HypernodeID& hn : this->hypergraph.nodes() ) {
    PartitionID part_id = this->partitioned_hypergraph.partID(hn);
    ASSERT(part_id >= 0 && part_id < this->context.partition.k);
    part_weight[part_id] += this->partitioned_hypergraph.nodeWeight(hn);
  }

  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_EQ(this->partitioned_hypergraph.partWeight(part_id), part_weight[part_id])
      << "Expected part weight of block " << part_id << " is " << part_weight[part_id]
      << ", but currently is " << this->partitioned_hypergraph.partWeight(part_id);
  }

  // Check that balance constraint is fullfilled
  for ( PartitionID part_id = 0; part_id < this->context.partition.k; ++part_id ) {
    ASSERT_LE(this->partitioned_hypergraph.partWeight(part_id),
              this->context.partition.max_part_weights[part_id])
      << "Block " << part_id << " violates the balance constraint (Part Weight = "
      << this->partitioned_hypergraph.partWeight(part_id) << ", Max Part Weight = "
      << this->context.partition.max_part_weights[part_id];
  }
}

}  // namespace mt_kahypar
