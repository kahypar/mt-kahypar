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

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/deep_multilevel.h"

using ::testing::Test;

namespace mt_kahypar {


template <Mode mode, PartitionID k>
struct TestConfig {
  static constexpr Mode MODE = mode;
  static constexpr PartitionID K = k;
};

template <typename Config>
class AInitialPartitionerTest : public Test {

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
    context.partition.objective = Objective::km1;
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
  }

  void assignCommunities() {
    std::vector<PartitionID> communities;
    io::readPartitionFile(context.partition.graph_community_filename, communities);

    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, communities[hn]);
    }
  }

  void runInitialPartitioning() {
    switch ( context.initial_partitioning.mode ) {
      case Mode::recursive_bipartitioning:
        recursive_bipartitioning::partition(partitioned_hypergraph, context); break;
      case Mode::deep_multilevel:
        deep_multilevel::partition(partitioned_hypergraph, context); break;
      case Mode::direct:
      case Mode::UNDEFINED:
        ERROR("Undefined initial partitioning algorithm.");
    }
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
};

template <typename Config>
size_t AInitialPartitionerTest<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<Mode::deep_multilevel, 2>,
                         TestConfig<Mode::deep_multilevel, 3>,
                         TestConfig<Mode::deep_multilevel, 4>,
                         TestConfig<Mode::recursive_bipartitioning, 2>,
                         TestConfig<Mode::recursive_bipartitioning, 3>,
                         TestConfig<Mode::recursive_bipartitioning, 4> > TestConfigs;

TYPED_TEST_CASE(AInitialPartitionerTest, TestConfigs);

TYPED_TEST(AInitialPartitionerTest, VerifiesComputedPartition) {
  this->runInitialPartitioning();

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
