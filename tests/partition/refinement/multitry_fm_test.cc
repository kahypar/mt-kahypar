/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"

using ::testing::Test;

namespace mt_kahypar {

template <typename TypeTraitsT, PartitionID k>
struct TestConfig {
  using TypeTraits = TypeTraitsT;
  static constexpr PartitionID K = k;
};

template<typename Config>
class MultiTryFMTest : public Test {

 public:
  using TypeTraits = typename Config::TypeTraits;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using Refiner = MultiTryKWayFM<TypeTraits, Km1GainTypes>;

  MultiTryFMTest() :
          hypergraph(),
          partitioned_hypergraph(),
          context(),
          gain_cache(),
          refiner(nullptr),
          metrics() {
    TBBInitializer::instance(std::thread::hardware_concurrency());
    context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = Mode::direct;
    context.partition.epsilon = 0.25;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.original_num_threads = std::thread::hardware_concurrency();
    context.shared_memory.num_threads = std::thread::hardware_concurrency();

    // Initial Partitioning
    context.initial_partitioning.mode = Mode::deep_multilevel;
    context.initial_partitioning.runs = 1;

    context.partition.k = Config::K;

    context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
    context.refinement.fm.multitry_rounds = 10;
    context.refinement.fm.num_seed_nodes = 5;
    context.refinement.fm.rollback_balance_violation_factor = 1.0;

    context.partition.objective = Objective::km1;
    context.partition.gain_policy = GainPolicy::km1;

    // Read hypergraph
    hypergraph = io::readInputFile<Hypergraph>(
      "../tests/instances/contracted_unweighted_ibm01.hgr", FileFormat::hMetis, true);
    partitioned_hypergraph = PartitionedHypergraph(
            context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    initialPartition();

    refiner = std::make_unique<Refiner>(hypergraph.initialNumNodes(),
      hypergraph.initialNumEdges(), context, gain_cache);
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hypergraph);
    refiner->initialize(phg);
  }

  void initialPartition() {
    Context ip_context(context);
    ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    InitialPartitioningDataContainer<TypeTraits> ip_data(partitioned_hypergraph, ip_context);
    ip_data_container_t* ip_data_ptr = ip::to_pointer(ip_data);
    BFSInitialPartitioner<TypeTraits> initial_partitioner(
      InitialPartitioningAlgorithm::bfs, ip_data_ptr, ip_context, 420, 0);
    initial_partitioner.partition();
    ip_data.apply();
    metrics.quality = metrics::quality(partitioned_hypergraph, context);
    metrics.imbalance = metrics::imbalance(partitioned_hypergraph, context);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  Km1GainCache gain_cache;
  std::unique_ptr<Refiner> refiner;
  Metrics metrics;
};


typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits, 2>,
                         TestConfig<StaticHypergraphTypeTraits, 4>,
                         TestConfig<StaticHypergraphTypeTraits, 8>,
                         TestConfig<StaticHypergraphTypeTraits, 2>,
                         TestConfig<StaticHypergraphTypeTraits, 4>,
                         TestConfig<StaticHypergraphTypeTraits, 8>
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 2>)
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 4>)
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 8>)
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 2>)
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 4>)
                         ENABLE_HIGHEST_QUALITY(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 8>) > TestConfigs;

TYPED_TEST_CASE(MultiTryFMTest, TestConfigs);

TYPED_TEST(MultiTryFMTest, UpdatesImbalanceCorrectly) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
}


TYPED_TEST(MultiTryFMTest, DoesNotViolateBalanceConstraint) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
}

TYPED_TEST(MultiTryFMTest, UpdatesMetricsCorrectly) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.quality);
}

TYPED_TEST(MultiTryFMTest, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);
}

TYPED_TEST(MultiTryFMTest, AlsoWorksWithNonDefaultFeatures) {
  this->context.refinement.fm.obey_minimal_parallelism = true;
  this->context.refinement.fm.rollback_parallel = false;
  this->context.refinement.fm.perform_moves_global = true;
  HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.quality);
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
}

TYPED_TEST(MultiTryFMTest, WorksWithRefinementNodes) {
  parallel::scalable_vector<HypernodeID> refinement_nodes;
  for (HypernodeID u = 0; u < this->partitioned_hypergraph.initialNumNodes(); ++u) {
    refinement_nodes.push_back(u);
  }
  HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, refinement_nodes, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.quality);
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);

  std::stringstream buffer;
  std::streambuf* old = std::cout.rdbuf(buffer.rdbuf());  // redirect std::cout to discard output
  this->refiner->printMemoryConsumption();
  std::cout.rdbuf(old);                                   // and reset again
}

/*TYPED_TEST(MultiTryFMTest, IncreasesTheNumberOfBlocks) {
  using PartitionedHypergraph = typename TestFixture::PartitionedHypergraph;
  HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);

  // Initialize partition with larger K
  const PartitionID old_k = this->context.partition.k;
  this->context.partition.k = 2 * old_k;
  this->context.setupPartWeights(this->hypergraph.totalWeight());
  PartitionedHypergraph phg_with_larger_k(
    this->context.partition.k, this->hypergraph, mt_kahypar::parallel_tag_t());
  utils::Randomize& rand = utils::Randomize::instance();
  vec<PartitionID> non_optimized_partition(this->hypergraph.initialNumNodes(), kInvalidPartition);
  this->partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
    const PartitionID block = this->partitioned_hypergraph.partID(hn);
    phg_with_larger_k.setOnlyNodePart(hn, rand.flipCoin(SCHED_GETCPU) ? 2 * block : 2 * block + 1);
    non_optimized_partition[hn] = phg_with_larger_k.partID(hn);
  });
  phg_with_larger_k.initializePartition();
  this->metrics.quality = metrics::quality(phg_with_larger_k, this->context);
  this->metrics.imbalance = metrics::imbalance(phg_with_larger_k, this->context);

  objective_before = metrics::quality(phg_with_larger_k, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg_larger_k = utils::partitioned_hg_cast(phg_with_larger_k);
  this->refiner->initialize(phg_larger_k);
  this->refiner->refine(phg_larger_k, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);
  ASSERT_EQ(metrics::quality(phg_with_larger_k, this->context.partition.objective),
            this->metrics.quality);
}*/

}  // namespace mt_kahypar
