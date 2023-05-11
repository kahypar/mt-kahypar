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
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/registries/register_refinement_algorithms.cpp"
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/jet/jet_refiner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/cast.h"

using ::testing::Test;

namespace mt_kahypar {
template <typename TypeTraitsT, PartitionID k, Objective objective, bool precomputed>
struct TestConfig { };

template <typename TypeTraitsT, PartitionID k, bool precomputed>
struct TestConfig<TypeTraitsT, k, Objective::km1, precomputed> {
  using TypeTraits = TypeTraitsT;
  using GainTypes = Km1GainTypes;
  using Refiner = JetRefiner<TypeTraits, Km1GainTypes, precomputed>;
  static constexpr PartitionID K = k;
  static constexpr Objective OBJECTIVE = Objective::km1;
  static constexpr JetAlgorithm JET_ALGO = JetAlgorithm::greedy_unordered;
};

template <typename TypeTraitsT, PartitionID k, bool precomputed>
struct TestConfig<TypeTraitsT, k, Objective::cut, precomputed> {
  using TypeTraits = TypeTraitsT;
  using GainTypes = CutGainTypes;
  using Refiner = JetRefiner<TypeTraits, CutGainTypes, precomputed>;
  static constexpr PartitionID K = k;
  static constexpr Objective OBJECTIVE = Objective::cut;
  static constexpr JetAlgorithm JET_ALGO = JetAlgorithm::greedy_unordered;
};

template <typename Config>
class AJetRefiner : public Test {
  static size_t num_threads;

 public:
  using TypeTraits = typename Config::TypeTraits;
  using GainTypes = typename Config::GainTypes;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using Refiner = typename Config::Refiner;
  using GainCache = typename GainTypes::GainCache;

  AJetRefiner() :
    hypergraph(),
    partitioned_hypergraph(),
    context(),
    gain_cache(),
    refiner(nullptr),
    metrics() {
    context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
    context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
    context.partition.mode = Mode::direct;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.gain_policy = context.partition.objective == Objective::km1 ? GainPolicy::km1 : GainPolicy::cut;
    context.partition.epsilon = 0.25;
    context.partition.k = Config::K;
    #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
    context.partition.preset_type = Hypergraph::is_static_hypergraph ?
     PresetType::default_preset : PresetType::quality_preset;
    #else
    context.partition.preset_type = PresetType::default_preset;
    #endif
    context.partition.instance_type = InstanceType::hypergraph;
    context.partition.partition_type = PartitionedHypergraph::TYPE;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;

    // Initial Partitioning
    context.initial_partitioning.mode = Mode::deep_multilevel;
    context.initial_partitioning.runs = 1;

    // Jet
    context.refinement.jet.algorithm = Config::JET_ALGO;
    context.initial_partitioning.refinement.jet.algorithm = Config::JET_ALGO;
    context.initial_partitioning.refinement.jet.vertex_locking = false;

    // Read hypergraph
    hypergraph = io::readInputFile<Hypergraph>(
      "../tests/instances/contracted_unweighted_ibm01.hgr", FileFormat::hMetis, true);
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    initialPartition();

    refiner = std::make_unique<Refiner>(hypergraph.initialNumNodes(), context, gain_cache);
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
  GainCache gain_cache;
  std::unique_ptr<Refiner> refiner;
  Metrics metrics;
};

template <typename Config>
size_t AJetRefiner<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits, 2, Objective::cut, false>,
                         TestConfig<StaticHypergraphTypeTraits, 4, Objective::cut, false>,
                         TestConfig<StaticHypergraphTypeTraits, 8, Objective::cut, false>,
                         TestConfig<StaticHypergraphTypeTraits, 2, Objective::km1, false>,
                         TestConfig<StaticHypergraphTypeTraits, 4, Objective::km1, false>,
                         TestConfig<StaticHypergraphTypeTraits, 8, Objective::km1, false>,
                         TestConfig<StaticHypergraphTypeTraits, 2, Objective::cut, true>,
                         TestConfig<StaticHypergraphTypeTraits, 4, Objective::cut, true>,
                         TestConfig<StaticHypergraphTypeTraits, 8, Objective::cut, true>,
                         TestConfig<StaticHypergraphTypeTraits, 2, Objective::km1, true>,
                         TestConfig<StaticHypergraphTypeTraits, 4, Objective::km1, true>,
                         TestConfig<StaticHypergraphTypeTraits, 8, Objective::km1, true>
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 2 COMMA Objective::cut COMMA false>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 4 COMMA Objective::cut COMMA false>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 8 COMMA Objective::cut COMMA false>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 2 COMMA Objective::km1 COMMA false>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 4 COMMA Objective::km1 COMMA false>)
                         ENABLE_N_LEVEL(COMMA TestConfig<DynamicHypergraphTypeTraits COMMA 8 COMMA Objective::km1 COMMA false>) > TestConfigs;

TYPED_TEST_CASE(AJetRefiner, TestConfigs);

TYPED_TEST(AJetRefiner, UpdatesImbalanceCorrectly) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
}

TYPED_TEST(AJetRefiner, DoesNotViolateBalanceConstraint) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon + EPS);
}

TYPED_TEST(AJetRefiner, UpdatesMetricsCorrectly) {
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context.partition.objective),
            this->metrics.quality);
}

TYPED_TEST(AJetRefiner, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
  ASSERT_LE(this->metrics.quality, objective_before);
}


TYPED_TEST(AJetRefiner, IncreasesTheNumberOfBlocks) {
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

  // Check if refiner has moved some nodes from new blocks
  bool has_moved_nodes = false;
  for ( const HypernodeID hn : phg_with_larger_k.nodes() ) {
    if ( non_optimized_partition[hn] >= old_k &&
         non_optimized_partition[hn] != phg_with_larger_k.partID(hn) ) {
      has_moved_nodes = true;
      break;
    }
  }
  ASSERT_TRUE(has_moved_nodes);
}

}  // namespace mt_kahypar
