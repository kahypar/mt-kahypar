/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Robert Krause <robert.krause@student.kit.edu>
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
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_jet_refiner.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/rebalancing/simple_rebalancer.h"
#include "mt-kahypar/partition/refinement/rebalancing/advanced_rebalancer.h"
#include "mt-kahypar/partition/refinement/rebalancing/deterministic_rebalancer.h"
#include "mt-kahypar/utils/cast.h"

using ::testing::Test;

namespace mt_kahypar {
template <PartitionID k, Objective objective, RebalancingAlgorithm rebalancing>
struct TestConfig {};

template < PartitionID k, RebalancingAlgorithm rebalancing>
struct TestConfig< k, Objective::km1, rebalancing> {
    using TypeTraits = StaticHypergraphTypeTraits;
    using GainTypes = Km1GainTypes;
    using Refiner = DeterministicJetRefiner<GraphAndGainTypes<TypeTraits, GainTypes>>;
    static constexpr PartitionID K = k;
    static constexpr Objective OBJECTIVE = Objective::km1;
    static constexpr RebalancingAlgorithm REBALANCER = rebalancing;
};

template <PartitionID k, RebalancingAlgorithm rebalancing>
struct TestConfig<k, Objective::cut, rebalancing> {
    using TypeTraits = StaticHypergraphTypeTraits;
    using GainTypes = CutGainTypes;
    using Refiner = DeterministicJetRefiner<GraphAndGainTypes<TypeTraits, GainTypes>>;
    static constexpr PartitionID K = k;
    static constexpr Objective OBJECTIVE = Objective::cut;
    static constexpr RebalancingAlgorithm REBALANCER = rebalancing;
};

template <typename Config>
class ADeterministicJetRefiner : public Test {
    static size_t num_threads;

public:
    using TypeTraits = typename Config::TypeTraits;
    using GainTypes = typename Config::GainTypes;
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Refiner = typename Config::Refiner;
    using GainCache = typename GainTypes::GainCache;

    ADeterministicJetRefiner() :
        hypergraph(),
        partitioned_hypergraph(),
        context(),
        gain_cache(),
        refiner(nullptr),
        rebalancer(nullptr),
        metrics() {
        context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
        context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
        context.partition.mode = Mode::direct;
        context.partition.objective = Config::OBJECTIVE;
        context.partition.gain_policy = context.partition.objective == Objective::km1 ? GainPolicy::km1 : GainPolicy::cut;
        context.partition.epsilon = 0.25;
        context.partition.k = Config::K;
        context.partition.preset_type = PresetType::deterministic;
        context.partition.instance_type = InstanceType::hypergraph;
        context.partition.partition_type = PartitionedHypergraph::TYPE;
        context.partition.verbose_output = false;

        // Shared Memory
        context.shared_memory.num_threads = num_threads;

        // Initial Partitioning
        context.initial_partitioning.mode = Mode::deep_multilevel;
        context.initial_partitioning.runs = 1;

        // Jet
        context.refinement.rebalancing.algorithm = Config::REBALANCER;
        context.refinement.jet.num_iterations = 12;
        context.refinement.jet.relative_improvement_threshold = 0.001;
        context.refinement.jet.initial_negative_gain_factor = 0.75;
        context.refinement.jet.final_negative_gain_factor = 0.0;
        context.refinement.jet.dynamic_rounds = 3;


        // Read hypergraph
        hypergraph = io::readInputFile<Hypergraph>(
            "../tests/instances/contracted_unweighted_ibm01.hgr", FileFormat::hMetis, true);
        partitioned_hypergraph = PartitionedHypergraph(
            context.partition.k, hypergraph, parallel_tag_t());
        context.setupPartWeights(hypergraph.totalWeight());
        initialPartition();
        if (context.refinement.rebalancing.algorithm == RebalancingAlgorithm::simple_rebalancer) {
            rebalancer = std::make_unique<SimpleRebalancer<GraphAndGainTypes<TypeTraits, GainTypes>>>(
                context);
        } else if (context.refinement.rebalancing.algorithm == RebalancingAlgorithm::advanced_rebalancer) {
            rebalancer = std::make_unique<AdvancedRebalancer<GraphAndGainTypes<TypeTraits, GainTypes>>>(
                hypergraph.initialNumNodes(), context, gain_cache);
        } else if (context.refinement.rebalancing.algorithm == RebalancingAlgorithm::deterministic) {
            rebalancer = std::make_unique<DeterministicRebalancer<GraphAndGainTypes<TypeTraits, GainTypes>>>(
                hypergraph.initialNumNodes(), context, gain_cache);
        }
        refiner = std::make_unique<Refiner>(
            hypergraph.initialNumNodes(), hypergraph.initialNumEdges(), context, gain_cache, *rebalancer);
        mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hypergraph);
        rebalancer->initialize(phg);
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
    std::unique_ptr<IRebalancer> rebalancer;
    Metrics metrics;
};

template <typename Config>
size_t ADeterministicJetRefiner<Config>::num_threads = HardwareTopology::instance().num_cpus();

static constexpr double EPS = 0.05;

typedef ::testing::Types<
    TestConfig<2, Objective::cut, RebalancingAlgorithm::deterministic>,
    TestConfig<4, Objective::cut, RebalancingAlgorithm::deterministic>,
    TestConfig<8, Objective::cut, RebalancingAlgorithm::deterministic>,
    TestConfig<2, Objective::km1, RebalancingAlgorithm::deterministic>,
    TestConfig<2, Objective::km1, RebalancingAlgorithm::advanced_rebalancer>,
    TestConfig<4, Objective::km1, RebalancingAlgorithm::deterministic>,
    TestConfig<4, Objective::km1, RebalancingAlgorithm::advanced_rebalancer>,
    TestConfig<8, Objective::km1, RebalancingAlgorithm::deterministic>,
    TestConfig<8, Objective::km1, RebalancingAlgorithm::advanced_rebalancer>
> TestConfigs;

TYPED_TEST_CASE(ADeterministicJetRefiner, TestConfigs);

TYPED_TEST(ADeterministicJetRefiner, UpdatesImbalanceCorrectly) {
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
    this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
}

TYPED_TEST(ADeterministicJetRefiner, DoesNotViolateBalanceConstraint) {
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
    this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon + EPS);
}

TYPED_TEST(ADeterministicJetRefiner, UpdatesMetricsCorrectly) {
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
    this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_EQ(metrics::quality(this->partitioned_hypergraph, this->context.partition.objective),
        this->metrics.quality);
}

TYPED_TEST(ADeterministicJetRefiner, DoesNotWorsenSolutionQuality) {
    HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
    this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.quality, objective_before);
}


TYPED_TEST(ADeterministicJetRefiner, IncreasesTheNumberOfBlocks) {
    using PartitionedHypergraph = typename TestFixture::PartitionedHypergraph;
    HyperedgeWeight objective_before = metrics::quality(this->partitioned_hypergraph, this->context.partition.objective);
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
    this->refiner->refine(phg, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.quality, objective_before);

    // Initialize partition with smaller K
    const PartitionID old_k = this->context.partition.k;
    this->context.partition.k = std::max(old_k / 2, 2);
    this->context.setupPartWeights(this->hypergraph.totalWeight());
    PartitionedHypergraph phg_with_new_k(
        this->context.partition.k, this->hypergraph, mt_kahypar::parallel_tag_t());
    vec<PartitionID> non_optimized_partition(this->hypergraph.initialNumNodes(), kInvalidPartition);
    this->partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        // create a semi-random partition
        const PartitionID block = this->partitioned_hypergraph.partID(hn);
        phg_with_new_k.setOnlyNodePart(hn, (block + hn) % this->context.partition.k);
        non_optimized_partition[hn] = phg_with_new_k.partID(hn);
    });
    phg_with_new_k.initializePartition();
    this->metrics.quality = metrics::quality(phg_with_new_k, this->context);
    this->metrics.imbalance = metrics::imbalance(phg_with_new_k, this->context);

    objective_before = metrics::quality(phg_with_new_k, this->context.partition.objective);
    mt_kahypar_partitioned_hypergraph_t phg_new_k = utils::partitioned_hg_cast(phg_with_new_k);
    this->gain_cache.reset(phg_with_new_k.initialNumNodes(), phg_with_new_k.k());
    this->refiner->initialize(phg_new_k);
    this->rebalancer->initialize(phg_new_k);
    this->refiner->refine(phg_new_k, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.quality, objective_before);
    ASSERT_EQ(metrics::quality(phg_with_new_k, this->context.partition.objective),
        this->metrics.quality);

    // Check if refiner has moved some nodes
    bool has_moved_nodes = false;
    for (const HypernodeID hn : phg_with_new_k.nodes()) {
        if (non_optimized_partition[hn] != phg_with_new_k.partID(hn)) {
            has_moved_nodes = true;
            break;
        }
    }
    ASSERT_TRUE(has_moved_nodes);
}

}  // namespace mt_kahypar