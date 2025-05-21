/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <functional>
#include <random>


#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/refinement/rebalancing/deterministic_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {

template <typename TypeTraitsT, typename GainTypesT, PartitionID k>
struct TestConfig {
  using TypeTraits = TypeTraitsT;
  using GainTypes = GainTypesT;
  static constexpr PartitionID K = k;
};

template<typename Config>
class DeterministicRebalancerTest : public Test {

 public:
  using TypeTraits = typename Config::TypeTraits;
  using GainTypes = typename Config::GainTypes;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;
  using GainCache = typename GainTypes::GainCache;
  using Rebalancer = DeterministicRebalancer<GraphAndGainTypes<TypeTraits, GainTypes>>;

  DeterministicRebalancerTest() :
          hypergraph(),
          partitioned_hypergraph(),
          context(),
          rebalancer(nullptr) {
    TBBInitializer::instance(std::thread::hardware_concurrency());
    context.partition.mode = Mode::direct;
    context.partition.epsilon = 0.05;
    context.partition.k = Config::K;

    context.partition.preset_type = PresetType::default_preset;
    context.partition.instance_type = InstanceType::hypergraph;
    context.partition.partition_type = PartitionedHypergraph::TYPE;
    context.partition.verbose_output = false;

    context.refinement.rebalancing.det_heavy_vertex_exclusion_factor = 1.5;
    context.refinement.rebalancing.det_relative_deadzone_size = 0.25;

    // Shared Memory
    context.shared_memory.original_num_threads = std::thread::hardware_concurrency();
    context.shared_memory.num_threads = std::thread::hardware_concurrency();

    context.partition.objective = Hypergraph::is_graph ? Objective::cut : Objective::km1;
    context.partition.gain_policy = Hypergraph::is_graph ? GainPolicy::cut_for_graphs : GainPolicy::km1;
  }

  void constructFromFile() {
    if constexpr ( Hypergraph::is_graph ) {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/delaunay_n10.graph", FileFormat::Metis, true);
    } else {
      hypergraph = io::readInputFile<Hypergraph>(
        "../tests/instances/contracted_unweighted_ibm01.hgr", FileFormat::hMetis, true);
    }
  }

  void constructFromValues(const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,
                           const vec<vec<HypernodeID>>& edge_vector, const vec<HypernodeWeight> hypernode_weight) {
    hypergraph = HypergraphFactory::construct(num_hypernodes, num_hyperedges, edge_vector, nullptr, hypernode_weight.data());
  }

  void setup() {
    partitioned_hypergraph = PartitionedHypergraph(context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());

    rebalancer = std::make_unique<Rebalancer>(hypergraph.initialNumNodes(), context);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<Rebalancer> rebalancer;
};


typedef ::testing::Types<TestConfig<StaticHypergraphTypeTraits, Km1GainTypes, 2>,
                         TestConfig<StaticHypergraphTypeTraits, Km1GainTypes, 4>
                         ENABLE_GRAPHS(COMMA TestConfig<StaticGraphTypeTraits COMMA CutGainForGraphsTypes COMMA 2>)
                         ENABLE_GRAPHS(COMMA TestConfig<StaticGraphTypeTraits COMMA CutGainForGraphsTypes COMMA 4>) > TestConfigs;

TYPED_TEST_CASE(DeterministicRebalancerTest, TestConfigs);


TYPED_TEST(DeterministicRebalancerTest, CanNotBeRebalanced) {
  this->constructFromValues(3, 1, { {0, 1} }, {6, 5, 4});
  this->setup();

  this->partitioned_hypergraph.setOnlyNodePart(0, 0);
  this->partitioned_hypergraph.setOnlyNodePart(1, 1);
  this->partitioned_hypergraph.setOnlyNodePart(2, 0);
  this->partitioned_hypergraph.initializePartition();
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->rebalancer->initialize(phg);

  Metrics metrics;
  metrics.quality = metrics::quality(this->partitioned_hypergraph, this->context);
  metrics.imbalance = metrics::imbalance(this->partitioned_hypergraph, this->context);
  this->rebalancer->refine(phg, {}, metrics, std::numeric_limits<double>::max());

  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), metrics.imbalance);
}


TYPED_TEST(DeterministicRebalancerTest, MustMoveHeavyVertex) {
  // tests a specific edge case which can only be solved if the rebalancer moves vertices which are usually excluded
  this->constructFromValues(7, 1, { {0, 1} }, {4, 4, 4, 3, 4, 3, 3});
  this->setup();

  this->partitioned_hypergraph.setOnlyNodePart(0, 0);
  this->partitioned_hypergraph.setOnlyNodePart(1, 1);
  this->partitioned_hypergraph.setOnlyNodePart(2, 0);
  this->partitioned_hypergraph.setOnlyNodePart(3, 1);
  this->partitioned_hypergraph.setOnlyNodePart(4, 2 % this->context.partition.k);
  this->partitioned_hypergraph.setOnlyNodePart(5, 2 % this->context.partition.k);
  this->partitioned_hypergraph.setOnlyNodePart(6, 3 % this->context.partition.k);
  this->partitioned_hypergraph.initializePartition();
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->rebalancer->initialize(phg);

  Metrics metrics;
  metrics.quality = metrics::quality(this->partitioned_hypergraph, this->context);
  metrics.imbalance = metrics::imbalance(this->partitioned_hypergraph, this->context);
  this->rebalancer->refine(phg, {}, metrics, std::numeric_limits<double>::max());

  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), metrics.imbalance);
  for (PartitionID part = 0; part < this->context.partition.k; ++part) {
    ASSERT_LE(this->partitioned_hypergraph.partWeight(part), this->context.partition.max_part_weights[part]);
  }
}


TYPED_TEST(DeterministicRebalancerTest, ProducesBalancedResult) {
  this->constructFromFile();
  this->setup();

  this->partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
    PartitionID block = 0;
    for (PartitionID p = 1; p < this->context.partition.k; ++p) {
      if (utils::Randomize::instance().flipCoin(THREAD_ID)) {
        block++;
      }
    }
    this->partitioned_hypergraph.setOnlyNodePart(hn, block);
  });

  this->partitioned_hypergraph.initializePartition();
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(this->partitioned_hypergraph);
  this->rebalancer->initialize(phg);

  Metrics metrics;
  metrics.quality = metrics::quality(this->partitioned_hypergraph, this->context);
  metrics.imbalance = metrics::imbalance(this->partitioned_hypergraph, this->context);
  this->rebalancer->refine(phg, {}, metrics, std::numeric_limits<double>::max());

  ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), metrics.imbalance);
  for (PartitionID part = 0; part < this->context.partition.k; ++part) {
    ASSERT_LE(this->partitioned_hypergraph.partWeight(part), this->context.partition.max_part_weights[part]);
  }
}

}
