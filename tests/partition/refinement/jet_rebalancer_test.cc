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

#include <memory>

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/partition/refinement/rebalancing/jet_rebalancer.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {

namespace {
  using TypeTraits = StaticHypergraphTypeTraits;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename Km1GainTypes::GainCache;
}

template<template <class, class> typename RebalancerT, bool use_weak_v, bool use_greedy_v, bool use_deadzone_v>
struct Config {
  using Rebalancer = RebalancerT<TypeTraits, Km1GainTypes>;
  static constexpr bool use_weak = use_weak_v;
  static constexpr bool use_greedy = use_greedy_v;
  static constexpr bool use_deadzone = use_deadzone_v;
};

template <typename Config>
struct ARebalancerTest : public Test {
  using Rebalancer = typename Config::Rebalancer;

  void setupContext(Context& context) {
    context.refinement.jet_rebalancing.num_weak_iterations = Config::use_weak ? 2 : 0;
    context.refinement.jet_rebalancing.use_greedy_balanced_instead_of_strong_iteration = Config::use_greedy;
    context.refinement.jet_rebalancing.greedy_balanced_use_deadzone = Config::use_deadzone;
  }

  Rebalancer& getRebalancer(HypernodeID num_nodes, const Context& context, GainCache& gain_cache) {
    rebalancer = std::make_unique<Rebalancer>(num_nodes, context, gain_cache);
    return *rebalancer;
  }

  std::unique_ptr<Rebalancer> rebalancer;
};

typedef ::testing::Types< Config<JetRebalancer, true, false, true>,
                          Config<JetRebalancer, true, true, true>,
                          Config<JetRebalancer, true, true, false>,
                          Config<JetRebalancer, false, false, true>,
                          Config<JetRebalancer, false, true, true>,
                          Config<JetRebalancer, false, true, false> > TestConfigs;

TYPED_TEST_CASE(ARebalancerTest, TestConfigs);

TYPED_TEST(ARebalancerTest, BalancesBlockWithCompleteGraph) {
  Hypergraph hg = io::readInputFile<Hypergraph>(
    "../tests/instances/contracted_ibm01.hgr", FileFormat::hMetis,
    true /* enable stable construction */);

  for (PartitionID k: {2, 4, 8, 16}) {
    Context context;
    context.partition.k = k;
    context.partition.epsilon = (k <= 8) ? 0.01 : 0.02;
    context.partition.objective = Objective::km1;
    context.setupPartWeights(hg.totalWeight());
    this->setupContext(context);
    PartitionedHypergraph phg = PartitionedHypergraph(k, hg);

    for (HypernodeID u = 0; u < hg.initialNumNodes(); ++u) {
      phg.setOnlyNodePart(u, 0);
    }
    phg.initializePartition();
    Km1GainCache gain_cache;
    gain_cache.initializeGainCache(phg);

    Metrics current_metrics {metrics::quality(phg, context), metrics::imbalance(phg, context)};
    auto& rebalancer = this->getRebalancer(hg.initialNumNodes(), context, gain_cache);

    mt_kahypar_partitioned_hypergraph_t tmp_phg = utils::partitioned_hg_cast(phg);
    rebalancer.initialize(tmp_phg);
    rebalancer.refine(tmp_phg, {}, current_metrics, 0.0);
    ASSERT_LE(metrics::imbalance(phg, context), context.partition.epsilon);
  }
}

TYPED_TEST(ARebalancerTest, BalancesRandomAssignment) {
  Hypergraph hg = io::readInputFile<Hypergraph>(
    "../tests/instances/contracted_ibm01.hgr", FileFormat::hMetis,
    true /* enable stable construction */);
  utils::Randomize& rand = utils::Randomize::instance();
  rand.setSeed(0);

  for (PartitionID k: {2, 4, 8, 16}) {
    Context context;
    context.partition.k = k;
    context.partition.epsilon = (k <= 8) ? 0.01 : 0.02;
    context.partition.objective = Objective::km1;
    context.setupPartWeights(hg.totalWeight());
    PartitionedHypergraph phg = PartitionedHypergraph(k, hg);

    for (HypernodeID u = 0; u < hg.initialNumNodes(); ++u) {
      phg.setOnlyNodePart(u, rand.getRandomInt(0, k - 1, SCHED_GETCPU));
    }
    phg.initializePartition();
    Km1GainCache gain_cache;
    gain_cache.initializeGainCache(phg);

    Metrics current_metrics {metrics::quality(phg, context), metrics::imbalance(phg, context)};
    auto& rebalancer = this->getRebalancer(hg.initialNumNodes(), context, gain_cache);

    mt_kahypar_partitioned_hypergraph_t tmp_phg = utils::partitioned_hg_cast(phg);
    rebalancer.initialize(tmp_phg);
    rebalancer.refine(tmp_phg, {}, current_metrics, 0.0);
    ASSERT_LE(metrics::imbalance(phg, context), context.partition.epsilon);
  }
}

}
