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

#include <random>
#include "gmock/gmock.h"

#include "mt-kahypar/io/hypergraph_io.h"

#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"


using ::testing::Test;

namespace mt_kahypar {


template<typename Strategy>
vec<Gain> insertAndExtractAllMoves(Strategy& strat, PartitionedHypergraph& phg) {
  Move m;
  vec<Gain> gains;
  for (HypernodeID u : phg.nodes()) {
    strat.insertIntoPQ(phg, u, 0);
  }

  while (strat.findNextMove(phg, m)) {
    gains.push_back(m.gain);
  }
  strat.clearPQs(0);
  return gains;
}

TEST(StrategyTests, FindNextMove) {
  PartitionID k = 8;
  Context context;
  context.partition.k = k;
  context.partition.epsilon = 0.03;
  Hypergraph hg = io::readHypergraphFile("../tests/instances/contracted_ibm01.hgr", true);
  context.setupPartWeights(hg.totalWeight());
  PartitionedHypergraph phg = PartitionedHypergraph(k, hg);
  for (PartitionID i = 0; i < k; ++i) {
    context.partition.max_part_weights[i] = std::numeric_limits<HypernodeWeight>::max();
  }

  std::mt19937 rng(420);
  std::uniform_int_distribution<PartitionID> distr(0, k - 1);
  for (HypernodeID u : hg.nodes()) {
    phg.setOnlyNodePart(u, distr(rng));
  }
  phg.initializePartition();
  phg.initializeGainCache();


  context.refinement.fm.algorithm = FMAlgorithm::fm_gain_delta; // use this one because it allocates the most memory in shared data!

  FMSharedData sd(hg.initialNumNodes(), context);
  FMStats fm_stats;
  fm_stats.moves = 1;

  GainCacheStrategy gain_caching(context, hg.initialNumNodes(), sd, fm_stats);
  GainDeltaStrategy gain_deltas(context, hg.initialNumNodes(), sd, fm_stats);
  RecomputeGainStrategy recompute_gain(context, hg.initialNumNodes(), sd, fm_stats);
  GainCacheOnDemandStrategy gain_caching_on_demand(context, hg.initialNumNodes(), sd, fm_stats);

  vec<Gain> gains_from_deltas = insertAndExtractAllMoves(gain_deltas, phg);
  ASSERT_TRUE(std::is_sorted(gains_from_deltas.begin(), gains_from_deltas.end(), std::greater<Gain>()));

  vec<Gain> gains_cached = insertAndExtractAllMoves(gain_caching, phg);
  ASSERT_TRUE(std::is_sorted(gains_cached.begin(), gains_cached.end(), std::greater<Gain>()));
  ASSERT_EQ(gains_from_deltas, gains_cached);

  vec<Gain> gains_cached_on_demand = insertAndExtractAllMoves(gain_caching_on_demand, phg);
  ASSERT_TRUE(std::is_sorted(gains_cached_on_demand.begin(), gains_cached_on_demand.end(), std::greater<Gain>()));
  ASSERT_EQ(gains_cached_on_demand, gains_cached);

  vec<Gain> gains_recomputed = insertAndExtractAllMoves(recompute_gain, phg);
  ASSERT_TRUE(std::is_sorted(gains_recomputed.begin(), gains_recomputed.end(), std::greater<Gain>()));
  ASSERT_EQ(gains_recomputed, gains_cached);
}

TEST(StrategyTests, DeltaUpdatesWork) {
  PartitionID k = 8;
  Context context;
  context.partition.k = k;
  context.partition.epsilon = 0.03;
  Hypergraph hg = io::readHypergraphFile("../tests/instances/contracted_ibm01.hgr", true);
  context.setupPartWeights(hg.totalWeight());
  PartitionedHypergraph phg = PartitionedHypergraph(k, hg);
  for (PartitionID i = 0; i < k; ++i) {
    context.partition.max_part_weights[i] = std::numeric_limits<HypernodeWeight>::max();
  }

  std::mt19937 rng(420);
  std::uniform_int_distribution<PartitionID> distr(0, k - 1);
  for (HypernodeID u : hg.nodes()) {
    phg.setOnlyNodePart(u, distr(rng));
  }
  phg.initializePartition();
  phg.initializeGainCache();

  context.refinement.fm.algorithm = FMAlgorithm::fm_gain_delta; // use this one because it allocates the most memory in shared data!

  FMSharedData sd(hg.initialNumNodes(), context);
  FMStats fm_stats;

  GainDeltaStrategy strat(context, hg.initialNumNodes(), sd, fm_stats);
  for (HypernodeID u : hg.nodes())
    strat.insertIntoPQ(phg, u, 0);

  Move m;

  auto delta_func = [&](HyperedgeID e, HyperedgeWeight edge_weight, HypernodeID, HypernodeID pin_count_in_from_part_after,
          HypernodeID pin_count_in_to_part_after) {
    strat.deltaGainUpdates(phg, e, edge_weight, m.from, pin_count_in_from_part_after, m.to, pin_count_in_to_part_after);
  };

  auto check_gains = [&] {
    strat.doParallelForAllEntries([&](PartitionID to, HypernodeID u, Gain gain) {
      Gain re_gain = 0;
      PartitionID from = phg.partID(u);
      for (HyperedgeID e : phg.incidentEdges(u)) {
        if (phg.pinCountInPart(e, from) == 1) re_gain += phg.edgeWeight(e);
        if (phg.pinCountInPart(e, to) == 0) re_gain -= phg.edgeWeight(e);
      }
      ASSERT_EQ(gain, re_gain);
    });
  };

  check_gains();
  while (strat.findNextMove(phg, m)) {
    phg.changeNodePart(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func);
    fm_stats.moves++;
    check_gains();
  }
}


}