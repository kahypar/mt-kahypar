/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
  Hypergraph hg = io::readHypergraphFile("../tests/instances/contracted_ibm01.hgr", 0, true);
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
  phg.initializePartition(0);
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
  Hypergraph hg = io::readHypergraphFile("../tests/instances/contracted_ibm01.hgr", 0, true);
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
  phg.initializePartition(0);
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