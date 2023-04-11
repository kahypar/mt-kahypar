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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"

#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_cache.h"

using ::testing::Test;

namespace mt_kahypar {

namespace {
  using Hypergraph = typename StaticHypergraphTypeTraits::Hypergraph;
  using PartitionedHypergraph = typename StaticHypergraphTypeTraits::PartitionedHypergraph;
}


template<typename Strategy>
vec<Gain> insertAndExtractAllMoves(Strategy& strat,
                                   PartitionedHypergraph& phg,
                                   Km1GainCache& gain_cache) {
  Move m;
  vec<Gain> gains;
  for (HypernodeID u : phg.nodes()) {
    strat.insertIntoPQ(phg, gain_cache, u);
  }

  while (strat.findNextMove(phg, gain_cache, m)) {
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
  Hypergraph hg = io::readInputFile<Hypergraph>(
    "../tests/instances/contracted_ibm01.hgr", FileFormat::hMetis, true);
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
  Km1GainCache gain_cache;
  gain_cache.initializeGainCache(phg);


  context.refinement.fm.algorithm = FMAlgorithm::kway_fm;

  FMSharedData sd(hg.initialNumNodes());
  FMStats fm_stats;
  fm_stats.moves = 1;

  GainCacheStrategy gain_caching(context, sd, fm_stats);
  vec<Gain> gains_cached = insertAndExtractAllMoves(gain_caching, phg, gain_cache);
  ASSERT_TRUE(std::is_sorted(gains_cached.begin(), gains_cached.end(), std::greater<Gain>()));
}

}