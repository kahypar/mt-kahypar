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

#include <functional>
#include <random>

#include "gmock/gmock.h"

#include "mt-kahypar/macros.h"


#include "mt-kahypar/io/hypergraph_io.h"

#include "mt-kahypar/partition/refinement/fm/global_rollback.h"

#include "mt-kahypar/partition/metrics.h"

using ::testing::Test;

namespace mt_kahypar {

#ifndef KAHYPAR_USE_N_LEVEL_PARADIGM

TEST(RollbackTests, GainRecalculationAndRollsbackCorrectly) {
  Hypergraph hg = io::readHypergraphFile("../tests/instances/twocenters.hgr", 0);
  PartitionID k = 2;


  PartitionedHypergraph phg(k, hg);
  phg.setNodePart(0, 1);
  phg.setNodePart(1, 1);
  for (HypernodeID u = 4; u < 12; ++u) {
    phg.setNodePart(u, 0);
  }
  phg.setNodePart(2, 0);
  phg.setNodePart(3, 0);
  for (HypernodeID u = 12; u < 20; ++u) {
    phg.setNodePart(u, 1);
  }
    phg.initializeGainCache();

  Context context;
  context.partition.k = k;
  context.setupPartWeights(phg.totalWeight());
  context.partition.max_part_weights = { std::numeric_limits<HypernodeWeight>::max(), std::numeric_limits<HypernodeWeight>::max()};
  context.refinement.fm.rollback_balance_violation_factor = 0.0;


  FMSharedData sharedData(hg.initialNumNodes(), context);

  GlobalRollback grb(hg, context, k);
  grb.setRemainingOriginalPins(phg);
  auto performMove = [&](Move m) {
    if (phg.changeNodePartWithGainCacheUpdate(m.node, m.from, m.to)) {
      sharedData.moveTracker.insertMove(m);
    }
  };

  ASSERT_EQ(3, phg.km1Gain(0, 1, 0));
  performMove({1, 0,  0,  3});
  ASSERT_EQ(3, phg.km1Gain(1, 1, 0));
  performMove({1, 0,  1,  3});
  ASSERT_EQ(1, phg.km1Gain(2, 0, 1));
  performMove({0, 1,  2,  1});
  ASSERT_EQ(1, phg.km1Gain(3, 0, 1));
  performMove({0, 1,  3,  1});

  ASSERT_EQ(phg.km1Gain(4, 0, 1), -1);
  performMove({0, 1, 4, -1});
  ASSERT_EQ(phg.km1Gain(5, 0, 1), 0);
  performMove({0, 1, 5, 0});

  vec<HypernodeWeight> dummy_part_weights(k, 0);
  grb.revertToBestPrefix<true>(phg, sharedData, dummy_part_weights);
  // revert last two moves
  ASSERT_EQ(phg.partID(4), 0);
  ASSERT_EQ(phg.partID(5), 0);
  ASSERT_EQ(metrics::km1(phg, false), 2);
}


TEST(RollbackTests, GainRecalculation2) {
  Hypergraph hg = io::readHypergraphFile("../tests/instances/twocenters.hgr", 0);
  PartitionID k = 2;
  PartitionedHypergraph phg(k, hg);
  phg.setNodePart(0, 1);
  phg.setNodePart(1, 1);
  for (HypernodeID u = 4; u < 12; ++u) {
    phg.setNodePart(u, 0);
  }
  phg.setNodePart(2, 0);
  phg.setNodePart(3, 0);
  for (HypernodeID u = 12; u < 20; ++u) {
    phg.setNodePart(u, 1);
  }
    phg.initializeGainCache();

  Context context;
  context.partition.k = k;
  context.setupPartWeights(phg.totalWeight());
  context.partition.max_part_weights = { std::numeric_limits<HypernodeWeight>::max(), std::numeric_limits<HypernodeWeight>::max()};
  context.refinement.fm.rollback_balance_violation_factor = 0.0;

  FMSharedData sharedData(hg.initialNumNodes(), context);

  GlobalRollback grb(hg, context, k);
  grb.setRemainingOriginalPins(phg);

  auto performUpdates = [&](Move& m) {
   sharedData.moveTracker.insertMove(m);
  };

  vec<Gain> expected_gains = { 3, 1 };

  ASSERT_EQ(phg.km1Gain(2, 0, 1), 3);
  Move move_2 = { 0, 1, 2, 3 };
    phg.changeNodePartWithGainCacheUpdate(move_2.node, move_2.from, move_2.to);

  ASSERT_EQ(phg.km1Gain(0, 1, 0), 1);
  Move move_0 = { 1, 0, 0, 1 };
    phg.changeNodePartWithGainCacheUpdate(move_0.node, move_0.from, move_0.to);

  performUpdates(move_0);
  performUpdates(move_2);

  grb.recalculateGains(phg, sharedData);
  grb.verifyGains<true>(phg, sharedData);
}

#endif

}   // namespace mt_kahypar
