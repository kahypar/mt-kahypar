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

#include <mt-kahypar/partition/refinement/fm/global_rollback.h>
#include <mt-kahypar/io/hypergraph_io.h>

using ::testing::Test;

namespace mt_kahypar {
namespace refinement {
/*
TEST(RollbackTests, FindsBestPrefix) {

  vec<Gain> gains = { -42, 5, 4, -20, 1, 99, -100, 50 };
  BestIndexReduceBody b(gains);
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b, tbb::static_partitioner());
  ASSERT_EQ(b.best_sum,  -42 + 5 + 4 - 20 + 1 + 99);
  ASSERT_EQ(b.sum, -42 + 5 + 4 - 20 + 1 + 99 - 100 + 50);
  ASSERT_EQ(b.best_index, 5);
}

TEST(RollbackTests, FindsBestPrefixLargeRandom) {
  bool display_timing = false;

  tbb::task_scheduler_init tsi(4);
  size_t n = 1000 * 100;
  vec<Gain> gains(n, 0);
  std::mt19937 rng(420);
  std::uniform_int_distribution<Gain> distr(-5000, 5000);

  auto start_init = tbb::tick_count::now();
  for (MoveID i = 0; i < n; ++i) {
    gains[i] = distr(rng);
  }
  if (display_timing) LOG << "Finish init in " << (tbb::tick_count::now() - start_init).seconds() << "seconds";

  auto start_reduce_sequential = tbb::tick_count::now();
  Gain sum = 0, best_sum = 0;
  MoveID best_index = 0;
  for (MoveID i = 0; i < n; ++i) {
    sum += gains[i];
    if (sum > best_sum) {
      best_sum = sum;
      best_index = i;
    }
  }
  if (display_timing) LOG << "Finish sequential  reduce in " << (tbb::tick_count::now() - start_reduce_sequential).seconds() << "seconds";

  auto start_reduce = tbb::tick_count::now();
  BestIndexReduceBody b(gains);
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b);
  if (display_timing) LOG << "Finish reduce in " << (tbb::tick_count::now() - start_reduce).seconds() << "seconds";

  auto start_reduce_static = tbb::tick_count::now();
  BestIndexReduceBody b2(gains);
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b2, tbb::static_partitioner());
  if (display_timing) LOG << "Finish reduce with static partitioner in " << (tbb::tick_count::now() - start_reduce_static).seconds() << "seconds";

  ASSERT_EQ(best_sum, b.best_sum);
  ASSERT_EQ(sum, b.sum);
  ASSERT_EQ(best_index, b.best_index);
}
*/


TEST(RollbackTests, GainRecalculationAndRollsbackCorrectly) {
  Hypergraph hg = io::readHypergraphFile<Hypergraph, HypergraphFactory>("../test_instances/twocenters.hgr", 0);
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
  phg.initializeGainInformation();
  FMSharedData sharedData(hg.initialNumNodes(), hg.initialNumEdges(), k, 4);
  sharedData.setRemainingOriginalPins(phg);

  auto performMove = [&](Move m) {
    phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(), []{});
    MoveID move_id = sharedData.moveTracker.insertMove(m);
    for (HyperedgeID e : phg.incidentEdges(m.node)) {
      sharedData.performHyperedgeSpecificMoveUpdates(m, move_id, e);
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

  GlobalRollBack grb(hg.initialNumNodes());
  grb.recalculateGains(phg, sharedData);
  for (MoveID round_local_move_id = 0; round_local_move_id < 4; ++round_local_move_id) {
    ASSERT_EQ(sharedData.moveTracker.globalMoveOrder[round_local_move_id].gain, grb.gains[round_local_move_id]);
  }

  ASSERT_EQ(phg.km1Gain(4, 0, 1), -1);
  performMove({0, 1, 4, -1});
  ASSERT_EQ(phg.km1Gain(5, 0, 1), 0);
  performMove({0, 1, 5, 0});

  grb.globalRollbackToBestPrefix(phg, sharedData);
  // revert last two moves
  ASSERT_EQ(phg.partID(4), 0);
  ASSERT_EQ(phg.partID(5), 0);
}


TEST(RollbackTests, GainRecalculation2) {
  Hypergraph hg = io::readHypergraphFile<Hypergraph, HypergraphFactory>("../test_instances/twocenters.hgr", 0);
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
  phg.initializeGainInformation();
  FMSharedData sharedData(hg.initialNumNodes(), hg.initialNumEdges(), k, 4);
  sharedData.setRemainingOriginalPins(phg);

  auto performUpdates = [&](Move& m) {
    MoveID move_id = sharedData.moveTracker.insertMove(m);
    for ( HyperedgeID e : phg.incidentEdges(m.node) ){
      sharedData.performHyperedgeSpecificMoveUpdates(m, move_id, e);
    }
  };

  vec<Gain> expected_gains = { 3, 1 };

  ASSERT_EQ(phg.km1Gain(2, 0, 1), 3);
  Move move_2 = { 0, 1, 2, 3 };
  phg.changeNodePartFullUpdate(move_2.node, move_2.from, move_2.to, std::numeric_limits<HypernodeWeight>::max(), []{});

  ASSERT_EQ(phg.km1Gain(0, 1, 0), 1);
  Move move_0 = { 1, 0, 0, 1 };
  phg.changeNodePartFullUpdate(move_0.node, move_0.from, move_0.to, std::numeric_limits<HypernodeWeight>::max(), []{});

  performUpdates(move_0);
  performUpdates(move_2);

  GlobalRollBack grb(hg.initialNumNodes());
  grb.recalculateGains(phg, sharedData);
  for (MoveID round_local_move_id = 0; round_local_move_id < sharedData.moveTracker.numPerformedMoves(); ++round_local_move_id) {
    ASSERT_EQ(grb.gains[round_local_move_id], expected_gains[round_local_move_id]);
  }

}

}   // namespace refinement
}   // namespace mt_kahypar
