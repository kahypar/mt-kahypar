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
#include <mt-kahypar/definitions.h>
#include <mt-kahypar/io/hypergraph_io.h>

using ::testing::Test;

namespace mt_kahypar {
namespace refinement {

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
  size_t n = 1000 * 1000 * 2;
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
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b);//, tbb::static_partitioner());
  if (display_timing) LOG << "Finish reduce in " << (tbb::tick_count::now() - start_reduce).seconds() << "seconds";

  auto start_reduce_static = tbb::tick_count::now();
  BestIndexReduceBody b2(gains);
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b2, tbb::static_partitioner());
  if (display_timing) LOG << "Finish reduce with static partitioner in " << (tbb::tick_count::now() - start_reduce_static).seconds() << "seconds";

  ASSERT_EQ(best_sum, b.best_sum);
  ASSERT_EQ(sum, b.sum);
  ASSERT_EQ(best_index, b.best_index);
}

TEST(RollbackTests, GainRecalculation) {
  Hypergraph hg = io::readHypergraphFile<Hypergraph, HypergraphFactory>("../test_instances/twocenters.hgr", 0);
  PartitionID k = 2;
  PartitionedHypergraph phg(k, hg);

  FMSharedData sharedData(hg.initialNumNodes(), hg.initialNumEdges(), k, 4);



}

}   // namespace refinement
}   // namespace mt_kahypar
