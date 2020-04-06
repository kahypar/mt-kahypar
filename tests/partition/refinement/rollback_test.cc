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


using ::testing::Test;

namespace mt_kahypar {
namespace refinement {

TEST(RollBackTests, FindsBestPrefix) {

  vec<Gain> gains = { -42, 5, 4, -20, 1, 99, -100, 50 };
  BestIndexReduceBody b(gains);
  tbb::parallel_reduce(tbb::blocked_range<MoveID>(0, gains.size()), b, tbb::static_partitioner());
  ASSERT_EQ(b.best_sum,  -42 + 5 + 4 - 20 + 1 + 99);
  ASSERT_EQ(b.sum, -42 + 5 + 4 - 20 + 1 + 99 - 100 + 50);
  ASSERT_EQ(b.best_index, 5);
}

}   // namespace refinement
}   // namespace mt_kahypar
