/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>
#include <cstdlib>
#include <mt-kahypar/macros.h>

#include "gmock/gmock.h"
#include "tbb/task_group.h"

#include "mt-kahypar/datastructures/sparse_map.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(ADynamicSparseMap, AddsSeveralElements) {
  DynamicSparseMap<size_t, size_t> map;
  map[4] = 5;
  map[8] = 1;
  map[1] = 4;
  ASSERT_EQ(3, map.size());
  ASSERT_EQ(5, map[4]);
  ASSERT_EQ(1, map[8]);
  ASSERT_EQ(4, map[1]);
}

TEST(ADynamicSparseMap, ModifiesAnExistingValue) {
  DynamicSparseMap<size_t, size_t> map;
  map[4] = 5;
  map[8] = 1;
  map[1] = 4;
  ++map[1];
  ASSERT_EQ(3, map.size());
  ASSERT_EQ(5, map[4]);
  ASSERT_EQ(1, map[8]);
  ASSERT_EQ(5, map[1]);
}

TEST(ADynamicSparseMap, IsForcedToGrow) {
  DynamicSparseMap<size_t, size_t> map;
  const size_t n = map.capacity();
  for ( size_t i = 0; i < n / 3; ++i ) {
    map[i] = i;
  }
  const size_t initial_capacity = DynamicSparseMap<size_t, size_t>::MAP_SIZE;
  ASSERT_EQ(initial_capacity, map.capacity());
  ASSERT_EQ(n / 3, map.size());

  // Forces map to dynamically grow
  map[n] = n;

  ASSERT_EQ(2 * initial_capacity, map.capacity());
  ASSERT_EQ(n / 3 + 1, map.size());
  ASSERT_EQ(n, map[n]++);
  for ( size_t i = 0; i < n / 3; ++i ) {
    ASSERT_EQ(i, map[i]++);
  }

  ASSERT_EQ(n + 1, map[n]);
  for ( size_t i = 0; i < n / 3; ++i ) {
    ASSERT_EQ(i + 1, map[i]);
  }
}

}  // namespace ds
}  // namespace mt_kahypar
