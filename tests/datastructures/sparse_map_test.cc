/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <atomic>
#include <cstdlib>
#include <mt-kahypar/macros.h>

#include "gmock/gmock.h"
#include "tbb/task_group.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/datastructures/concurrent_flat_map.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

template<typename MapType>
struct ADynamicSparseMap : public Test {
  MapType map;
};

using DynamicSparseMapTestTypes = ::testing::Types<DynamicSparseMap<size_t, size_t>, DynamicFlatMap<size_t, size_t>>;

TYPED_TEST_CASE(ADynamicSparseMap, DynamicSparseMapTestTypes);

TYPED_TEST(ADynamicSparseMap, AddsSeveralElements) {
  auto& map = this->map;
  map.initialize(1);
  map[4] = 5;
  map[8] = 1;
  map[1] = 4;
  ASSERT_EQ(3, map.size());
  ASSERT_EQ(5, map[4]);
  ASSERT_EQ(1, map[8]);
  ASSERT_EQ(4, map[1]);
}

TYPED_TEST(ADynamicSparseMap, ModifiesAnExistingValue) {
  auto& map = this->map;
  map.initialize(1);
  map[4] = 5;
  map[8] = 1;
  map[1] = 4;
  ++map[1];
  ASSERT_EQ(3, map.size());
  ASSERT_EQ(5, map[4]);
  ASSERT_EQ(1, map[8]);
  ASSERT_EQ(5, map[1]);
}

TYPED_TEST(ADynamicSparseMap, IsForcedToGrow) {
  const size_t initial_capacity = 256;
  auto& map = this->map;
  map.initialize(initial_capacity);
  const size_t n = map.capacity();
  for ( size_t i = 0; i < (2 * n) / 5; ++i ) {
    map[i] = i;
  }
  ASSERT_EQ(initial_capacity, map.capacity());
  ASSERT_EQ((2 * n) / 5, map.size());

  // Forces map to dynamically grow
  map[n] = n;

  ASSERT_EQ(2 * initial_capacity, map.capacity());
  ASSERT_EQ((2 * n) / 5 + 1, map.size());
  ASSERT_EQ(n, map[n]++);
  for ( size_t i = 0; i < (2 * n) / 5; ++i ) {
    ASSERT_EQ(i, map[i]++);
  }

  ASSERT_EQ(n + 1, map[n]);
  for ( size_t i = 0; i < (2 * n) / 5; ++i ) {
    ASSERT_EQ(i + 1, map[i]);
  }
}

}  // namespace ds
}  // namespace mt_kahypar
