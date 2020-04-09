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

#include <numeric>
#include <algorithm>

#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/vector.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(AVector, WritesAnValueToStrippedVector1) {
  Vector<int> vec(256, 0);
  vec[0] = 31;
  ASSERT_EQ(31, vec[0]);
}

TEST(AVector, WritesAnValueToStrippedVector2) {
  Vector<int> vec(256, 0);
  vec[65] = 35;
  ASSERT_EQ(35, vec[65]);
}

TEST(AVector, WritesAnValueToStrippedVector3) {
  Vector<int> vec(256, 0);
  vec[127] = 42;
  ASSERT_EQ(42, vec[127]);
}

TEST(AVector, WritesAnValueToStrippedVector4) {
  Vector<int> vec(256, 0);
  vec[128] = 43;
  ASSERT_EQ(43, vec[128]);
}

TEST(AVector, WritesValuesToWholeVector) {
  Vector<int> vec(256, 0);

  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] = i;
  }
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}

TEST(AVector, IsInitializedWithNonDefaultValues) {
  Vector<int> vec(256, 42);

  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
}

TEST(AVector, AssignValuesToAlreadyInitializedVector1) {
  Vector<int> vec(256, 0);

  const size_t count = 31;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AVector, AssignValuesToAlreadyInitializedVector2) {
  Vector<int> vec(256, 0);

  const size_t count = 42;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AVector, AssignValuesToAlreadyInitializedVector3) {
  Vector<int> vec(256, 0);

  const size_t count = 127;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AVector, AssignValuesToAlreadyInitializedVector4) {
  Vector<int> vec(256, 0);

  const size_t count = 128;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AVector, AssignValuesToAlreadyInitializedVector5) {
  Vector<int> vec(256, 0);

  const size_t count = 256;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AVector, FilledWithNumbersFromZeroToN) {
  Vector<int> vec(256, 0);
  std::iota(vec.begin(), vec.end(), 0);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}

TEST(AVector, ChecksDistanceBetweenTwoPointers1) {
  Vector<int> vec(256, 0);
  ASSERT_EQ(5, std::distance(vec.begin(), vec.begin() + 5));
}

TEST(AVector, ChecksDistanceBetweenTwoPointers2) {
  Vector<int> vec(256, 0);
  ASSERT_EQ(42, std::distance(vec.begin() + 24, vec.begin() + 66));
}

TEST(AVector, ChecksDistanceBetweenTwoPointers3) {
  Vector<int> vec(256, 0);
  ASSERT_EQ(256, std::distance(vec.begin(), vec.end()));
}

TEST(AVector, CanBeSorted) {
  Vector<int> vec(256, 0);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] = (vec.size() - 1) - i;
  }
  std::sort(vec.begin(), vec.end());
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}

TEST(AVector, MemcopiesContentToVector) {
  Vector<int> vec(256, 0);
  std::vector<int> vec2(256, 5);
  memcpy(vec.data(), vec2.data(), sizeof(int) * 256);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(5, vec[i]);
  }
}

}  // namespace ds
}  // namespace mt_kahypar
