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

#include "mt-kahypar/parallel/stl/scalable_vector.h"
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

TEST(AVector, IsInitializedWithMemoryChunkFromMemoryPool) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Vector<size_t> vec("TEST_GROUP", "TEST_CHUNK", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AVector, IsInitializedWithSeveralMemoryChunksFromMemoryPool) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK_1", 5, sizeof(size_t));
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK_2", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Vector<size_t> vec_1("TEST_GROUP", "TEST_CHUNK_1", 5);
  Vector<size_t> vec_2("TEST_GROUP", "TEST_CHUNK_2", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK_1"), (char *) vec_1.data());
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK_2"), (char *) vec_2.data());
  ASSERT_NE(vec_1.data(), vec_2.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AVector, ReleasesMemoryChunkFromMemoryPoolInDestructor) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  {
    Vector<size_t> vec("TEST_GROUP", "TEST_CHUNK", 5);
    ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());
    ASSERT_EQ(nullptr, parallel::MemoryPool::instance().request_mem_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t)));
  }

  ASSERT_NE(nullptr, parallel::MemoryPool::instance().request_mem_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t)));

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AVector, AllocatesOwnMemoryIfNotAvailableInMemoryPool) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Vector<size_t> vec("TEST_GROUP", "OTHER_CHUNK", 5);
  ASSERT_NE(nullptr, vec.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AVector, AllocatesOwnMemoryOnOverAllocationInMemoryPool) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Vector<size_t> vec("TEST_GROUP", "TEST_CHUNK", 10);
  ASSERT_NE(nullptr, vec.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AVector, AllocatesOwnMemoryIfAlreadyRequestedInMemoryPool) {
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Vector<size_t> vec_1("TEST_GROUP", "TEST_CHUNK", 5);
  Vector<size_t> vec_2("TEST_GROUP", "TEST_CHUNK", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec_1.data());
  ASSERT_NE(nullptr, vec_2.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec_2.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

}  // namespace ds
}  // namespace mt_kahypar
