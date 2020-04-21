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

#include "gmock/gmock.h"
#include "tbb/task_group.h"

#include "mt-kahypar/parallel/memory_pool.h"

using ::testing::Test;

namespace mt_kahypar {
namespace parallel {

template <class F, class K>
void executeConcurrent(F f1, K f2) {
  std::atomic<int> cnt(0);
  tbb::task_group group;

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f1();
      });

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f2();
      });

  group.wait();
}

static void setupMemoryPool() {
  MemoryPool::instance().register_memory_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t));
  MemoryPool::instance().register_memory_chunk("TEST_GROUP_1", "TEST_CHUNK_2", 5, sizeof(int));
  MemoryPool::instance().register_memory_chunk("TEST_GROUP_2", "TEST_CHUNK_1", 5, sizeof(double));
  MemoryPool::instance().register_memory_chunk("TEST_GROUP_2", "TEST_CHUNK_2", 5, sizeof(size_t));
  MemoryPool::instance().allocate_memory_chunks();
}

TEST(AMemoryPool, AllocatesMemory) {
 setupMemoryPool();

  ASSERT_NE(nullptr, MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"));
  ASSERT_NE(nullptr, MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_2"));
  ASSERT_NE(nullptr, MemoryPool::instance().mem_chunk("TEST_GROUP_2", "TEST_CHUNK_1"));
  ASSERT_NE(nullptr, MemoryPool::instance().mem_chunk("TEST_GROUP_2", "TEST_CHUNK_2"));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, CanHandleARequest) {
 setupMemoryPool();

  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, CanHandleSeveralRequests) {
 setupMemoryPool();

  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));
  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_2", "TEST_CHUNK_2"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_2", "TEST_CHUNK_2", 5, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, ReturnsNullptrIfMemoryChunkAlreadyRequested) {
 setupMemoryPool();

  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));
  ASSERT_EQ(nullptr,
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}


TEST(AMemoryPool, ReturnsNullptrOnOverallocation) {
 setupMemoryPool();

  ASSERT_EQ(nullptr,
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 10, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, ReturnsNullptrIfMemoryChunkIsNotAvailable) {
 setupMemoryPool();

  ASSERT_EQ(nullptr,
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_3", 5, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, RequestMemoryAfterRelease) {
 setupMemoryPool();

  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));
  MemoryPool::instance().release_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1");
  ASSERT_EQ(MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1"),
            MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t)));

  MemoryPool::instance().free_memory_chunks();
}

TEST(AMemoryPool, OnlyOneRequestSucceedsOnConcurrentAccess) {
 setupMemoryPool();

  char* chunk_1 = nullptr;
  char* chunk_2 = nullptr;
  executeConcurrent([&] {
    chunk_1 = MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t));
  }, [&] {
    chunk_2 = MemoryPool::instance().request_mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1", 5, sizeof(size_t));
  });

  char* expected_chunk = MemoryPool::instance().mem_chunk("TEST_GROUP_1", "TEST_CHUNK_1");
  if ( chunk_1 ) {
    ASSERT_EQ(nullptr, chunk_2);
    ASSERT_EQ(expected_chunk, chunk_1);
  } else {
    ASSERT_EQ(nullptr, chunk_1);
    ASSERT_EQ(expected_chunk, chunk_2);
  }

  MemoryPool::instance().free_memory_chunks();
}


}  // namespace parallel
}  // namespace mt_kahypar
