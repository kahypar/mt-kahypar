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

#include "gmock/gmock.h"

#include <mt-kahypar/parallel/work_stack.h>
#include <thread>

using ::testing::Test;

namespace mt_kahypar {
namespace parallel {

size_t n = 100000;


TEST(WorkContainer, HasCorrectSizeAfterParallelInsertionAndDeletion) {
  int m = 75000;
  WorkContainer<int> cdc(std::thread::hardware_concurrency());
  tbb::parallel_for(0, m, [&](int i) {
    cdc.safe_push(i, tbb::this_task_arena::current_thread_index());
  });
  ASSERT_EQ(cdc.unsafe_size(), m);

  tbb::enumerable_thread_specific<int> counters;
  tbb::task_group tg;
  int num_tasks = 7;
  for (int i = 0; i < num_tasks; ++i) {
    tg.run([&]() {
      int res = 0;
      int& lc = counters.local();
      while (cdc.try_pop(res, tbb::this_task_arena::current_thread_index())) {
        lc++;
      }
    });
  }
  tg.wait();

  int overall = counters.combine(std::plus<int>());
  ASSERT_EQ(overall, m);

  cdc.clear();
  ASSERT_EQ(cdc.unsafe_size(), 0);
}

TEST(WorkContainer, ClearWorks) {
  WorkContainer<int> cdc(std::thread::hardware_concurrency());
  cdc.safe_push(5, tbb::this_task_arena::current_thread_index());
  cdc.safe_push(420, tbb::this_task_arena::current_thread_index());
  ASSERT_EQ(cdc.unsafe_size(), 2);
  cdc.clear();
  ASSERT_TRUE(cdc.unsafe_size() == 0);
}


TEST(WorkContainer, WorkStealingWorks) {
  WorkContainer<int> cdc(std::thread::hardware_concurrency());

  std::atomic<size_t> stage { 0 };
  size_t steals = 0;
  size_t own_pops = 0;

  int m = 99999;

  std::thread producer([&] {
    int thread_id = 0;
    for (int i = 0; i < m; ++i) {
      cdc.safe_push(i, thread_id);
    }

    stage.fetch_add(1, std::memory_order_acq_rel);

    int own_element;
    while (cdc.try_pop(own_element, thread_id)) {
      own_pops++;
    }
  });

  std::thread consumer([&] {
    int thread_id = 1;
    while (stage.load(std::memory_order_acq_rel) < 1) { } //spin

    int stolen_element;
    while (cdc.try_pop(stolen_element, thread_id)) {
      steals++;
    }
  });

  consumer.join();
  producer.join();

  ASSERT_EQ(steals + own_pops, m);
}

}  // namespace parallel
}  // namespace mt_kahypar
