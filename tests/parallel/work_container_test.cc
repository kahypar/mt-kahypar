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
  WorkContainer<int> cdc(n);
  tbb::parallel_for(0, m, [&](int i) {
    cdc.template push_back<true>(i);
  });
  ASSERT_EQ(cdc.unsafe_size(), m);

  tbb::enumerable_thread_specific<int> counters;
  tbb::task_group tg;
  int num_tasks = 7;
  for (int i = 0; i < num_tasks; ++i) {
    tg.run([&]() {
      int res = 0;
      int& lc = counters.local();
      while (cdc.try_pop(res)) {
        lc++;
      }
    });
  }
  tg.wait();

  int overall = counters.combine(std::plus<int>());
  ASSERT_EQ(overall, m);
  ASSERT_TRUE(cdc.unsafe_size() == 0);
}

TEST(WorkContainer, ClearWorks) {
  WorkContainer<int> cdc(n);
  cdc.template push_back<true>(5);
  cdc.template push_back<true>(420);
  ASSERT_EQ(cdc.unsafe_size(), 2);
  cdc.clear();
  ASSERT_TRUE(cdc.unsafe_size() == 0);
}

TEST(WorkContainer, WorkStealingWorks) {
  WorkContainer<int> cdc(n);

  std::atomic<size_t> stage { 0 };
  size_t steals = 0;
  size_t own_pops = 0;

  int m = 99999;

  std::thread producer([&] {
    for (int i = 0; i < m; ++i) {
      cdc.template push_back<false>(i);
    }

    stage.fetch_add(1, std::memory_order_acq_rel);

    int own_element;
    while (cdc.try_pop(own_element)) {
      own_pops++;
    }
  });

  std::thread consumer([&] {
    while (stage.load(std::memory_order_acq_rel) < 1) { } //spin

    int stolen_element;
    while (cdc.try_pop(stolen_element)) {
      steals++;
    }
  });

  consumer.join();
  producer.join();

  ASSERT_GE(steals, 1);   // this can fail. but it is unlikely --> if it fails for you, just remove it
  ASSERT_EQ(steals + own_pops, m);
  ASSERT_EQ(cdc.unsafe_size(), 0);
}


TEST(WorkContainer, QueueBlocksOnReallocation) {
  SPMCQueue<int> q;

  std::atomic<size_t> stage { 0 };

  std::thread producer([&] {
    for (int i = 0; i < (1 << 13); ++i) {
      q.template push_back<true>(i);
    }
    ASSERT_TRUE(q.unsafe_size() == (1 << 13));
    ASSERT_TRUE(q.next_push_causes_reallocation());

    // this one causes the reallocation
    stage.fetch_add(1, std::memory_order_acq_rel);
    q.template push_back<false>(420);
    stage.fetch_add(1, std::memory_order_acq_rel); // races the consumer thread for getting stage 2
  });

  std::thread consumer([&] {
    while (stage.load(std::memory_order_acq_rel) < 1) { } //spin
    int front_element;
    const bool hit_the_realloc = !q.try_pop_front(front_element);
    size_t s = stage.fetch_add(1, std::memory_order_acq_rel); // races the producer thread for getting stage 2
    if (s == 1) {
      ASSERT_TRUE(hit_the_realloc);
    } else {
      LOG << "could not verify whether queue blocked. realloc was too fast. try again if you'd like";
    }
  });

  consumer.join();
  producer.join();
}


}  // namespace parallel
}  // namespace mt_kahypar
