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

using ::testing::Test;

namespace mt_kahypar {
namespace parallel {

size_t n = 100000;
size_t nSockets = 9;

TEST(NumaWorkQueueContainer, HasCorrectSizeAfterParallelInsertion) {
  int m = 75000;
  WorkStack<int> cdc(n);
  tbb::parallel_for(0, m, [&](int i) {
    cdc.push_back(i);
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
  ASSERT_TRUE(cdc.unsafe_size() == static_cast<size_t>(-num_tasks));
  ASSERT_TRUE(cdc.empty());
}

TEST(NumaWorkQueueContainer, ClearWorks) {
  WorkStack<int> cdc(n);
  cdc.push_back(5);
  cdc.push_back(420);
  ASSERT_EQ(cdc.unsafe_size(), 2);
  cdc.clear();
  ASSERT_TRUE(cdc.empty());
  cdc.shrink_to_fit();
  ASSERT_EQ(cdc.capacity(), 0);
}


/*TEST(NumaWorkQueue, WorkStealingWorks) {
  NumaWorkQueue<int> wq(nSockets, n);
  wq.push(5, 2);
  wq.push(4, 0);
  wq.push(919, 0);
  wq.push(9, 1);
  wq.push(1613, 2);


  int res = -1;
  ASSERT_TRUE(wq.tryPop(res, 0));
  ASSERT_EQ(res, 919);

  ASSERT_TRUE(wq.tryPop(res, 0));
  ASSERT_EQ(res, 4);

  ASSERT_TRUE(wq.tryPop(res, 0));  // this one should steal from socket 2
  ASSERT_EQ(res, 1613);

  ASSERT_TRUE(wq.tryPop(res, 0));  // steal from socket 1
  ASSERT_EQ(res, 9);

  ASSERT_TRUE(wq.tryPop(res, 0));  // steal again from socket 2
  ASSERT_EQ(res, 5);

  ASSERT_FALSE(wq.tryPop(res, 0));
}*/


}  // namespace parallel
}  // namespace mt_kahypar
