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

#include "mt-kahypar/datastructures/connectivity_set.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {
using PartitionID = int32_t;

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

void add(ConnectivitySet& conn_set, const std::set<PartitionID>& ids) {
  for (const PartitionID& id : ids) {
    conn_set.add(id);
  }
}

void remove(ConnectivitySet& conn_set, const std::set<PartitionID>& ids) {
  for (const PartitionID& id : ids) {
    conn_set.remove(id);
  }
}

void verify(const ConnectivitySet& conn_set,
            const PartitionID k,
            const std::set<PartitionID>& contained) {
  // Verify bitset in connectivity set
  ASSERT_EQ(contained.size(), conn_set.size());
  for (PartitionID i = 0; i < k; ++i) {
    if (contained.find(i) != contained.end()) {
      ASSERT_TRUE(conn_set.contains(i)) << V(i);
    } else {
      ASSERT_FALSE(conn_set.contains(i)) << V(i);
    }
  }

  // Verify iterator
  size_t connectivity = 0;
  for (const PartitionID id : conn_set.connectivitySet()) {
    ASSERT_TRUE(contained.find(id) != contained.end()) << V(id);
    ++connectivity;
  }
  ASSERT_EQ(contained.size(), connectivity);
}

TEST(AConnectivitySet, IsCorrectInitialized) {
  ConnectivitySet conn_set(32);
  verify(conn_set, 32, { });
}

TEST(AConnectivitySet, AddOnePartition1) {
  ConnectivitySet conn_set(32);
  conn_set.add(2);
  verify(conn_set, 32, { 2 });
}

TEST(AConnectivitySet, AddOnePartition2) {
  ConnectivitySet conn_set(32);
  conn_set.add(14);
  verify(conn_set, 32, { 14 });
}

TEST(AConnectivitySet, AddOnePartition3) {
  ConnectivitySet conn_set(32);
  conn_set.add(23);
  verify(conn_set, 32, { 23 });
}

TEST(AConnectivitySet, AddOnePartition4) {
  ConnectivitySet conn_set(32);
  conn_set.add(30);
  verify(conn_set, 32, { 30 });
}

TEST(AConnectivitySet, AddTwoPartitions1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 5, 31 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitions2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 14, 24 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitions3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 7, 16 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitions1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 0, 1, 5, 14, 24, 27, 31 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitions2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 5, 6, 7, 11, 13, 15, 24, 28, 30 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitions3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
  add(conn_set, added);
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitionsAndRemoveOne1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 5, 31 };
  std::set<PartitionID> removed = { 31 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitionsAndRemoveOne2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 16, 17 };
  std::set<PartitionID> removed = { 16 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitionsAndRemoveOne3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 7, 21 };
  std::set<PartitionID> removed = { 7 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddTwoPartitionsAndRemoveOne4) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 25, 27 };
  std::set<PartitionID> removed = { 27 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitionsAndRemoveSeveral1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 1, 13, 15, 23, 24, 30 };
  std::set<PartitionID> removed = { 13, 15, 23 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitionsAndRemoveSeveral2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 2, 5, 6, 14, 15, 21, 23, 29 };
  std::set<PartitionID> removed = { 5, 14, 21, 29 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddSeveralPartitionsAndRemoveSeveral3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 0, 1, 2, 3, 4, 5, 6, 7, 24, 25, 26, 27, 28, 29 };
  std::set<PartitionID> removed = { 5, 6, 7, 24, 25, 26, 27 };
  add(conn_set, added);
  remove(conn_set, removed);
  for (const PartitionID id : removed) {
    added.erase(id);
  }
  verify(conn_set, 32, added);
}


TEST(AConnectivitySet, AddConcurrentTwoPartitions1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 7, 16 };
  executeConcurrent([&] {
        conn_set.add(7);
      }, [&] {
        conn_set.add(16);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddConcurrentTwoPartitions2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 4, 5 };
  executeConcurrent([&] {
        conn_set.add(4);
      }, [&] {
        conn_set.add(5);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddConcurrentTwoPartitions3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 12, 14 };
  executeConcurrent([&] {
        conn_set.add(12);
      }, [&] {
        conn_set.add(14);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddConcurrentTwoPartitions4) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 30, 31 };
  executeConcurrent([&] {
        conn_set.add(30);
      }, [&] {
        conn_set.add(31);
      });
  verify(conn_set, 32, added);
}


TEST(AConnectivitySet, AddConcurrentSeveralPartitions1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 1, 13, 15, 23, 24, 30 };
  executeConcurrent([&] {
        add(conn_set, { 1, 15, 24 });
      }, [&] {
        add(conn_set, { 13, 23, 30 });
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddConcurrentSeveralPartitions2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 0, 4, 15, 23, 24, 25, 28, 30, 31 };
  executeConcurrent([&] {
        add(conn_set, { 0, 15, 24, 28, 31 });
      }, [&] {
        add(conn_set, { 4, 23, 25, 30 });
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddConcurrentSeveralPartitions3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 6, 7, 8, 9, 10, 14, 15, 16, 17, 23, 24, 25, 26 };
  executeConcurrent([&] {
        add(conn_set, { 6, 8, 10, 15, 17, 24, 26 });
      }, [&] {
        add(conn_set, { 7, 9, 14, 16, 23, 25 });
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveOnePartitionConcurrently1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 3 };
  conn_set.add(4);
  executeConcurrent([&] {
        conn_set.add(3);
      }, [&] {
        conn_set.remove(4);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveOnePartitionConcurrently2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 24 };
  conn_set.add(9);
  executeConcurrent([&] {
        conn_set.add(24);
      }, [&] {
        conn_set.remove(9);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveOnePartitionConcurrently3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 15 };
  conn_set.add(31);
  executeConcurrent([&] {
        conn_set.add(15);
      }, [&] {
        conn_set.remove(31);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveOnePartitionConcurrently4) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 13 };
  conn_set.add(14);
  executeConcurrent([&] {
        conn_set.add(13);
      }, [&] {
        conn_set.remove(14);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveSeveralPartitionsConcurrently1) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 1, 4, 5, 15, 16, 21, 30 };
  add(conn_set, { 2, 3, 14, 17, 23, 24 });
  executeConcurrent([&] {
        conn_set.add(1);
        conn_set.remove(3);
        conn_set.add(5);
        conn_set.remove(14);
        conn_set.add(16);
        conn_set.remove(23);
        conn_set.add(30);
      }, [&] {
        conn_set.remove(2);
        conn_set.add(4);
        conn_set.remove(17);
        conn_set.add(15);
        conn_set.remove(24);
        conn_set.add(21);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveSeveralPartitionsConcurrently2) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 4, 5, 8, 10, 18, 21, 22, 27, 31 };
  add(conn_set, { 1, 2, 5, 17, 19, 20, 23, 30 });
  executeConcurrent([&] {
        conn_set.add(4);
        conn_set.remove(1);
        conn_set.add(8);
        conn_set.remove(5);
        conn_set.add(18);
        conn_set.remove(19);
        conn_set.add(22);
        conn_set.remove(23);
        conn_set.add(31);
      }, [&] {
        conn_set.remove(2);
        conn_set.add(5);
        conn_set.remove(17);
        conn_set.add(10);
        conn_set.remove(20);
        conn_set.add(21);
        conn_set.remove(30);
        conn_set.add(27);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveSeveralPartitionsConcurrently3) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 13, 14, 15, 17, 18, 22, 23, 24, 25, 26, 27 };
  add(conn_set, { 6, 7, 8, 9, 10, 19, 20, 21, 28, 29, 30 });
  executeConcurrent([&] {
        conn_set.add(13);
        conn_set.remove(6);
        conn_set.add(15);
        conn_set.remove(8);
        conn_set.add(18);
        conn_set.remove(10);
        conn_set.add(23);
        conn_set.remove(20);
        conn_set.add(25);
        conn_set.remove(28);
        conn_set.add(27);
        conn_set.remove(30);
      }, [&] {
        conn_set.remove(7);
        conn_set.add(14);
        conn_set.remove(9);
        conn_set.add(17);
        conn_set.remove(19);
        conn_set.add(22);
        conn_set.remove(21);
        conn_set.add(24);
        conn_set.remove(29);
        conn_set.add(26);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently1) {
  ConnectivitySet conn_set(32);
  executeConcurrent([&] {
        conn_set.add(15);
      }, [&] {
        conn_set.remove(15);
      });
  verify(conn_set, 32, { });
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently2) {
  ConnectivitySet conn_set(32);
  conn_set.add(15);
  executeConcurrent([&] {
        conn_set.add(15);
      }, [&] {
        conn_set.remove(15);
      });
  verify(conn_set, 32, { 15 });
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently3) {
  ConnectivitySet conn_set(32);
  executeConcurrent([&] {
        conn_set.add(6);
      }, [&] {
        conn_set.remove(6);
      });
  verify(conn_set, 32, { });
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently4) {
  ConnectivitySet conn_set(32);
  conn_set.add(6);
  executeConcurrent([&] {
        conn_set.add(6);
      }, [&] {
        conn_set.remove(6);
      });
  verify(conn_set, 32, { 6 });
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently5) {
  ConnectivitySet conn_set(32);
  executeConcurrent([&] {
        conn_set.add(22);
      }, [&] {
        conn_set.remove(22);
      });
  verify(conn_set, 32, { });
}

TEST(AConnectivitySet, AddAndRemoveSamePartitionConcurrently6) {
  ConnectivitySet conn_set(32);
  conn_set.add(22);
  executeConcurrent([&] {
        conn_set.add(22);
      }, [&] {
        conn_set.remove(22);
      });
  verify(conn_set, 32, { 22 });
}

TEST(AConnectivitySet, AddAndRemoveSeveralSamePartitionsConcurrently) {
  ConnectivitySet conn_set(32);
  std::set<PartitionID> added = { 13, 14, 15, 16, 18, 22, 23 };
  add(conn_set, { 13, 14, 15, 17, 18, 22, 23, 24, 25, 26, 27 });
  executeConcurrent([&] {
        conn_set.add(16);
        conn_set.remove(13);
        conn_set.add(24);
        conn_set.remove(24);
        conn_set.remove(26);
        conn_set.remove(17);
      }, [&] {
        conn_set.add(13);
        conn_set.remove(16);
        conn_set.add(16);
        conn_set.remove(24);
        conn_set.remove(25);
        conn_set.add(26);
        conn_set.remove(26);
        conn_set.remove(27);
      });
  verify(conn_set, 32, added);
}

TEST(AConnectivitySet, IteratesThroughPartitionsAndSimultanouslyAddElements) {
  ConnectivitySet conn_set(32);
  add(conn_set, { 1, 2, 6, 10, 15, 22, 24, 31 });

  std::atomic<size_t> cnt(0);

  executeConcurrent([&] {
        while (cnt < 1);
        // Add 5 -> since iterator already traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be not visible in the iterator
        conn_set.add(5);
        ++cnt;

        while (cnt < 3);
        // Add 12 -> since iterator currently traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be visible in the current iterator
        conn_set.add(12);
        ++cnt;

        while (cnt < 5);
        // Add 14 -> since iterator already traverse
        // bitset containing partition ids 8 .. 16,
        // the change will be not visible in the iterator
        conn_set.add(14);
        ++cnt;

        while (cnt < 7);
        // Add 30 -> since iterator currently traverse
        // bitset containing partition ids 16 .. 24,
        // the change will be visible in the current iterator
        conn_set.add(30);
        ++cnt;
      }, [&] {
        std::vector<PartitionID> expected = { 1, 2, 6, 10, 12, 15, 22, 24, 30, 31 };
        size_t i = 0;
        for (const PartitionID& id : conn_set.connectivitySet()) {
          if (i == 1) {
            ++cnt;
            while (cnt < 2);
            // Wait until thread 1 adds 5
          }

          if (i == 2) {
            ++cnt;
            while (cnt < 4);
            // Wait until thread 1 adds 12
          }

          if (i == 4) {
            ++cnt;
            while (cnt < 6);
            // Wait until thread 1 adds 14
          }

          if (i == 6) {
            ++cnt;
            while (cnt < 8);
            // Wait until thread 1 adds 30
          }

          ASSERT_EQ(expected[i++], id);
        }
        ASSERT_EQ(expected.size(), i);
      });

  std::vector<PartitionID> expected = { 1, 2, 5, 6, 10, 12, 14, 15, 22, 24, 30, 31 };
  size_t i = 0;
  for (const PartitionID& id : conn_set.connectivitySet()) {
    ASSERT_EQ(expected[i++], id);
  }
}

TEST(AConnectivitySet, IteratesThroughPartitionsAndSimultanouslyRemoveElements) {
  ConnectivitySet conn_set(32);
  add(conn_set, { 1, 2, 6, 10, 15, 22, 24, 31 });

  std::atomic<size_t> cnt(0);

  executeConcurrent([&] {
        while (cnt < 1);
        // Remove 6 -> since iterator already traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be not visible in the iterator
        conn_set.remove(6);
        ++cnt;

        while (cnt < 3);
        // Remove 15 -> since iterator currently traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be visible in the current iterator
        conn_set.remove(15);
        ++cnt;

        while (cnt < 5);
        // Remove 31 -> since iterator currently traverse
        // bitset containing partition ids 16 .. 24,
        // the change will be visible in the current iterator
        conn_set.remove(31);
        ++cnt;
      }, [&] {
        std::vector<PartitionID> expected = { 1, 2, 6, 10, 22, 24 };
        size_t i = 0;
        for (const PartitionID& id : conn_set.connectivitySet()) {
          if (i == 1) {
            ++cnt;
            while (cnt < 2);
            // Wait until thread 1 removes 6
          }

          if (i == 2) {
            ++cnt;
            while (cnt < 4);
            // Wait until thread 1 removes 15
          }

          if (i == 4) {
            ++cnt;
            while (cnt < 6);
            // Wait until thread 1 removes 31
          }

          ASSERT_EQ(expected[i++], id);
        }
        ASSERT_EQ(expected.size(), i);
      });

  std::vector<PartitionID> expected = { 1, 2, 10, 22, 24 };
  size_t i = 0;
  for (const PartitionID& id : conn_set.connectivitySet()) {
    ASSERT_EQ(expected[i++], id);
  }
}

TEST(AConnectivitySet, IteratesThroughPartitionsAndSimultanouslyAddAndRemoveElements) {
  ConnectivitySet conn_set(32);
  add(conn_set, { 1, 2, 6, 10, 15, 22, 24, 31 });

  std::atomic<size_t> cnt(0);

  executeConcurrent([&] {
        while (cnt < 1);
        // Add 5 -> since iterator already traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be not visible in the iterator
        conn_set.add(5);
        ++cnt;

        while (cnt < 3);
        // Remove 15 and Add 16 -> since iterator currently traverse
        // bitset containing partition ids 0 .. 8,
        // the change will be visible in the current iterator
        conn_set.remove(15);
        conn_set.add(12);
        ++cnt;

        while (cnt < 5);
        // Add 14 -> since iterator already traverse
        // bitset containing partition ids 8 .. 16,
        // the change will be not visible in the iterator
        conn_set.add(14);
        // Remove 24 -> since iterator currently traverse
        // bitset containing partition ids 8 .. 16,
        // the change will be visible in the current iterator
        conn_set.remove(24);
        ++cnt;

        while (cnt < 7);
        // Add 29 and 30 and remove 31 ->
        // since iterator currently traverse
        // bitset containing partition ids 16 .. 24,
        // the change will be visible in the current iterator
        conn_set.add(29);
        conn_set.add(30);
        conn_set.remove(31);
        ++cnt;
      }, [&] {
        std::vector<PartitionID> expected = { 1, 2, 6, 10, 12, 22, 29, 30 };
        size_t i = 0;
        for (const PartitionID& id : conn_set.connectivitySet()) {
          if (i == 1) {
            ++cnt;
            while (cnt < 2);
            // Wait until thread 1 adds 5
          }

          if (i == 2) {
            ++cnt;
            while (cnt < 4);
            // Wait until thread 1 adds 12
          }

          if (i == 4) {
            ++cnt;
            while (cnt < 6);
            // Wait until thread 1 adds 14
          }

          if (i == 5) {
            ++cnt;
            while (cnt < 8);
            // Wait until thread 1 adds 30
          }

          ASSERT_EQ(expected[i++], id);
        }
        ASSERT_EQ(expected.size(), i);
      });

  std::vector<PartitionID> expected = { 1, 2, 5, 6, 10, 12, 14, 22, 29, 30 };
  size_t i = 0;
  for (const PartitionID& id : conn_set.connectivitySet()) {
    ASSERT_EQ(expected[i++], id);
  }
}
}  // namespace ds
}  // namespace mt_kahypar
