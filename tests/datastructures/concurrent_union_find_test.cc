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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/concurrent_union_find.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

class HypergraphMock {

 public:
  explicit HypergraphMock(const std::vector<HypernodeWeight>& weights) :
    _weights(weights) { }

  HypernodeID initialNumNodes() const {
    return _weights.size();
  }

  bool nodeIsEnabled(const HypernodeID) const {
    return true;
  }

  bool isHighDegreeVertex(const HypernodeID) const {
    return false;
  }

  HypernodeID globalNodeID(const HypernodeID hn) const {
    return hn;
  }

  HypernodeWeight nodeWeight(const HypernodeID hn) const {
    return _weights[hn];
  }

 private:
  std::vector<HypernodeWeight> _weights;
};

using UnionFind = ConcurrentUnionFind<HypergraphMock>;


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

TEST(AConcurrentUnionFind, LinksTwoVertices1) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  uf.link(0, 1);
  ASSERT_EQ(7, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(0, 1));
  ASSERT_EQ(1, uf.find(0));
  ASSERT_EQ(1, uf.find(1));
  ASSERT_EQ(3, uf.weight(0));
  ASSERT_EQ(3, uf.weight(1));
}

TEST(AConcurrentUnionFind, LinksTwoVertices2) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  uf.link(5, 7);
  ASSERT_EQ(7, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(5, 7));
  ASSERT_EQ(5, uf.find(5));
  ASSERT_EQ(5, uf.find(7));
  ASSERT_EQ(4, uf.weight(5));
  ASSERT_EQ(4, uf.weight(7));
}

TEST(AConcurrentUnionFind, LinksTwoVertices3) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  uf.link(4, 6);
  ASSERT_EQ(7, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(4, 6));
  ASSERT_EQ(4, uf.find(4));
  ASSERT_EQ(4, uf.find(6));
  ASSERT_EQ(2, uf.weight(4));
  ASSERT_EQ(2, uf.weight(6));
}

TEST(AConcurrentUnionFind, LinksSeveralVertices1) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  uf.link(1, 2);
  uf.link(4, 6);
  ASSERT_EQ(6, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(1, 2));
  ASSERT_TRUE(uf.isSameSet(4, 6));
  ASSERT_FALSE(uf.isSameSet(2, 4));
  ASSERT_EQ(1, uf.find(2));
  ASSERT_EQ(3, uf.weight(2));
  ASSERT_EQ(4, uf.find(6));
  ASSERT_EQ(2, uf.weight(6));
}

TEST(AConcurrentUnionFind, LinksSeveralVertices2) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  uf.link(3, 5);
  uf.link(5, 7);
  ASSERT_EQ(6, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(3, 5));
  ASSERT_TRUE(uf.isSameSet(5, 7));
  ASSERT_EQ(3, uf.find(5));
  ASSERT_EQ(6, uf.weight(5));
  ASSERT_EQ(3, uf.find(7));
  ASSERT_EQ(6, uf.weight(7));
}

TEST(AConcurrentUnionFind, LinksSeveralVerticesConcurrent1) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  executeConcurrent([&] {
        uf.link(1, 2);
      }, [&] {
        uf.link(4, 6);
      });

  ASSERT_EQ(6, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(1, 2));
  ASSERT_TRUE(uf.isSameSet(4, 6));
  ASSERT_FALSE(uf.isSameSet(2, 4));
  ASSERT_EQ(1, uf.find(2));
  ASSERT_EQ(3, uf.weight(2));
  ASSERT_EQ(4, uf.find(6));
  ASSERT_EQ(2, uf.weight(6));
}

TEST(AConcurrentUnionFind, LinksSeveralVerticesConcurrent2) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  executeConcurrent([&] {
        uf.link(3, 5);
      }, [&] {
        uf.link(5, 7);
      });

  ASSERT_EQ(6, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(3, 5));
  ASSERT_TRUE(uf.isSameSet(5, 7));

  const HypernodeID root = uf.find(3);
  ASSERT_EQ(root, uf.find(3));
  ASSERT_EQ(6, uf.weight(3));
  ASSERT_EQ(root, uf.find(5));
  ASSERT_EQ(6, uf.weight(5));
  ASSERT_EQ(root, uf.find(7));
  ASSERT_EQ(6, uf.weight(7));
}

TEST(AConcurrentUnionFind, LinksSeveralVerticesConcurrent3) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  executeConcurrent([&] {
        uf.link(1, 2);
        uf.link(3, 5);
      }, [&] {
        uf.link(1, 4);
        uf.link(5, 7);
      });

  ASSERT_EQ(4, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(3, 5));
  ASSERT_TRUE(uf.isSameSet(5, 7));
  ASSERT_TRUE(uf.isSameSet(1, 2));
  ASSERT_TRUE(uf.isSameSet(2, 4));
  ASSERT_FALSE(uf.isSameSet(4, 7));

  HypernodeID root = uf.find(3);
  ASSERT_EQ(root, uf.find(3));
  ASSERT_EQ(6, uf.weight(3));
  ASSERT_EQ(root, uf.find(5));
  ASSERT_EQ(6, uf.weight(5));
  ASSERT_EQ(root, uf.find(7));
  ASSERT_EQ(6, uf.weight(7));

  root = uf.find(1);
  ASSERT_EQ(root, uf.find(1));
  ASSERT_EQ(4, uf.weight(1));
  ASSERT_EQ(root, uf.find(2));
  ASSERT_EQ(4, uf.weight(2));
  ASSERT_EQ(root, uf.find(4));
  ASSERT_EQ(4, uf.weight(4));
}

TEST(AConcurrentUnionFind, LinksSeveralVerticesConcurrent4) {
  HypergraphMock mock({1, 2, 1, 2, 1, 2, 1, 2});
  UnionFind uf(mock);

  executeConcurrent([&] {
        uf.link(1, 2);
        uf.link(3, 5);
        uf.link(0, 6);
      }, [&] {
        uf.link(1, 4);
        uf.link(5, 7);
        uf.link(1, 5);
      });

  ASSERT_EQ(2, uf.numDistinctSets());
  ASSERT_TRUE(uf.isSameSet(0, 6));
  ASSERT_TRUE(uf.isSameSet(1, 2));
  ASSERT_TRUE(uf.isSameSet(2, 3));
  ASSERT_TRUE(uf.isSameSet(4, 5));
  ASSERT_TRUE(uf.isSameSet(5, 7));
  ASSERT_FALSE(uf.isSameSet(6, 7));

  HypernodeID root = uf.find(1);
  ASSERT_EQ(root, uf.find(1));
  ASSERT_EQ(10, uf.weight(1));
  ASSERT_EQ(root, uf.find(2));
  ASSERT_EQ(10, uf.weight(2));
  ASSERT_EQ(root, uf.find(3));
  ASSERT_EQ(10, uf.weight(3));
  ASSERT_EQ(root, uf.find(4));
  ASSERT_EQ(10, uf.weight(4));
  ASSERT_EQ(root, uf.find(7));
  ASSERT_EQ(10, uf.weight(7));

  root = uf.find(0);
  ASSERT_EQ(root, uf.find(0));
  ASSERT_EQ(2, uf.weight(0));
  ASSERT_EQ(root, uf.find(6));
  ASSERT_EQ(2, uf.weight(6));
}

TEST(AConcurrentUnionFind, SmokeTest) {
  const size_t N = 100000;
  const size_t Q = 1000;

  std::vector<HypernodeWeight> weights;
  for ( size_t i = 0; i < N; ++i ) {
    weights.emplace_back(utils::Randomize::instance().getRandomInt(1, 100, sched_getcpu()));
  }

  HypergraphMock mock(weights);
  UnionFind parallel_uf(mock);
  UnionFind sequential_uf(mock);

  std::vector<std::pair<HypernodeID, HypernodeID>> queries;
  for ( size_t i = 0; i < Q; ++i ) {
    queries.emplace_back(std::make_pair(
      utils::Randomize::instance().getRandomInt(0, N - 1, sched_getcpu()),
      utils::Randomize::instance().getRandomInt(0, N - 1, sched_getcpu())));
    sequential_uf.link(queries.back().first, queries.back().second);
  }

  tbb::parallel_for(0UL, Q, [&](const size_t q) {
    parallel_uf.link(queries[q].first, queries[q].second);
  });

  ASSERT_EQ(sequential_uf.numDistinctSets(), parallel_uf.numDistinctSets());
  for ( size_t i = 0; i < Q; ++i ) {
    HypernodeID u = queries[i].first;
    HypernodeID v = queries[i].second;

    ASSERT_TRUE(parallel_uf.isSameSet(u, v));
    ASSERT_EQ(sequential_uf.weight(u), parallel_uf.weight(u));
    ASSERT_EQ(sequential_uf.weight(v), parallel_uf.weight(v));
  }
}

}  // namespace ds
}  // namespace mt_kahypar
