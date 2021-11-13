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
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/datastructures/shared_mapping.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using SharedMap = ds::SharedMapping<size_t, size_t, size_t>;

template <class F1, class F2>
void executeConcurrent(const F1& f1, const F2& f2) {
  std::atomic<int> cnt(0);
  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < 2) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < 2) { }
    f2();
  });
}

template <class F1, class F2, class F3>
void executeConcurrent(const F1& f1, const F2& f2, const F3& f3) {
  std::atomic<int> cnt(0);
  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < 3) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < 3) { }
    f2();
  }, [&] {
    cnt++;
    while (cnt < 3) { }
    f3();
  });
}

TEST(ASharedMapping, InsertsElementsWithDistinctKeys) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,1,2);
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,1));
  ASSERT_EQ(2, *sm.get_if_contained(1,1));
}

TEST(ASharedMapping, InsertsElementsWithDistinctKeysConcurrently) {
  SharedMap sm(10);
  executeConcurrent([&] {
    sm.add(0,0,1);
  }, [&] {
    sm.add(1,1,2);
  });
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,1));
  ASSERT_EQ(2, *sm.get_if_contained(1,1));
}

TEST(ASharedMapping, InsertsElementsWithSameKeys1) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
}

TEST(ASharedMapping, InsertsElementsWithSameKeys2) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_TRUE(sm.contains(2,0));
  ASSERT_EQ(3, *sm.get_if_contained(2,0));
}

TEST(ASharedMapping, InsertsElementsWithSameKeysConcurrently1) {
  SharedMap sm(10);
  executeConcurrent([&] {
    sm.add(0,0,1);
  }, [&] {
    sm.add(1,0,2);
  });
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
}

TEST(ASharedMapping, InsertsElementsWithSameKeysConcurrently2) {
  SharedMap sm(10);
  executeConcurrent([&] {
    sm.add(0,0,1);
  }, [&] {
    sm.add(1,0,2);
  }, [&] {
    sm.add(2,0,3);
  });
  ASSERT_TRUE(sm.contains(0,0));
  ASSERT_EQ(1, *sm.get_if_contained(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_TRUE(sm.contains(2,0));
  ASSERT_EQ(3, *sm.get_if_contained(2,0));
}

TEST(ASharedMapping, RemoveElements1) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.remove(0, 0);
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
}

TEST(ASharedMapping, RemoveElements2) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  sm.remove(0, 0);
  sm.remove(2, 0);
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_FALSE(sm.contains(2,0));
}

TEST(ASharedMapping, RemoveElementsConcurrently2) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  executeConcurrent([&] {
    sm.remove(0, 0);
  }, [&] {
    sm.remove(2, 0);
  });
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_FALSE(sm.contains(2,0));
}

TEST(ASharedMapping, RemoveElementsAndAddAgain1) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  sm.remove(0, 0);
  sm.remove(2, 0);
  sm.add(4,0,4);
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_FALSE(sm.contains(2,0));
  ASSERT_TRUE(sm.contains(4,0));
  ASSERT_EQ(4, *sm.get_if_contained(4,0));
}

TEST(ASharedMapping, RemoveElementsAndAddAgain2) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  sm.remove(0, 0);
  sm.remove(2, 0);
  sm.add(4,0,4);
  sm.add(5,0,5);
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_FALSE(sm.contains(2,0));
  ASSERT_TRUE(sm.contains(4,0));
  ASSERT_EQ(4, *sm.get_if_contained(4,0));
  ASSERT_TRUE(sm.contains(5,0));
  ASSERT_EQ(5, *sm.get_if_contained(5,0));
}

TEST(ASharedMapping, RemoveElementsAndAddAgainConcurrently) {
  SharedMap sm(10);
  sm.add(0,0,1);
  sm.add(1,0,2);
  sm.add(2,0,3);
  sm.remove(0, 0);
  sm.remove(2, 0);
  sm.add(4,0,4);
  executeConcurrent([&] {
    sm.add(5,0,5);
  }, [&] {
    sm.add(6,0,6);
  });
  ASSERT_FALSE(sm.contains(0,0));
  ASSERT_TRUE(sm.contains(1,0));
  ASSERT_EQ(2, *sm.get_if_contained(1,0));
  ASSERT_FALSE(sm.contains(2,0));
  ASSERT_TRUE(sm.contains(4,0));
  ASSERT_EQ(4, *sm.get_if_contained(4,0));
  ASSERT_TRUE(sm.contains(5,0));
  ASSERT_EQ(5, *sm.get_if_contained(5,0));
  ASSERT_TRUE(sm.contains(6,0));
  ASSERT_EQ(6, *sm.get_if_contained(6,0));
}

struct Element {
  size_t key;
  bool contained;
};

TEST(ASharedMapping, StressTest) {
  const size_t N = 1000;
  SharedMap sm(N);

  tbb::parallel_for(0U, std::thread::hardware_concurrency(), [&](const uint32_t i) {
    const size_t id = i;
    vec<Element> elems;
    for ( size_t j = 0; j < N; ++j ) {
      elems.push_back(Element { j, utils::Randomize::instance().flipCoin(i) });
    }
    std::random_shuffle(elems.begin(), elems.end());
    for ( const Element& elem : elems ) {
      if ( elem.contained ) {
        sm.add(id, elem.key, id);
      }
    }

    for ( const Element& elem : elems ) {
      if ( elem.contained ) {
        ASSERT_TRUE(sm.contains(id, elem.key));
        ASSERT_EQ(id, *sm.get_if_contained(id,elem.key));
      } else {
        ASSERT_FALSE(sm.contains(id, elem.key));
      }
    }
  });
}

}  // namespace ds
}  // namespace mt_kahypar
