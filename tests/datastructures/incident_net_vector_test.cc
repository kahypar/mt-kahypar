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
#include <mt-kahypar/macros.h>

#include "gmock/gmock.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/datastructures/incident_net_vector.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

template<typename T>
void compare_equal(IncidentNetVector<T>& actual,
                   const parallel::scalable_vector<T>& expected) {
  std::sort(actual.begin(), actual.end());
  ASSERT_EQ(actual.size(), expected.size());
  for ( size_t i = 0; i < actual.size(); ++i ) {
    ASSERT_EQ(actual[i], expected[i]);
  }
}

TEST(AIncidentVector, UpdatesActiveIterators1) {
  IncidentNetVector<int> vec;
  auto a = vec.cbegin();
  ASSERT_EQ(1, vec.active_iterators());
  auto b = vec.cbegin();
  ASSERT_EQ(2, vec.active_iterators());
}

TEST(AIncidentVector, UpdatesActiveIterators2) {
  IncidentNetVector<int> vec;
  {
    auto a = vec.cbegin();
    ASSERT_EQ(1, vec.active_iterators());
    {
      auto b = vec.begin();
      ASSERT_EQ(2, vec.active_iterators());
    }
    ASSERT_EQ(1, vec.active_iterators());
  }
  ASSERT_EQ(0, vec.active_iterators());
}

TEST(AIncidentVector, UpdatesActiveIterators3) {
  IncidentNetVector<int> vec;
  {
    auto a = vec.cbegin();
    ASSERT_EQ(1, vec.active_iterators());
    {
      auto b = vec.cbegin();
      ASSERT_EQ(2, vec.active_iterators());
    }
    ASSERT_EQ(1, vec.active_iterators());

    {
      auto b = vec.end();
      ASSERT_EQ(2, vec.active_iterators());
    }
    ASSERT_EQ(1, vec.active_iterators());
  }
  ASSERT_EQ(0, vec.active_iterators());
}

TEST(AIncidentVector, UpdatesActiveIteratorsInParallel1) {
  IncidentNetVector<int> vec;
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    auto a = vec.cbegin();
    ++cnt;
    while ( cnt < 2 ) { }
    ASSERT_EQ(2, vec.active_iterators());
    ++cnt;
    while ( cnt < 4 ) { }
  }, [&] {
    auto b = vec.cend();
    ++cnt;
    while ( cnt < 2 ) { }
    ASSERT_EQ(2, vec.active_iterators());
    ++cnt;
    while ( cnt < 4 ) { }
  });
  ASSERT_EQ(0, vec.active_iterators());
}

TEST(AIncidentVector, BulkInsert1) {
  IncidentNetVector<int> vec;
  vec.bulk_insert({2,1,3,4,5});
  compare_equal(vec, {1,2,3,4,5});
}

TEST(AIncidentVector, BulkInsert2) {
  IncidentNetVector<int> vec;
  vec.bulk_insert({2,1,3,4,5});
  vec.bulk_insert({10,9,7,8,6});
  compare_equal(vec, {1,2,3,4,5,6,7,8,9,10});
  ASSERT_EQ(0, vec.active_iterators());
}

TEST(AIncidentVector, BulkInsertInParallel1) {
  IncidentNetVector<int> vec;
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while ( cnt < 2 ) { }
    vec.bulk_insert({2,1,3,4,5});
  }, [&] {
    ++cnt;
    while ( cnt < 2 ) { }
    vec.bulk_insert({10,9,7,8,6});
  });
  compare_equal(vec, {1,2,3,4,5,6,7,8,9,10});
  ASSERT_EQ(0, vec.active_iterators());
}

TEST(AIncidentVector, BulkInsertInParallel2) {
  IncidentNetVector<int> vec;
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while ( cnt < 2 ) { }
    vec.bulk_insert({2,1,3,4,5});
    vec.bulk_insert({10,9,7,8,6});
  }, [&] {
    ++cnt;
    while ( cnt < 2 ) { }
    vec.bulk_insert({10,9,7,8,6});
    vec.bulk_insert({2,1,3,4,5});
  });
  compare_equal(vec, {1,1,2,2,3,3,4,4,5,5,
                      6,6,7,7,8,8,9,9,10,10});
  ASSERT_EQ(0, vec.active_iterators());
}

#ifndef KAHYPAR_TRAVIS_BUILD
TEST(AIncidentVector, BulkInsertWhileIteratorsAreActive) {
  IncidentNetVector<int> vec;
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while ( cnt < 2 ) { }
    {
      auto a = vec.cbegin();
      ++cnt;
      compare_equal(vec, {});
    }
    while ( cnt < 5 ) { }
    vec.bulk_insert({1,2,3,4,5});
    compare_equal(vec, {1,1,2,2,3,3,4,4,5,5});
    ++cnt;
  }, [&] {
    ++cnt;
    while ( cnt < 3 ) { }
    vec.bulk_insert({1,2,3,4,5});
    ++cnt;
    {
      auto a = vec.cbegin();
      ++cnt;
      compare_equal(vec, {1,2,3,4,5});
    }
  });
  ASSERT_EQ(0, vec.active_iterators());
}
#endif

}  // namespace ds
}  // namespace mt_kahypar
