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
#pragma once

#include <algorithm>
#include <thread>
#include <cmath>

#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Concurrent data structures to distribute key-value pairs
 * into buckets such that values with the same key reside in
 * the same bucket. Main use case for that data structure is
 * within the parallel hyperedge detection inside the multilevel
 * contractions. There, hyperedges are inserted with its footprint
 * as key into this data structure and afterwards all hyperedges
 * with the same footprint reside in the same bucket. Finally,
 * parallel hyperedges can be detected by processing each bucket
 * sequential, but several buckets in parallel.
 * If we insert a key-value pair, we compute its corresponding
 * bucket by computing key % num_buckets. To insert the key-value
 * pair, we acquire a lock on the corresponding bucket. Note,
 * key must be of type uint64_t.
 */
template <typename Value>
class ConcurrentBucketMap {

  static constexpr bool debug = false;
  static constexpr size_t BUCKET_FACTOR = 128;

  using Bucket = parallel::scalable_vector<Value>;

 public:

  ConcurrentBucketMap() :
    _num_buckets(align_to_next_power_of_two(
      BUCKET_FACTOR * std::thread::hardware_concurrency())),
    _mod_mask(_num_buckets - 1),
    _spin_locks(_num_buckets),
    _buckets(_num_buckets) { }

  ConcurrentBucketMap(const ConcurrentBucketMap&) = delete;
  ConcurrentBucketMap & operator= (const ConcurrentBucketMap &) = delete;

  ConcurrentBucketMap(ConcurrentBucketMap&& other) :
    _num_buckets(align_to_next_power_of_two(
      BUCKET_FACTOR * std::thread::hardware_concurrency())),
    _mod_mask(_num_buckets - 1),
    _spin_locks(_num_buckets),
    _buckets(std::move(other._buffer)) { }

  // ! Returns the number of buckets
  size_t numBuckets() const {
    return _num_buckets;
  }

  // ! Returns the corresponding bucket
  Bucket& getBucket(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    return _buckets[bucket];
  }

  // ! Reserves memory in each bucket such that the estimated number of insertions
  // ! can be handled without the need (with high probability) of expensive bucket resizing.
  void reserve_for_estimated_number_of_insertions(const size_t estimated_num_insertions) {
    // ! Assumption is that keys are evenly distributed among buckets (with a small buffer)
    const size_t estimated_bucket_size = std::max(
      static_cast<size_t>( 1.5 * estimated_num_insertions ) / _num_buckets, 1UL);
    tbb::parallel_for(0UL, _num_buckets, [&](const size_t i) {
      _buckets[i].reserve(estimated_bucket_size);
    });
  }

  // ! Inserts a key-value pair
  void insert(const size_t& key, Value&& value) {
    size_t bucket = key & _mod_mask;
    ASSERT(bucket < _num_buckets);
    _spin_locks[bucket].lock();
    _buckets[bucket].emplace_back( std::move(value) );
    _spin_locks[bucket].unlock();
  }

  // ! Frees the memory of all buckets
  void free() {
    parallel::parallel_free(_buckets);
  }

  // ! Frees the memory of the corresponding bucket
  void free(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    parallel::free(_buckets[bucket]);
  }

  // ! Clears the corresponding bucket
  void clear(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    _buckets[bucket].clear();
  }

 private:
  size_t align_to_next_power_of_two(const size_t size) const {
    return std::pow(2.0, std::ceil(std::log2(static_cast<double>(size))));
  }

  const size_t _num_buckets;
  const size_t _mod_mask;
  std::vector<SpinLock> _spin_locks;
  parallel::scalable_vector<Bucket> _buckets;
};
}  // namespace ds
}  // namespace mt_kahypar
