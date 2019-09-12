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

#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

namespace kahypar {
namespace ds {

/**
 * Allows to insert key-value pairs concurrently. Internally, 
 * the hash of a key is used to map a key-value pair to a bucket.
 * Once the bucket is determined a lock is acquired and the key-value pair
 * is inserted.
 * By using the copy function one can map the key-value pairs to a unique
 * bucket in parallel (see internal implementation).
 * 
 * For example:
 * Bucket 0: (1, 2), (1, 3), (2, 1), (1, 4)
 * Bucket 1: (3, 5), (4, 1), (3, 2), (4, 2)
 * 
 * After copy:
 * Bucket 0: 
 * Bucket 1: 2 -> 3 -> 4
 * Bucket 2: 1
 * Bucket 3: 5 -> 2
 * Bucket 4: 1 -> 2
 */
template< typename Key, typename Value >
class StreamingMap {

  static constexpr bool debug = false;
  static constexpr size_t BUCKET_FACTOR = 16;

  using KeyValuePair = std::pair<Key, Value>;

 public:
  StreamingMap() :
    _size(BUCKET_FACTOR * std::thread::hardware_concurrency()),
    _mutex(_size),
    _buffer(_size) { }

  StreamingMap(const StreamingMap&) = delete;
  StreamingMap& operator= (const StreamingMap&) = delete;

  StreamingMap(StreamingMap&& other) :
    _size(BUCKET_FACTOR * std::thread::hardware_concurrency()),
    _mutex(_size),
    _buffer(std::move(other._buffer)) { }

  StreamingMap& operator= (StreamingMap&&) = default;

  ~StreamingMap() = default;

  void stream(const Key& key, const Value& value) {
    size_t bucket = std::hash<Key>()(key) % _size;
    ASSERT(bucket < _size);
    std::lock_guard<std::mutex> lock(_mutex[bucket]);
    _buffer[bucket].emplace_back(key, value);
  }

  template < class F >
  void copy(tbb::task_arena& arena,
            std::vector<std::vector<Value>>& destination,
            F&& key_extractor) {
    tbb::task_group group;
    arena.execute([&] {
      for ( size_t bucket = 0; bucket < _buffer.size(); ++bucket ) {
        group.run([&, bucket] {
          ASSERT(bucket < _buffer.size());
          for ( const KeyValuePair& p : _buffer[bucket] ) {
            const Key& key = key_extractor(p.first);
            ASSERT(key < destination.size());
            destination[key].emplace_back(p.second);
          }
        });
      }
    });

    arena.execute([&] {
      group.wait();
    });
  }

  void clear() {
    for ( size_t i = 0; i < _buffer.size(); ++i ) {
      std::vector<KeyValuePair> tmp_buffer;
      _buffer[i] = std::move(tmp_buffer);
    }
  }

 private:
  const size_t _size;
  std::vector<std::mutex> _mutex;
  std::vector<std::vector<KeyValuePair>> _buffer;
};

} // namespace ds
} // namespace kahypar