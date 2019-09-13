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
#include <type_traits>
#include <thread>

#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

namespace mt_kahypar {
namespace ds {

/**
 * Vector that allows to insert values concurrently. Internally,
 * a buffer is allocated for each cpu. A stream operation will insert
 * into the buffer of the cpu, where the calling thread is scheduled.
 * One can use copy to memcpy the data in the buffers to a global
 * vector in parallel.
 * Note, to ensure that this class is thread safe one have to make sure
 * that the calling threads are all scheduled on an unique CPU.
 */
template< typename Value >
class StreamingVector {

  static_assert( std::is_trivially_copyable<Value>::value, "Value must be trivially copyable" );

  static constexpr bool debug = false;

 public:
  StreamingVector() :
    _cpu_buffer(std::thread::hardware_concurrency()),
    _prefix_sum(std::thread::hardware_concurrency()) { }

  StreamingVector(const StreamingVector&) = delete;
  StreamingVector& operator= (const StreamingVector&) = delete;

  StreamingVector(StreamingVector&& other) = default;
  StreamingVector& operator= (StreamingVector&&) = default;

  ~StreamingVector() = default;

  template< class... Args >
  void stream( Args&&... args ) {
    int cpu_id = sched_getcpu();
    _cpu_buffer[cpu_id].emplace_back( std::forward<Args>(args)... );
  }

  std::vector<Value> copy(tbb::task_arena& arena) {
    std::vector<Value> values;
    
    size_t total_size = 0;
    for ( size_t i = 0; i < _cpu_buffer.size(); ++i) {
      _prefix_sum[i] = total_size;
      total_size += _cpu_buffer[i].size();
    }
    values.resize(total_size);

    tbb::task_group group;
    arena.execute([&] {
      for ( int cpu_id = 0; cpu_id < (int)_cpu_buffer.size(); ++cpu_id ) {
        group.run([&, cpu_id] {
          memcpy_from_cpu_buffer_to_destination(values, cpu_id, _prefix_sum[cpu_id]);
        });
      }
    });

    arena.execute([&] {
      group.wait();
    });

    return values;
  }

  const Value& value(const size_t cpu_id, const size_t idx) {
    ASSERT(cpu_id < _cpu_buffer.size());
    ASSERT(idx < _cpu_buffer[cpu_id].size());
    return _cpu_buffer[cpu_id][idx];
  }

  size_t num_buffers() const {
    return _cpu_buffer.size();
  }

  size_t size() const {
    size_t size = 0;
    for ( size_t i = 0; i < _cpu_buffer.size(); ++i ) {
      size += _cpu_buffer[i].size();
    }
    return size;
  }

  size_t size(const size_t cpu_id) const {
    ASSERT(cpu_id < _cpu_buffer.size());
    return _cpu_buffer[cpu_id].size();
  }

  size_t prefix_sum(const size_t cpu_id) const {
    ASSERT(cpu_id < _prefix_sum.size());
    return _prefix_sum[cpu_id];
  }

  void clear() {
    for ( size_t i = 0; i < _cpu_buffer.size(); ++i ) {
      std::vector<Value> tmp_value;
      _cpu_buffer[i] = std::move(tmp_value);
    }
    _prefix_sum.assign(_cpu_buffer.size(), 0);
  }

 private:

  void memcpy_from_cpu_buffer_to_destination(std::vector<Value>& destination,
                                             const int cpu_id,
                                             const size_t position) {
    DBG << "Copy buffer of cpu" << cpu_id << "of size" << _cpu_buffer[cpu_id].size()
        << "to position" << position << "in dest ( CPU =" << sched_getcpu() << " )";
    memcpy(destination.data() + position, _cpu_buffer[cpu_id].data(), 
           _cpu_buffer[cpu_id].size() * sizeof(Value)); 
  }

  std::vector<std::vector<Value>> _cpu_buffer;
  std::vector<size_t> _prefix_sum;

};

} // namespace ds
} // namespace mt_kahypar