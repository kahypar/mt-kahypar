/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <vector>
#include <atomic>
#include <tbb/scalable_allocator.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/macros.h"

namespace mt_kahypar::ds {

template<typename T>
class BufferedVector {
public:
  using vec_t = std::vector<T, tbb::scalable_allocator<T>>;

  BufferedVector(size_t max_size) :
    data(max_size, T()),
    buffers([&] { vec_t x; x.reserve(MAX_BUFFER_SIZE); return x; })
  { }

  void clear() {
    back.store(0, std::memory_order_relaxed);
    ASSERT(std::all_of(buffers.begin(), buffers.end(), [&](vec_t& x) { return x.empty(); }));
  }

  size_t size() const {
    return back.load(std::memory_order_relaxed);
  }

  size_t capacity() const {
    return data.size();
  }

  void adapt_capacity(size_t sz) {
    if (sz > data.size()) {
      data.resize(sz, T());
    }
  }

  void push_back_atomic(const T& element) {
    size_t pos = back.fetch_add(1, std::memory_order_relaxed);
    ASSERT(pos < data.size());
    data[pos] = element;
  }

  void push_back_buffered(const T& element) {
    vec_t& buffer = buffers.local();
    buffer.push_back(element);
    if (buffer.size() == MAX_BUFFER_SIZE) {
      flush_buffer(buffer);
    }
  }

  void finalize() {
    for (vec_t& buffer : buffers) {
      flush_buffer(buffer);
    }
  }

  auto begin() { return data.begin(); }
  auto end() { return data.begin() + size(); }
  T& operator[](size_t pos) { return data[pos]; }
  const T& operator[](size_t pos) const { return data[pos]; }


  struct RandomAccessRange {
    size_t actual_size;
    const vec_t& data_ref;
    const T& operator[](size_t i) const { return data_ref[i]; }
    size_t size() const { return actual_size; }
  };
  RandomAccessRange range() const { return { size(), data }; }

  const vec_t& getData() { return data; }

private:

  void flush_buffer(vec_t& buffer) {
    if (!buffer.empty()) {
      size_t pos = back.fetch_add(buffer.size(), std::memory_order_relaxed);
      ASSERT(pos + buffer.size() - 1 < data.size());
      std::copy_n(buffer.begin(), buffer.size(), data.begin() + pos);
      buffer.clear();
    }
  }

  vec_t data;
  std::atomic<size_t> back{0};
  tbb::enumerable_thread_specific<vec_t> buffers;
  static constexpr size_t MAX_BUFFER_SIZE = 1024;
};
}