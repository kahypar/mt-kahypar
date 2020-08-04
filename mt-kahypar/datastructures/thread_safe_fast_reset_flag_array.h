/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

// based on http://upcoder.com/9/fast-resettable-flag-vector/

namespace mt_kahypar {
namespace ds {
template <typename Type = std::uint16_t>
class ThreadSafeFastResetFlagArray {

  using UnderlyingType = CAtomic<Type>;

  public:
  explicit ThreadSafeFastResetFlagArray(const size_t size) :
    _v(std::make_unique<UnderlyingType[]>(size)),
    _threshold(1),
    _size(size) {
    init();
  }

  ThreadSafeFastResetFlagArray() :
    _v(nullptr),
    _threshold(1),
    _size(0) { }

  ThreadSafeFastResetFlagArray(const ThreadSafeFastResetFlagArray&) = delete;
  ThreadSafeFastResetFlagArray& operator= (const ThreadSafeFastResetFlagArray&) = delete;

  ThreadSafeFastResetFlagArray(ThreadSafeFastResetFlagArray&&) = default;
  ThreadSafeFastResetFlagArray& operator= (ThreadSafeFastResetFlagArray&&) = default;

  ~ThreadSafeFastResetFlagArray() = default;

  void swap(ThreadSafeFastResetFlagArray& other) {
    using std::swap;
    swap(_v, other._v);
    swap(_threshold, other._threshold);
  }

  bool operator[] (const size_t i) const {
    return isSet(i);
  }

  // ! Changes value of entry i from false to true and returns true, if the value
  // ! hold on position i was false and was successfully set to true
  bool compare_and_set_to_true(const size_t i) {
    Type expected = _v[i].load(std::memory_order_relaxed);
    Type desired = _threshold;
    if ( expected != _threshold && _v[i].compare_exchange_strong(expected, desired) ) {
      // Value was successfully set from false to true
      return true;
    } else {
      // Either expected == _threshold or compare_exchange_strong failed, which means that
      // an other thread set _v[i] to true before.
      return false;
    }
  }

  void set(const size_t i, const bool value) {
    _v[i].store(value ? _threshold : 0, std::memory_order_relaxed);
  }

  void reset() {
    if (_threshold == std::numeric_limits<UnderlyingType>::max()) {
      init();
      _threshold = 0;
    }
    ++_threshold;
  }

  void setSize(const size_t size, const bool initialiser = false) {
    ASSERT(_v == nullptr, "Error");
    _v = std::make_unique<UnderlyingType[]>(size);
    _size = size;
    init(initialiser);
  }

 private:
  bool isSet(size_t i) const {
    return _v[i].load(std::memory_order_relaxed) == _threshold;
  }

  void init(const bool initialiser = false) {
    const Type init_value = initialiser ? _threshold : 0;
    for ( size_t i = 0; i < _size; ++i ) {
      _v[i].store(init_value, std::memory_order_relaxed);
    }
  }

  std::unique_ptr<UnderlyingType[]> _v;
  Type _threshold;
  size_t _size;
};

template <typename UnderlyingType>
void swap(ThreadSafeFastResetFlagArray<UnderlyingType>& a,
          ThreadSafeFastResetFlagArray<UnderlyingType>& b) {
  a.swap(b);
}
}  // namespace ds
}  // namespace mt_kahypar
