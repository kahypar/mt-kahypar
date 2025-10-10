/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <type_traits>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/weight/hypernode_weight_base.h"

namespace mt_kahypar {
namespace weight {

class HypernodeWeightArray {

  template<bool is_const>
  class WeightArrayIterator {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::conditional_t<is_const, HNWeightConstRef, HNWeightRef>;
    using reference = value_type;
    using pointer = const value_type*;
    using difference_type = std::ptrdiff_t;

    explicit WeightArrayIterator(HNWeightScalar* ptr, Dimension dimension) :
      _ptr(ptr), _buffer(ptr, dimension) { }

    WeightArrayIterator(const WeightArrayIterator& other) :
      _ptr(other._ptr), _buffer(weight::copy(other._buffer)) { }

    WeightArrayIterator& operator=(const WeightArrayIterator& other) {
        _ptr = other._ptr;
        weight::replace(_buffer, weight::copy(other._buffer));
        return *this;
      }

    reference operator*() const {
      ASSERT(_ptr == _buffer.get_raw_data());
      return weight::copy(_buffer);
    }

    pointer operator->() const {
      ASSERT(_ptr == _buffer.get_raw_data());
      return &_buffer;
    }

    WeightArrayIterator& operator++() {
      return *this += 1;
    }

    WeightArrayIterator& operator--() {
      return *this -= 1;
    }

    WeightArrayIterator operator++(int) {
      WeightArrayIterator tmp_it(*this);
      *this += 1;
      return tmp_it;
    }

    WeightArrayIterator operator--(int) {
      WeightArrayIterator tmp_it(*this);
      *this -= 1;
      return tmp_it;
    }

    WeightArrayIterator operator+(const difference_type& n) const {
      WeightArrayIterator tmp_it(*this);
      tmp_it += n;
      return tmp_it;
    }

    WeightArrayIterator& operator+=(const difference_type& n) {
      _ptr += n * _buffer.dimension();
      _buffer.set_raw_data(_ptr, _buffer.dimension());
      return *this;
    }

    WeightArrayIterator operator-(const difference_type& n) const {
      WeightArrayIterator tmp_it(*this);
      tmp_it -= n;
      return tmp_it;
    }

    WeightArrayIterator& operator-=(const difference_type& n) {
      _ptr -= n * _buffer.dimension();
      _buffer.set_raw_data(_ptr, _buffer.dimension());
      return *this;
    }

    reference operator[](const difference_type& n) const {
      ASSERT(_ptr == _buffer.get_raw_data());
      return reference(_ptr + n * _buffer.dimension(), _buffer.dimension());
    }

    bool operator==(const WeightArrayIterator& other) const {
      return _ptr == other._ptr;
    }

    bool operator!=(const WeightArrayIterator& other) const {
      return _ptr != other._ptr;
    }

    bool operator<(const WeightArrayIterator& other) const {
      return _ptr < other._ptr;
    }

    bool operator>(const WeightArrayIterator& other) const {
      return _ptr > other._ptr;
    }

    bool operator<=(const WeightArrayIterator& other) const {
      return _ptr <= other._ptr;
    }

    bool operator>=(const WeightArrayIterator& other) const {
      return _ptr >= other._ptr;
    }

    difference_type operator-(const WeightArrayIterator& other) const {
      return (_ptr - other._ptr) / _buffer.dimension();
    }

    friend WeightArrayIterator operator+(const difference_type& n, const WeightArrayIterator& it) {
      return it + n;
    }

   private:
    HNWeightScalar* _ptr;
    reference _buffer;
  };

 public:
  // Type Traits
  // using value_type      = ??;
  using size_type       = size_t;
  using reference       = HNWeightRef;
  using const_reference = HNWeightConstRef;
  using iterator        = WeightArrayIterator<false>;
  using const_iterator  = WeightArrayIterator<true>;

  HypernodeWeightArray() :
    _data(),
    _dimension(0) { }

  // note: uses parallel assignment
  HypernodeWeightArray(const size_type size,
                       const Dimension dimension,
                       const weight::HNWeightScalar init_value,
                       bool parallel = true) :
    _data(),
    _dimension(0) {
      resize(size, dimension, init_value, parallel);
    }

  HypernodeWeightArray(const HypernodeWeightArray&) = delete;
  HypernodeWeightArray & operator= (const HypernodeWeightArray&) = delete;

  HypernodeWeightArray(HypernodeWeightArray&& other) :
    _data(std::move(other._data)),
    _dimension(other._dimension) {
    other._dimension = 0;
  }

  HypernodeWeightArray& operator=(HypernodeWeightArray&& other) {
    _data = std::move(other._data);
    _dimension = other._dimension;
    other._dimension = 0;
    return *this;
  }

  // ####################### Access Operators #######################

  // ! Returns a reference to the element at specified location pos
  reference operator[](const size_type pos) {
    ASSERT(_dimension > 0);
    ASSERT(pos < size() && pos * _dimension < _data.size());
    return weight::HNWeightRef(_data.data() + (pos * _dimension), _dimension);
  }

  // ! Returns a reference to the element at specified location pos
  const_reference operator[](const size_type pos) const {
    ASSERT(_dimension > 0);
    ASSERT(pos < size() && pos * _dimension < _data.size());
    return weight::HNWeightConstRef(_data.data() + (pos * _dimension), _dimension);
  }

  weight::HNWeightScalar* data() {
    return _data.data();
  }

  const weight::HNWeightScalar* data() const {
    return _data.data();
  }

  ds::Array<weight::HNWeightScalar>& underlying_array() {
    return _data;
  }

  const ds::Array<weight::HNWeightScalar>& underlying_array() const {
    return _data;
  }

  // ####################### Iterators #######################

  iterator begin() {
    return iterator(data(), _dimension);
  }

  const_iterator begin() const {
    return cbegin();
  }

  const_iterator cbegin() const {
    return const_iterator(const_cast<weight::HNWeightScalar*>(data()), _dimension);
  }

  iterator end() {
    ASSERT(_data.size() % _dimension == 0);
    return iterator(data() + _data.size(), _dimension);
  }

  const_iterator end() const {
    return cend();
  }

  const_iterator cend() const {
    ASSERT(_data.size() % _dimension == 0);
    return const_iterator(const_cast<weight::HNWeightScalar*>(data()) + _data.size(), _dimension);
  }

  // ####################### Capacity #######################

  bool empty() const {
    return _data.empty();
  }

  size_type size() const {
    ASSERT(_dimension > 0);
    // TODO: it would probably be better to store the size, but is
    // not really compatible with array
    ASSERT(_data.size() % _dimension == 0);
    return _data.size() / _dimension;
  }

  Dimension dimension() const {
    return _dimension;
  }

  // ####################### Initialization #######################

  void resize(const size_type size,
              const Dimension dimension,
              const weight::HNWeightScalar init_value = 0,
              const bool assign_parallel = true) {
    _data.resize(size * dimension, init_value, assign_parallel);
    _dimension = dimension;
  }

  void replaceWith(const size_type size,
                   const Dimension dimension,
                   const weight::HNWeightScalar init_value,
                   const bool assign_parallel = true) {
    HypernodeWeightArray new_array;
    new_array.resize(size, dimension, init_value, assign_parallel);
    *this = std::move(new_array);
  }

  bool changeDimension(const Dimension new_dimension) {
    bool success = _data.size() % new_dimension == 0;
    if (success) {
      _dimension = new_dimension;
    }
    return success;
  }

  void assign(const size_type count,
              const weight::HNWeightScalar value = 0,
              const bool assign_parallel = true) {
    ASSERT(_dimension > 0 && count <= size());
    _data.assign(count * _dimension, value, assign_parallel);
  }

  // ####################### COPY #######################

  // sequential!
  HypernodeWeightArray copy() const {
    if (empty()) return HypernodeWeightArray();

    HypernodeWeightArray result;
    result._data.resize(_data.size(), 0, false);
    std::memcpy(result._data.data(), _data.data(), sizeof(weight::HNWeightScalar) * _data.size());
    result._dimension = _dimension;
    return result;
  }

 private:
  ds::Array<weight::HNWeightScalar> _data;
  Dimension _dimension;
};


class CopyableHypernodeWeightArray: public HypernodeWeightArray {
 public:
  CopyableHypernodeWeightArray() : HypernodeWeightArray() { }

  // note: uses parallel assignment
  CopyableHypernodeWeightArray(const size_type size,
                               const Dimension dimension,
                               const weight::HNWeightScalar init_value = 0) :
    HypernodeWeightArray(size, dimension, init_value) { }

  CopyableHypernodeWeightArray(const CopyableHypernodeWeightArray& other) : HypernodeWeightArray(other.copy()) { }

  CopyableHypernodeWeightArray& operator=(const CopyableHypernodeWeightArray& other) {
    static_cast<HypernodeWeightArray&>(*this) = other.copy();
    return *this;
  }

  CopyableHypernodeWeightArray(CopyableHypernodeWeightArray&& other) = default;
  CopyableHypernodeWeightArray& operator=(CopyableHypernodeWeightArray&& other) = default;
};

}  // namespace weight
}  // namespace mt_kahypar
