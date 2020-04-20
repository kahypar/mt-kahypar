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

#include <thread>
#include <memory>

#include "tbb/parallel_for.h"
#include "tbb/scalable_allocator.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"

namespace mt_kahypar {
namespace ds {

template <typename T>
class Vector {

  class VectorIterator : public std::iterator<std::random_access_iterator_tag, T> {

    using Base = std::iterator<std::random_access_iterator_tag, T>;

    public:
      using value_type = typename Base::value_type;
      using reference = typename Base::reference;
      using pointer = typename Base::pointer;
      using difference_type = typename Base::difference_type;

      VectorIterator() : _ptr(nullptr) { }
      VectorIterator(T* ptr) : _ptr(ptr) { }
      VectorIterator(const VectorIterator& other) : _ptr(other._ptr) { }

      reference operator*() const {
        return *_ptr;
      }

      pointer operator->() const {
        return _ptr;
      }

      VectorIterator& operator++() {
        ++_ptr;
        return *this;
      }

      VectorIterator& operator--() {
        --_ptr;
        return *this;
      }

      VectorIterator operator++(int) {
        VectorIterator tmp_it(_ptr);
        ++_ptr;
        return tmp_it;
      }

      VectorIterator operator--(int) {
        VectorIterator tmp_it(_ptr);
        --_ptr;
        return tmp_it;
      }

      VectorIterator operator+(const difference_type& n) const {
        return VectorIterator(_ptr + n);
      }

      VectorIterator& operator+=(const difference_type& n) {
        _ptr += n;
        return *this;
      }

      VectorIterator operator-(const difference_type& n) const {
        return VectorIterator(_ptr - n);
      }

      VectorIterator& operator-=(const difference_type& n) {
        _ptr -= n;
        return *this;
      }

      reference operator[](const difference_type& n) const {
        return *_ptr[n];
      }

      bool operator==(const VectorIterator& other) const {
        return _ptr == other._ptr;
      }

      bool operator!=(const VectorIterator& other) const {
        return _ptr != other._ptr;
      }

      bool operator<(const VectorIterator& other) const {
        return _ptr < other._ptr;
      }

      bool operator>(const VectorIterator& other) const {
        return _ptr > other._ptr;
      }

      bool operator<=(const VectorIterator& other) const {
        return _ptr <= other._ptr;
      }

      bool operator>=(const VectorIterator& other) const {
        return _ptr >= other._ptr;
      }

      difference_type operator+(const VectorIterator& other) const {
        return ( _ptr + other._ptr );
      }

      difference_type operator-(const VectorIterator& other) const {
        return (_ptr - other._ptr);
      }

    private:
      T* _ptr;

  };

 public:

  // Type Traits
  using value_type      = T;
  using size_type       = size_t;
  using reference       = T&;
  using const_reference = const T&;
  using iterator        = VectorIterator;
  using const_iterator  = const VectorIterator;

  Vector() :
    _size(0),
    _data(nullptr),
    _underlying_data(nullptr) { }

  Vector(const size_type size,
                 const value_type init_value = value_type()) :
    _size(0),
    _data(nullptr),
    _underlying_data(nullptr) {
    resize(size, init_value);
  }

  // ####################### Access Operators #######################

  // ! Returns a reference to the element at specified location pos
  reference operator[](const size_type pos) {
    ASSERT(pos < _size);
    return _underlying_data[pos];
  }

  // ! Returns a reference to the element at specified location pos
  const_reference operator[](const size_type pos) const {
    ASSERT(pos < _size);
    return _underlying_data[pos];
  }

  reference back() {
    ASSERT(_underlying_data && _size > 0);
    return _underlying_data[_size - 1];
  }

  const_reference back() const {
    ASSERT(_underlying_data && _size > 0);
    return _underlying_data[_size - 1];
  }

  value_type* data() {
    ASSERT(_underlying_data);
    return _underlying_data;
  }

  const value_type* data() const {
    ASSERT(_underlying_data);
    return _underlying_data;
  }

  // ####################### Iterators #######################

  iterator begin() {
    ASSERT(_underlying_data);
    return iterator(_underlying_data);
  }

  const_iterator cbegin() const {
    ASSERT(_underlying_data);
    return const_iterator(_underlying_data);
  }

  iterator end() {
    ASSERT(_underlying_data);
    return iterator(_underlying_data + _size);
  }

  const_iterator cend() const {
    ASSERT(_underlying_data);
    return const_iterator(_underlying_data + _size);
  }

  // ####################### Capacity #######################

  bool empty() const {
    return _size == 0;
  }

  size_type size() const {
    return _size;
  }

  // ####################### Initialization #######################

  void resize(const size_type size,
              const value_type init_value = value_type()) {
    allocate_data(size);
    assign(size, init_value);
  }

  // ! Replaces the contents of the container
  void assign(const size_type count, const value_type value) {
    if ( _underlying_data ) {
      ASSERT(count <= _size);
      tbb::parallel_for(0UL, count, [&](const size_type i) {
        _underlying_data[i] = value;
      });
    } else {
      resize(count, value);
    }
  }

 private:
  void allocate_data(const size_type size) {
    ASSERT(!_data, "Memory already allocated");
    _data = parallel::make_unique<value_type>(size);
    _underlying_data = _data.get();
    _size = size;
  }

  size_type _size;
  parallel::tbb_unique_ptr<value_type> _data;
  value_type* _underlying_data;
};


}  // namespace ds
}  // namespace mt_kahypar