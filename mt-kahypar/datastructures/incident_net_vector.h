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
#include <type_traits>
#include <shared_mutex>

#include "tbb/atomic.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

/**
 * Special vector that allows to perform bulk inserts and concurrently iterate over
 * the vector without invalidating pointers. Each iterator increments a reference count
 * and bulk insert is only allowed if corresponding thread hold a unique lock and that
 * reference counter is equal to zero.
 *
 * Note, that this data structure supports all operation of a std::vector, but they are not
 * thread-safe and should be used carefully.
 */
template <typename T>
class IncidentNetVector : public parallel::scalable_vector<T> {

  static constexpr bool debug = false;
  using Base = parallel::scalable_vector<T>;

  template<typename V,
           typename BaseIterator,
           typename IncidentNetVectorType>
  class IncidentNetIterator : public std::iterator<std::random_access_iterator_tag, V> {

    using Base = std::iterator<std::random_access_iterator_tag, V>;

    public:
      using value_type = typename Base::value_type;
      using reference = typename Base::reference;
      using pointer = typename Base::pointer;
      using difference_type = typename Base::difference_type;

      IncidentNetIterator(BaseIterator it, IncidentNetVectorType* vec) :
        _it(it),
        _vec(vec) {
        ++_vec->_ref_count;
      }

      IncidentNetIterator(const IncidentNetIterator& other) :
        _it(other._it),
        _vec(other._vec) {
        ++_vec->_ref_count;
      }

      IncidentNetIterator & operator= (const IncidentNetIterator & other) {
        _it = other._it;
        _vec = other._vec;
        return *this;
      }

      IncidentNetIterator(IncidentNetIterator&& other) :
        _it(std::move(other._it)),
        _vec(other._vec) {
        other._vec = nullptr;
      }

      IncidentNetIterator & operator=(IncidentNetIterator&& other) {
        _it = std::move(other._it);
        _vec = other._vec;
        other._vec = nullptr;
        return *this;
      }

      ~IncidentNetIterator() {
        if ( _vec ) {
          --_vec->_ref_count;
        }
      }

      reference operator*() const {
        return *_it;
      }

      pointer operator->() const {
        return _it;
      }

      IncidentNetIterator& operator++() {
        ++_it;
        return *this;
      }

      IncidentNetIterator& operator--() {
        --_it;
        return *this;
      }

      IncidentNetIterator operator++(int) {
        IncidentNetIterator tmp_it(_it, _vec);
        ++_it;
        return tmp_it;
      }

      IncidentNetIterator operator--(int) {
        IncidentNetIterator tmp_it(_it, _vec);
        --_it;
        return tmp_it;
      }

      IncidentNetIterator operator+(const difference_type& n) const {
        return IncidentNetIterator(_it + n, _vec);
      }

      IncidentNetIterator& operator+=(const difference_type& n) {
        _it += n;
        return *this;
      }

      IncidentNetIterator operator-(const difference_type& n) const {
        return IncidentNetIterator(_it - n, _vec);
      }

      IncidentNetIterator& operator-=(const difference_type& n) {
        _it -= n;
        return *this;
      }

      reference operator[](const difference_type& n) const {
        return *_it[n];
      }

      bool operator==(const IncidentNetIterator& other) const {
        return _it == other._it;
      }

      bool operator!=(const IncidentNetIterator& other) const {
        return _it != other._it;
      }

      bool operator<(const IncidentNetIterator& other) const {
        return _it < other._it;
      }

      bool operator>(const IncidentNetIterator& other) const {
        return _it > other._it;
      }

      bool operator<=(const IncidentNetIterator& other) const {
        return _it <= other._it;
      }

      bool operator>=(const IncidentNetIterator& other) const {
        return _it >= other._it;
      }

      difference_type operator+(const IncidentNetIterator& other) const {
        return ( _it + other._it );
      }

      difference_type operator-(const IncidentNetIterator& other) const {
        return (_it - other._it);
      }

    private:
      BaseIterator _it;
      IncidentNetVectorType* _vec;
  };

 public:
  using iterator        = IncidentNetIterator<T, typename Base::iterator, IncidentNetVector>;
  using const_iterator  = const IncidentNetIterator<const T, typename Base::const_iterator, const IncidentNetVector>;

  IncidentNetVector() :
    Base(),
    _rw_mutex(),
    _ref_count(0) { }

  IncidentNetVector(const size_t size) :
    Base(size),
    _rw_mutex(),
    _ref_count(0) {  }

  IncidentNetVector(const size_t size, const T init) :
    Base(size, init),
    _rw_mutex(),
    _ref_count(0) {  }

  IncidentNetVector(IncidentNetVector&& other) :
    Base(std::move(other)),
    _rw_mutex(),
    _ref_count(0) { }

  IncidentNetVector & operator=(IncidentNetVector&& other) {
    _rw_mutex = std::shared_timed_mutex();
    _ref_count = parallel::IntegralAtomicWrapper<size_t>(0);
    Base::operator=(std::move(other));
    return *this;
  }

  void bulk_insert(const parallel::scalable_vector<T>& values) {
    std::lock_guard<std::shared_timed_mutex> unique_lock(_rw_mutex);
    // Wait until destructors of all iterators are called
    // Note, at this point bulk_insert function holds unique lock
    // all calls to iterators have to wait.
    while ( _ref_count.load() > 0 ) { }
    for ( const T& value : values ) {
      this->push_back(value);
    }
  }

  IteratorRange<const_iterator> it() {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return IteratorRange<iterator>(
      const_iterator(Base::begin(), this), const_iterator(Base::end(), this));
  }

  IteratorRange<const_iterator> c_it() const {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return IteratorRange<const_iterator>(
      const_iterator(Base::cbegin(), this), const_iterator(Base::cend(), this));
  }

  iterator begin() {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return iterator(Base::begin(), this);
  }

  const_iterator cbegin() const {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return const_iterator(Base::cbegin(), this);
  }

  iterator end() {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return iterator(Base::end(), this);
  }

  const_iterator cend() const {
    std::shared_lock<std::shared_timed_mutex> read_lock(_rw_mutex);
    return const_iterator(Base::cend(), this);
  }

  // ! Only for testing
  size_t active_iterators() const {
    return _ref_count;
  }

 private:
  // ! Read-Write lock to protect iterators from bulk inserts
  mutable std::shared_timed_mutex _rw_mutex;
  // ! Number of active iterators
  mutable parallel::IntegralAtomicWrapper<size_t> _ref_count;
};
}  // namespace ds
}  // namespace mt_kahypar