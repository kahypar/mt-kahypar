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

namespace mt_kahypar {
namespace ds {

/**
 *
 */
template <typename T>
class IncidentNetVector : public parallel::scalable_vector<T> {

  static constexpr bool debug = false;
  using Base = parallel::scalable_vector<T>;

  using base_iterator = typename Base::iterator;

  template<typename BaseIterator>
  class IncidentNetIterator : public BaseIterator {

    public:
      IncidentNetIterator(BaseIterator it, IncidentNetVector* vec) :
        BaseIterator(it),
        _vec(vec) {
        ++_vec->_ref_count;
      }

      ~IncidentNetIterator() {
        --_vec->_ref_count;
      }

    private:
      IncidentNetVector* _vec;
  };

 public:
  using iterator        = IncidentNetIterator<typename Base::iterator>;
  using const_iterator  = const IncidentNetIterator<typename Base::const_iterator>;

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

 private:
  // ! Read-Write lock to protect iterators from bulk inserts
  std::shared_timed_mutex _rw_mutex;
  // ! Number of active iterators
  parallel::IntegralAtomicWrapper<size_t> _ref_count;
};
}  // namespace ds
}  // namespace mt_kahypar