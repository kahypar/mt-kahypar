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

#include <memory>
#include <type_traits>

#include "tbb/scalable_allocator.h"

namespace mt_kahypar {
namespace parallel {

template<typename T>
struct tbb_deleter {
  void operator()(T *p) {
    scalable_free(p);
  }
};

template<typename T>
using tbb_unique_ptr = std::unique_ptr<T, tbb_deleter<T>>;

template<typename T>
static tbb_unique_ptr<T> make_unique(const size_t size) {
  constexpr bool zero_initialize = !std::is_trivially_copyable<T>::value;
  // For types that are trivially copyable, it is unproblematic to return
  // uninitialized (non-zeroed) memory. However, in case of complex types it can
  // easily happen that the destructor is called on uninitialized data, causing
  // UB in the process.
  // We guard against this by returning zeroed memory in such cases (note: this
  // still might be problematic if zero is not a valid instance of T).
  T* ptr = (T*) (zero_initialize ? scalable_calloc(size, sizeof(T)) : scalable_malloc(sizeof(T) * size));
  return tbb_unique_ptr<T>(ptr, parallel::tbb_deleter<T>());
}

}  // namespace parallel
}  // namespace mt_kahypar
