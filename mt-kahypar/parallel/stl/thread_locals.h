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

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace mt_kahypar {

  namespace internals {
    template<typename T>
    using ThreadLocal = tbb::enumerable_thread_specific<T>;

    template<typename T, typename F>
    struct ThreadLocalFree {
      using RangeType = typename ThreadLocal<T>::range_type;
      using Iterator = typename ThreadLocal<T>::iterator;

      explicit ThreadLocalFree(F&& free_func) :
              _free_func(free_func) { }

      void operator()(RangeType& range) const {
        for ( Iterator it = range.begin(); it < range.end(); ++it ) {
          _free_func(*it);
        }
      }

      F _free_func;
    };
  } // namespace

  namespace parallel {
    template<typename T, typename F>
    static void parallel_free_thread_local_internal_data(internals::ThreadLocal<T>& local,
                                                         F&& free_func) {
      internals::ThreadLocalFree<T, F> thread_local_free(std::move(free_func));
      tbb::parallel_for(local.range(), thread_local_free);
    }
  }

  template<typename T>
  using tls_enumerable_thread_specific = tbb::enumerable_thread_specific<T, tbb::cache_aligned_allocator<T>, tbb::ets_key_per_instance>;

}