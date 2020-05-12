/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <tbb/concurrent_queue.h>
#include "../definitions.h"

namespace mt_kahypar {

template<typename T>
struct WorkContainer {

  using TimestampT = uint32_t;
  using Queue = tbb::concurrent_queue<T>;

  WorkContainer(size_t maxNumElements) : timestamps(maxNumElements, 0) { }

  size_t unsafe_size() const {
    size_t sz = 0;
    /*
    for (const Queue& x : tls_queues) {
      sz += x.unsafe_size();
    }
     */
    sz = q.unsafe_size();
    return sz;
  }

  void push(const T el) {
    //tls_queues.local().push(el);
    q.push(el);
    timestamps[el] = current;
  }

  bool try_pop(T& dest) {
    return q.try_pop(dest);
    /*
    // use pop_front even on the thread local queue to avoid immediately reusing a just released node
    if (tls_queues.local().try_pop(dest)) {
      timestamps[dest] = current+1;
      return true;
    } else {

      // try stealing
      for (Queue& other_queue : tls_queues) {
        if (other_queue.try_pop(dest)) {
          timestamps[dest] = current+1;
          return true;
        }
      }
    }
    return false;
     */
  }

  bool was_pushed_and_removed(const T el) const {
    return timestamps[el] == current+1;
  }

  void shuffle() {

  }

  void clear() {
    if (current >= std::numeric_limits<TimestampT>::max() - 2) {
      tbb::parallel_for_each(timestamps, [](TimestampT& x) { x = 0; });
      current = 0;
    }
    q.clear();
    /*
    for (Queue& tlq : tls_queues) {
      tlq.clear();
    }
     */
    current += 2;
    steal_failures.store(0, std::memory_order_relaxed);
  }

  TimestampT current = 2;
  vec<TimestampT> timestamps;
  CAtomic<size_t> steal_failures { 0 };
  //tls_enumerable_thread_specific<Queue> tls_queues;
  Queue q;
};


}