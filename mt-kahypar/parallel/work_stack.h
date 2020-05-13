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
struct ThreadQueue {
  vec<T> elements;
  CAtomic<size_t> front { 0 };

  ThreadQueue() {
    elements.reserve(1 << 13);
  }

  void clear() {
    elements.clear();
    front.store(0);
  }

  bool try_pop(T& dest) {
    size_t slot = front.fetch_add(1, std::memory_order_acq_rel);
    if (slot < elements.size()) {
      dest= elements[slot];
      return true;
    }
    return false;
  }
};

template<typename T>
struct WorkContainer {

  using TimestampT = uint32_t;

  WorkContainer(size_t maxNumElements) :
          timestamps(maxNumElements, 0) { }

  size_t unsafe_size() const {
    size_t sz = 0;
    for (const ThreadQueue<T>& q : tls_queues) {
      sz += q.elements.size() - q.front.load(std::memory_order_relaxed);
    }
    sz += conc_queue.unsafe_size();
    return sz;
  }

  void concurrent_push(const T el) {
    conc_queue.push(el);
    timestamps[el] = current;
  }

  // assumes that no thread is currently calling try_pop
  void safe_push(const T el) {
    tls_queues.local().elements.push_back(el);
    timestamps[el] = current;
  }

  bool try_pop(T& dest) {
    return tls_queues.local().try_pop(dest) || conc_queue.try_pop(dest) || steal_work(dest);
  }

  bool steal_work(T& dest) {
    for (ThreadQueue<T>& q : tls_queues) {
      if (q.try_pop(dest)) {
        return true;
      }
    }
    return false;
  }

  bool was_pushed_and_removed(const T el) const {
    return timestamps[el] == current+1;
  }

  void shuffle() {
    tbb::parallel_for_each(tls_queues, [&](ThreadQueue<T>& q) {
      utils::Randomize::instance().shuffleVector(q.elements);
    });
  }

  void clear() {
    if (current >= std::numeric_limits<TimestampT>::max() - 2) {
      tbb::parallel_for_each(timestamps, [](TimestampT& x) { x = 0; });
      current = 0;
    }
    for (ThreadQueue<T>& q : tls_queues) {
      q.clear();
    }
    conc_queue.clear();
    current += 2;
  }

  TimestampT current = 2;
  vec<TimestampT> timestamps;
  tls_enumerable_thread_specific<ThreadQueue<T>> tls_queues;
  tbb::concurrent_queue<T> conc_queue;
};


}