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
  CAtomic<size_t> front;

  ThreadQueue() {
    elements.reserve(1 << 13);
    front.store(0);
  }

  void clear() {
    elements.clear();
    front.store(0);
  }

  bool try_pop(T& dest) {
    size_t slot = front.fetch_add(1, std::memory_order_acq_rel);
    if (slot < elements.size()) {
      dest = elements[slot];
      return true;
    }
    return false;
  }
};

template<typename T>
struct WorkContainer {

  using TimestampT = uint32_t;

  WorkContainer(size_t maxNumElements = 0, size_t maxNumThreads = 0) :
          timestamps(maxNumElements, 0),
          tls_queues(maxNumThreads)
  {

  }

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
    ASSERT(el < timestamps.size());
    timestamps[el] = current;
  }

  // assumes that no thread is currently calling try_pop
  void safe_push(const T el, size_t thread_id) {
    ASSERT(thread_id < tls_queues.size());
    tls_queues[thread_id].elements.push_back(el);
    ASSERT(tls_queues[thread_id].front.load() == 0);
    ASSERT(el < timestamps.size());
    timestamps[el] = current;
  }

  bool try_pop(T& dest, size_t thread_id) {
    ASSERT(thread_id < tls_queues.size());
    return tls_queues[thread_id].try_pop(dest) || conc_queue.try_pop(dest) || steal_work(dest);
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
    ASSERT(el < timestamps.size());
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
  vec<ThreadQueue<T>> tls_queues;
  tbb::concurrent_queue<T> conc_queue;

  using SubRange = IteratorRange< typename vec<T>::const_iterator >;
  using Range = ConcatenatedRange<SubRange>;

  Range safely_inserted_range() const {
    Range r;
    for (const ThreadQueue<T>& q : tls_queues) {
      r.concat( SubRange(q.elements.cbegin(), q.elements.cend()) );
    }
    return r;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* work_container_node = parent->addChild("Work Container");
    utils::MemoryTreeNode* local_work_queue_node = work_container_node->addChild("Local Work Queue");
    for (const ThreadQueue<T>& q : tls_queues) {
      local_work_queue_node->updateSize(q.elements.capacity() * sizeof(T));
    }
    utils::MemoryTreeNode* timestamps_node = work_container_node->addChild("Timestamps");
    timestamps_node->updateSize(timestamps.capacity() * sizeof(TimestampT));
  }
};


}