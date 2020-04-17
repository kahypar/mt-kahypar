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
#include <vector>

#include "../definitions.h"

namespace mt_kahypar {

// TODO better name
// Supports a sequence of concurrent push_back followed by a sequence of concurrent try_pop but not both at the same time
template<typename T>
class ConcurrentDataContainer {
public:

  ConcurrentDataContainer(size_t maxNumElements) : size(0), elements(maxNumElements, T()) { }

  void push_back(const T& el) {
    const size_t old_size = size.fetch_add(1, std::memory_order_acq_rel);
    assert(old_size < elements.size());
    elements[old_size] = el;
  }

  bool try_pop(T& dest) {
    const size_t old_size = size.fetch_sub(1, std::memory_order_acq_rel);
    if (old_size > 0 && old_size < capacity()) {
      dest = elements[old_size - 1];
      return true;
    }
    return false;
  }

  bool empty() const {
    const size_t s = unsafe_size();
    return s == 0 || s >= capacity();
  }

  size_t unsafe_size() const {
    return size.load(std::memory_order_acq_rel);
  }

  size_t capacity() const {
    return elements.size();
  }

  vec<T>& get_underlying_container() {
    return elements;
  }

  void clear() {
    size.store(0);
  }

  void shrink_to_fit() {
    elements.resize(size);
    elements.shrink_to_fit();
  }

private:
  CAtomic<size_t> size;
  vec<T> elements;
};

template<typename Work>
class NumaWorkQueue {
public:
  explicit NumaWorkQueue(size_t numSockets, size_t maxNumElements) : queues(numSockets, ConcurrentDataContainer<Work>(maxNumElements)) { }
  explicit NumaWorkQueue(size_t maxNumElements) : NumaWorkQueue(static_cast<size_t>(TBBNumaArena::instance().num_used_numa_nodes()), maxNumElements) { }

  bool empty() const {
    return std::all_of(queues.begin(), queues.end(), [](const auto& q) { return q.empty(); });
  }

  void push(const Work& w, const int socket) {
    queues[socket].push_back(w);
  }

  bool tryPop(Work& dest, int preferredSocket) {
    if (queues[preferredSocket].try_pop(dest)) {
      return true;
    }
    while (!empty()) {
      // steal from the largest queue
      size_t maxIndex = 0;
      size_t maxSize = 0;
      for (size_t i = 0; i < queues.size(); ++i) {
        size_t size = queues[i].unsafe_size();
        if (!queues[i].empty() && size > maxSize) {
          maxSize = size;
          maxIndex = i;
        }
      }
      if (queues[maxIndex].try_pop(dest)) {
        return true;
      }
    }

    return false;
  }

  bool tryPop(Work& dest) {
    int socket = HardwareTopology::instance().numa_node_of_cpu(sched_getcpu());
    return tryPop(dest, socket);
  }

  void clear() {
    for (auto& q : queues) q.clear();
  }

  size_t unsafe_size() const {
    size_t s = 0;
    for (size_t i = 0; i < queues.size(); ++i) s += queues[i].unsafe_size();
    return s;
  }

  void shuffleQueues() {
    tbb::parallel_for(0UL, queues.size(), [&](size_t i) {
      std::mt19937 rng(queues[i].unsafe_size() + i);
      auto& data = queues[i].get_underlying_container();
      std::shuffle(data.begin(), data.end(), rng);  // TODO parallelize?
    });
  }

private:
  vec<ConcurrentDataContainer<Work>> queues;

};

}