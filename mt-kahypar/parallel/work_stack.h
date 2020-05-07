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
struct SingleProducerMultipleConsumerDeque {
  static constexpr size_t in_reallocation = std::numeric_limits<size_t>::max() / 2;

  SingleProducerMultipleConsumerDeque() {
    clear();
    elements.reserve(1 << 13);
  }

  vec<T> elements;
  CAtomic<size_t> front, back;

  void clear() {
    front.store(0);
    back.store(0);
    elements.clear();
  }

  void push_back(const T& el) {
    if (elements.size() == elements.capacity()) {
      size_t old_front = front.load(std::memory_order_acq_rel);
      while (front.compare_exchange_weak(old_front, in_reallocation, std::memory_order_acq_rel)) { /* spin */ }
      elements.push_back(el); // causes reallocation
      back.fetch_add(1, std::memory_order_acq_rel);
      front.store(old_front, std::memory_order_acq_rel);
    } else {
      elements.push_back(el);
      back.fetch_add(1, std::memory_order_acq_rel);
    }
  }

  // reserved for the owning thread
  bool pop_back(T& el) {
    if (!empty()) {
      back.fetch_sub(1, std::memory_order_acq_rel);
      el = elements.pop_back();
      return true;
    }
    return false;
  }

  // designated for other threads
  // in our scenario it is not problematic if the same element gets popped by the producer and one consumer
  // since we have a secondary locking mechanism in place for them
  bool pop_front(T& el) {
    if (front.load(std::memory_order_acq_rel) < back.load(std::memory_order_acq_rel)) { // extra check observes the spin lock
      size_t slot = front.fetch_add(1, std::memory_order_acq_rel);
      if (slot < back.load(std::memory_order_acq_rel)) {
        el = elements[slot];
        return true;
      }
    }
    return false;
  }

  bool currently_blocked() const {
    return front.load(std::memory_order_acq_rel) >= std::numeric_limits<size_t>::max() / 2;
  }

  size_t unsafe_size() const {
    const size_t f = front.load(std::memory_order_acq_rel), b = back.load(std::memory_order_acq_rel);
    return b >= f ? b - f : 0;
  }

  bool empty() const {
    return unsafe_size() > 0;
  }
};

template<typename T>
struct WorkStealingContainer {

  using TimestampT = uint32_t;

  WorkStealingContainer(size_t maxNumElements) : timestamps(maxNumElements, 0) { }

  void push_back(const T el) {
    tls_deques.local().push_back(el);
    timestamps[el] = current;
  }

  bool try_pop(T& dest) {
    if (tls_deques.local().pop_back(dest)) {
      timestamps[dest] = current+1;
      return true;
    } else {
      // try stealing
      for (SingleProducerMultipleConsumerDeque<T>& other : tls_deques) {
        if (other.pop_front(dest)) {
          timestamps[dest] = current+1;
          return true;
        }
      }

      size_t fails = steal_failures.fetch_add(1, std::memory_order_relaxed);
      if (fails < 1024) {
        // stealing failed --> check if any of them are currently blocked. if so spin. otherwise return false
        for (SingleProducerMultipleConsumerDeque<T>& other : tls_deques) {
          if (other.currently_blocked()) {
            while (other.currently_blocked()) { /* spin */ }
            if (other.pop_front(dest)) {
              timestamps[dest] = current+1;
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  bool was_pushed_and_removed(const T el) const {
    return timestamps[el] == current+1;
  }

  void balance() {

  }

  void shuffle() {
    tbb::parallel_for_each(tls_deques, [&](SingleProducerMultipleConsumerDeque<T>& tlq) {
      std::shuffle(tlq.elements.begin() + tlq.front, tlq.elements.begin() + tlq.back, utils::Randomize::instance().getGenerator());
    });
  }

  void clear() {
    if (current >= std::numeric_limits<TimestampT>::max() - 2) {
      tbb::parallel_for_each(timestamps, [](TimestampT& x) { x = 0; });
      current = 0;
    }
    current += 2;
    steal_failures.store(0, std::memory_order_relaxed);
  }

  TimestampT current = 2;
  vec<TimestampT> timestamps;
  CAtomic<size_t> steal_failures { 0 };
  tls_enumerable_thread_specific<SingleProducerMultipleConsumerDeque<T>> tls_deques;
  //tls_enumerable_thread_specific<tbb::concurrent_queue<T>> tls_queues;

  /*

  WorkStealingContainer(size_t maxNumElements) :
          size(0),
          elements(maxNumElements, T())
  {

  }



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

  void shuffle() {
    utils::Randomize::instance().parallelShuffleVector(elements, 0, unsafe_size());
  }
  */
};


}