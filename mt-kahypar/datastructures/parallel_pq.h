/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/
#pragma once

#include <vector>
#include <queue>

#include <tbb/enumerable_thread_specific.h>

#include "pcg_random.hpp"

#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

namespace ds {

// a simple parallel priority queue with relaxed ordering guarantee,
// which assumes that insertion and deletion happen in separate phases
// (and no indexing/update functionality)
template<typename KeyT, typename IdT, typename Comparator = std::less<KeyT>>
class MultiQueue {
  using PQ = MaxNAHeap<KeyT, IdT>;

 public:
  MultiQueue(size_t num_threads):
    _comp(),
    _pqs(),
    _locks(),
    _local_rand([]() { return pcg32(555 + THREAD_ID); }),
    _dist(0, 2 * num_threads - 1) {
    const size_t num_pqs = 2 * num_threads;
    _pqs.resize(num_pqs);
    _locks.resize(num_pqs);
  }

  void insert(KeyT key, IdT id) {
    pcg32& rand = _local_rand.local();
    size_t pq_id;
    while (true) {
      pq_id = _dist(rand);
      if (_locks[pq_id].tryLock()) {
        break;
      }
    }
    _pqs[pq_id].insert(id, key);
    _locks[pq_id].unlock();
  }

  bool tryPop(IdT& out_param) {
    static constexpr size_t INVALID = std::numeric_limits<size_t>::max();
    static constexpr size_t NUM_TRIES = 32;

    pcg32& rand = _local_rand.local();
    for (size_t i = 0; i < NUM_TRIES; ++i) {
      size_t first_id = _dist(rand);
      size_t second_id = _dist(rand);
      while (first_id == second_id) {
        second_id = _dist(rand);
      }
      auto& first_pq = _pqs[first_id];
      auto& second_pq = _pqs[second_id];

      if (first_pq.empty() && second_pq.empty()) continue;
      size_t best_id = first_id;
      if (first_pq.empty() || ( !second_pq.empty() && _comp(first_pq.topKey(), second_pq.topKey()) )) {
        best_id = second_id;
      }

      if (!_locks[best_id].tryLock()) continue;
      if (_pqs[best_id].empty()) {
        _locks[best_id].unlock();
        continue;
      }
      out_param = _pqs[best_id].top();
      _pqs[best_id].deleteTop();
      _locks[best_id].unlock();
      return true;
    }

    while (true) {
      KeyT best_key;
      size_t best_id = INVALID;
      for (size_t i = 0; i < _pqs.size(); ++i) {
        if (!_pqs[i].empty() && ( best_id == INVALID || _comp(best_key, _pqs[i].topKey()) )) {
          best_key =_pqs[i].topKey();
          best_id = i;
        }
      }
      if (best_id == INVALID) return false;

      if (!_locks[best_id].tryLock()) continue;
      if (_pqs[best_id].empty()) {
        _locks[best_id].unlock();
        continue;
      }
      out_param = _pqs[best_id].top();
      _pqs[best_id].deleteTop();
      _locks[best_id].unlock();
      return true;
    }
  }

  void reset() {
    for (PQ& pq: _pqs) {
      pq.clear();
    }
  }

 private:
  Comparator _comp;
  vec<PQ> _pqs;
  vec<SpinLock> _locks;
  tbb::enumerable_thread_specific<pcg32> _local_rand;
  std::uniform_int_distribution<size_t> _dist;
};

}
}
