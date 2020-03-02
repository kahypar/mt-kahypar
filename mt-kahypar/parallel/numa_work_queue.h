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

template<typename Work>
class NumaWorkQueue {
public:
  explicit NumaWorkQueue(size_t numSockets) : queues(numSockets) { }
  explicit NumaWorkQueue() : queues(static_cast<size_t>(TBBNumaArena::instance().num_used_numa_nodes())) { }

  bool empty() const {
    return std::all_of(queues.begin(), queues.end(), [](const auto& q) { return q.empty(); });
  }

  void push(const Work& w, const int socket) {
    queues[socket].push(w);
  }

  bool tryPop(Work& dest, int preferredSocket) {
    if (queues[preferredSocket].try_pop(dest)) {
      return true;
    }
    size_t maxIndex = 0;
    size_t maxSize = 0;
    for (size_t i = 0; i < queues.size(); ++i) {
      size_t size = queues[i].unsafe_size();
      if (size > maxSize) {
        maxSize = size;
        maxIndex = i;
      }
    }
    return queues[maxIndex].try_pop(dest);
  }

  bool tryPop(Work& dest) {
    int socket = HardwareTopology::instance().numa_node_of_cpu(sched_getcpu());
    return tryPop(dest, socket);
  }

  void shuffleQueues() {
    // Implement Me

    // reimplement with own vectors and give up parallel insertion?
  }

private:
  std::vector<tbb::concurrent_queue<Work>> queues;
};

}