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

#include "tbb/task_scheduler_observer.h"

#include "kahypar/macros.h"

namespace kahypar {
namespace parallel {


template< typename HwTopology >
class GlobalThreadPinningObserver : public tbb::task_scheduler_observer {
  using Base = tbb::task_scheduler_observer;

  static constexpr bool debug = true;

 public:
  explicit GlobalThreadPinningObserver(size_t num_threads) :
    Base(true),
    _topology(HwTopology::instance()),
    _num_threads(num_threads),
    _mutex(),
    _cpus(),
    _free_cpus(num_threads) {

    for ( int node = 0; node < (int) _topology.num_numa_nodes(); ++node ) {
      std::vector<int> numa_cpus = _topology.get_cpus_of_numa_node(node);
      _cpus.insert(_cpus.begin(), numa_cpus.begin(), numa_cpus.end());
    }
    std::sort(_cpus.begin(), _cpus.end());
    while ( _cpus.size() > num_threads ) {
      _cpus.pop_back();
    }
    _free_cpus = _cpus.size();

    observe(true);
  }

  GlobalThreadPinningObserver(const GlobalThreadPinningObserver&) = delete;
  GlobalThreadPinningObserver& operator= (const GlobalThreadPinningObserver&) = delete;

  GlobalThreadPinningObserver(GlobalThreadPinningObserver&& other) = default;
  GlobalThreadPinningObserver& operator= (GlobalThreadPinningObserver&&) = delete;

  ~GlobalThreadPinningObserver() {
    observe(false);
  }

  void on_scheduler_entry(bool) override {
    if ( _free_cpus == 0 ) {
      LOG << "More threads than logical cpus";
      exit(EXIT_FAILURE);
    }

    std::lock_guard lock(_mutex);
    int cpu_id = _cpus.front();
    std::swap(_cpus[0], _cpus[--_free_cpus]);
    _topology.pin_thread_to_cpu(cpu_id);

    DBG << "Assigned thread with PID" << std::this_thread::get_id()
        << "to cpu" << cpu_id;
  }

  void on_scheduler_exit(bool) override {
    std::lock_guard lock(_mutex);
    int cpu_id = sched_getcpu();
    size_t pos = _free_cpus;
    while ( pos < _cpus.size() ) {
      if ( _cpus[pos] == cpu_id ) {
        break;
      }
    }
    ASSERT(pos >= _free_cpus && pos < _cpus.size());
    ASSERT(_cpus[pos] == cpu_id);
    std::swap(_cpus[pos], _cpus[_free_cpus++]);
    DBG << "Free thread with PID" << std::this_thread::get_id()
        << "on cpu" << cpu_id;
  }

 private:
  HwTopology& _topology;
  size_t _num_threads;
  std::mutex _mutex;
  std::vector<int> _cpus;
  size_t _free_cpus;
};

} // namespace parallel
} // namespace kahypar