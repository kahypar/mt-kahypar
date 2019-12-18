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

#include <unordered_map>
#include <thread>
#include <mutex>

#undef __TBB_ARENA_OBSERVER
#define __TBB_ARENA_OBSERVER true
#include "tbb/task_scheduler_observer.h"
#undef __TBB_ARENA_OBSERVER

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/global_thread_pinning.h"

namespace mt_kahypar {
namespace parallel {
/**
 * Pins threads of task arena to a NUMA node. Each time a thread
 * enters a task arena on_scheduler_entry(...) is called. Each time
 * a thread leaves a task arena on_scheduler_exit(...) is called.
 */
template <typename HwTopology>
class NumaThreadPinningObserver : public tbb::task_scheduler_observer {
  using Base = tbb::task_scheduler_observer;
  using GlobalThreadPinning = mt_kahypar::parallel::GlobalThreadPinning<HwTopology>;

  static constexpr bool debug = false;

 public:
  explicit NumaThreadPinningObserver(tbb::task_arena& arena,
                                     int numa_node,
                                     const bool oversubscription_allowed) :
    Base(arena),
    _topology(HwTopology::instance()),
    _numa_node(numa_node),
    _oversubscription_allowed(oversubscription_allowed),
    _unsuccessful_pinning_mutex(),
    _is_unsuccessful_pinned_to_numa_node(),
    _cpus(HwTopology::instance().get_cpus_of_numa_node(numa_node)),
    _free_cpus(_cpus.size()) {
    observe(true);
  }

  NumaThreadPinningObserver(const NumaThreadPinningObserver&) = delete;
  NumaThreadPinningObserver & operator= (const NumaThreadPinningObserver &) = delete;

  NumaThreadPinningObserver(NumaThreadPinningObserver&& other) :
    _topology(other._topology),
    _numa_node(other._numa_node),
    _oversubscription_allowed(other._oversubscription_allowed),
    _unsuccessful_pinning_mutex(),
    _is_unsuccessful_pinned_to_numa_node(std::move(other._is_unsuccessful_pinned_to_numa_node)),
    _cpus(std::move(other._cpus)),
    _free_cpus(other._free_cpus) { }

  NumaThreadPinningObserver & operator= (NumaThreadPinningObserver &&) = delete;

  ~NumaThreadPinningObserver() {
    observe(false);
  }

  void on_scheduler_entry(bool) override {
    if( !_topology.pin_thread_to_numa_node(_numa_node) ) {
      // In case, pinning the threads to a numa node fails, we have
      // oversubsribed the numa node. This can happen during recursive initial
      // partitioning, because we split there the TBB Numa Arena and guarantee
      // that at least one slot is available on each numa node. To resolve
      // that conflict, we pin the thread to an other cpu here (we possibly
      // oversubscribe that here, but this should happen only rarely and on
      // small subproblems)

      if ( !_oversubscription_allowed ) {
        ERROR("Numa node" << _numa_node << "is oversubscribed with too many threads");
      }

      std::lock_guard<std::mutex> lock(_unsuccessful_pinning_mutex);
      if ( _free_cpus == 0 ) {
        _free_cpus = _cpus.size();
      }
      int cpu_id = _cpus.front();
      std::swap(_cpus.front(), _cpus[--_free_cpus]);
      _is_unsuccessful_pinned_to_numa_node[std::this_thread::get_id()] = true;
      GlobalThreadPinning::instance().pin_thread_to_numa_node(_numa_node, cpu_id);
    }
  }

  void on_scheduler_exit(bool) override {
    std::lock_guard<std::mutex> lock(_unsuccessful_pinning_mutex);
    if ( _is_unsuccessful_pinned_to_numa_node[std::this_thread::get_id()] ) {
      GlobalThreadPinning::instance().unpin_thread_from_numa_node(_numa_node);
    } else {
      _topology.unpin_thread_from_numa_node(_numa_node);
    }
  }

 private:
  HwTopology& _topology;
  int _numa_node;
  const bool _oversubscription_allowed;

  std::mutex _unsuccessful_pinning_mutex;
  std::unordered_map<std::thread::id, bool> _is_unsuccessful_pinned_to_numa_node;
  std::vector<int> _cpus;
  size_t _free_cpus;
};
}  // namespace parallel
}  // namespace mt_kahypar
