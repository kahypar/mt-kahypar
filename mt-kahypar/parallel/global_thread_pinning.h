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

#include <mutex>
#include <unordered_map>
#include <thread>
#include <algorithm>
#include <numeric>

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace parallel {

template< typename HwTopology >
class GlobalThreadPinning {

  static constexpr bool debug = false;

 public:
  static GlobalThreadPinning& instance(const int num_threads = std::thread::hardware_concurrency()) {
    static GlobalThreadPinning instance(num_threads);
    return instance;
  }

  void pin_thread() {
    std::lock_guard<std::mutex> lock(_pinning_mutex);
    std::thread::id thread_id = std::this_thread::get_id();
    int cpu_id = register_thread(thread_id);
    if ( !_is_pinned_to_numa_node[thread_id] ) {
      DBG << "Assign thread with PID" << thread_id << "to cpu" << cpu_id;
      pin_thread_to_cpu(cpu_id);
    }
  }

  void unregister_thread() {
    std::lock_guard<std::mutex> lock(_pinning_mutex);
    int cpu_id = sched_getcpu();
    std::thread::id thread_id = std::this_thread::get_id();
    ASSERT(std::find(_free_cpus.begin(), _free_cpus.end(), cpu_id) == _free_cpus.end(), "CPU" << cpu_id << "is already free");
    ASSERT(_pinned_threads.find(thread_id) != _pinned_threads.end(), "Thread wit PID" << thread_id << "is not registered");
    DBG << "Unregister thread with PID" << thread_id << "on cpu" << cpu_id;
    _free_cpus.push_back(cpu_id);
    _pinned_threads.erase(thread_id);
  }

  void pin_thread_to_numa_node( const int node, const int cpu_id ) {
    unused(node);
    std::lock_guard<std::mutex> lock(_pinning_mutex);
    std::thread::id thread_id = std::this_thread::get_id();
    ASSERT(!_is_pinned_to_numa_node[thread_id], "Thread already pinned to a numa node");
    ASSERT(HwTopology::instance().numa_node_of_cpu(cpu_id) == node,
      "CPU" << cpu_id << "is not on numa node" << node << ", actually it is on"
        << HwTopology::instance().numa_node_of_cpu(sched_getcpu()));
    _is_pinned_to_numa_node[thread_id] = true;
    DBG << "Assign thread with PID" << thread_id << "to cpu" << cpu_id
        << "on numa node" << node;
    pin_thread_to_cpu(cpu_id);
  }

  void unpin_thread_from_numa_node( const int node ) {
    unused(node);
    std::lock_guard<std::mutex> lock(_pinning_mutex);
    std::thread::id thread_id = std::this_thread::get_id();
    ASSERT(_is_pinned_to_numa_node[thread_id], "Thread was not pinned to a numa node");
    ASSERT(HwTopology::instance().numa_node_of_cpu(sched_getcpu()) == node,
      "CPU" << sched_getcpu() << "is not on numa node" << node << ", actually it is on"
        << HwTopology::instance().numa_node_of_cpu(sched_getcpu()));
    DBG << "Unassign thread with PID" << thread_id << "on cpu" << sched_getcpu()
        << "from numa node" << node;
    _is_pinned_to_numa_node[thread_id] = false;
    if (is_registered(thread_id)) {
      pin_thread_to_cpu(_pinned_threads[thread_id]);
    }
  }

 private:
  explicit GlobalThreadPinning(const int num_threads) :
    _num_cpus(std::thread::hardware_concurrency()),
    _num_threads(num_threads),
    _pinning_mutex(),
    _free_cpus(_num_cpus),
    _pinned_threads(),
    _is_pinned_to_numa_node() {
    std::iota(_free_cpus.begin(), _free_cpus.end(), 0);

    // Sort cpus in the following order
    // 1.) Non-hyperthread first
    // 2.) Increasing order of numa node
    // 3.) Increasing order of cpu id
    // ...
    HwTopology& topology = HwTopology::instance();
    std::sort(_free_cpus.begin(), _free_cpus.end(),
              [&](const int& lhs, const int& rhs) {
      int node_lhs = topology.numa_node_of_cpu(lhs);
      int node_rhs = topology.numa_node_of_cpu(rhs);
      bool is_hyperthread_lhs = topology.is_hyperthread(lhs);
      bool is_hyperthread_rhs = topology.is_hyperthread(rhs);
      return is_hyperthread_lhs < is_hyperthread_rhs ||
        (is_hyperthread_lhs == is_hyperthread_rhs && node_lhs < node_rhs) ||
        (is_hyperthread_lhs == is_hyperthread_rhs && node_lhs == node_rhs && lhs < rhs);
    });
    // ... this ensure that we first pop nodes in hyperthreading
    while ( (int)_free_cpus.size() > _num_threads ) {
      _free_cpus.pop_back();
    }
  }

  int register_thread(const std::thread::id thread_id) {
    if ( !is_registered(thread_id) ) {
      ASSERT(!_free_cpus.empty(), "There are more threads than CPUs");
      _pinned_threads[thread_id] = _free_cpus.back();
      _free_cpus.pop_back();
      DBG << "Thread with PID" << std::this_thread::get_id()
          << "successfully registered on CPU" << _pinned_threads[thread_id];
    }
    return _pinned_threads[thread_id];
  }

  bool is_registered(const std::thread::id thread_id) {
    return _pinned_threads.find(thread_id) != _pinned_threads.end();
  }

  void pin_thread_to_cpu(const int cpu_id) {
		const size_t size = CPU_ALLOC_SIZE( _num_cpus );
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu_id, &mask);
    const int err = sched_setaffinity(0, size, &mask);

		if ( err ) {
			LOG << "Failed to set thread affinity";
			exit( EXIT_FAILURE );
		}

    ASSERT(sched_getcpu() == cpu_id);
    DBG << "Thread with PID" << std::this_thread::get_id()
        << "successfully pinned to CPU" << cpu_id
        << "( Currently =" << V(sched_getcpu()) << ")";
  }

  static std::mutex _mutex;

  const int _num_cpus;
  const int _num_threads;
  std::mutex _pinning_mutex;
  std::vector<int> _free_cpus;
  std::unordered_map<std::thread::id, int> _pinned_threads;
  std::unordered_map<std::thread::id, bool> _is_pinned_to_numa_node;
};

template< typename HwTopology >
std::mutex GlobalThreadPinning<HwTopology>::_mutex;

} // namespace parallel
} // namespace mt_kahypar
