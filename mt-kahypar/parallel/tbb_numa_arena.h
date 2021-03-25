/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <hwloc.h>
#include <mutex>
#include <memory>
#include <shared_mutex>
#include <functional>

#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/thread_pinning_observer.h"

namespace mt_kahypar {
namespace parallel {
/**
 * Creates number of NUMA nodes TBB task arenas. Each task arena is pinned
 * to a unique NUMA node. Each task arena can then be used to execute tasks
 * on specific NUMA node.
 */
template <typename HwTopology, bool is_numa_aware>
class TBBNumaArena {

  static constexpr bool debug = false;

  struct MovableTaskGroup {
    MovableTaskGroup() :
      task_group() { }

    MovableTaskGroup(MovableTaskGroup&&) :
      task_group() { }

    tbb::task_group task_group;
  };

  using Self = TBBNumaArena<HwTopology, is_numa_aware>;
  using ThreadPinningObserver = mt_kahypar::parallel::ThreadPinningObserver<HwTopology>;
  using NumaTaskGroups = std::vector<MovableTaskGroup>;

 public:
  using TaskGroupID = size_t;
  using RecursionTaskGroups = std::pair<TaskGroupID, TaskGroupID>;
  static TaskGroupID GLOBAL_TASK_GROUP;

  TBBNumaArena(const TBBNumaArena&) = delete;
  TBBNumaArena & operator= (const TBBNumaArena &) = delete;

  TBBNumaArena(TBBNumaArena&&) = delete;
  TBBNumaArena & operator= (TBBNumaArena &&) = delete;

  static TBBNumaArena& instance(const size_t num_threads = std::thread::hardware_concurrency()) {
    static TBBNumaArena instance(num_threads);
    return instance;
  }

  int total_number_of_threads() const {
    return _num_threads;
  }

  int number_of_used_cpus_on_numa_node(const int node) const {
    ASSERT(static_cast<size_t>(node) < _numa_node_to_cpu_id.size());
    return _numa_node_to_cpu_id[node].size();
  }

  int num_used_numa_nodes() const {
    return _numa_node_to_cpu_id.size();
  }

  int num_numa_arenas() const {
    return _arenas.size();
  }

  hwloc_cpuset_t used_cpuset() const {
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    for ( const auto& numa_node : _numa_node_to_cpu_id ) {
      for ( const int cpu_id : numa_node ) {
        hwloc_bitmap_set(cpuset, cpu_id);
      }
    }
    return cpuset;
  }

  tbb::task_arena& numa_task_arena(const int node) {
    ASSERT(static_cast<size_t>(node) < _arenas.size(), V(node) << V(_arenas.size()));
    return _arenas[node];
  }

  tbb::task_group& numa_task_group(const TaskGroupID task_group_id, const int node) {
    std::shared_lock<std::shared_timed_mutex> read_lock(_task_group_read_write_mutex);
    ASSERT(static_cast<size_t>(node) <= _groups[task_group_id].size());
    return _groups[task_group_id][node].task_group;
  }

  RecursionTaskGroups create_tbb_task_groups_for_recursion() {
    std::lock_guard<std::shared_timed_mutex> write_lock(_task_group_read_write_mutex);
    TaskGroupID task_group_1 = _groups.size();
    _groups.emplace_back();
    TaskGroupID task_group_2 = _groups.size();
    _groups.emplace_back();
    for ( int node = 0; node < num_numa_arenas(); ++node ) {
      _groups[task_group_1].emplace_back();
      _groups[task_group_2].emplace_back();
    }
    return std::make_pair(task_group_1, task_group_2);
  }

  template <typename F>
  void execute_parallel_on_all_numa_nodes(const TaskGroupID task_group_id, F&& func) {
    for (int node = 0; node < num_numa_arenas(); ++node) {
      numa_task_arena(node).execute([&, node] {
            numa_task_group(task_group_id, node).run([&, node] {
              func(node);
            });
          });
    }
    wait(task_group_id);
  }

  void wait(const TaskGroupID task_group_id) {
    for (int node = 0; node < num_numa_arenas(); ++node) {
      _arenas[node].execute([&, node] {
            numa_task_group(task_group_id, node).wait();
      });
    }
  }

  void wait(const int node, tbb::task_group& group) {
    ASSERT(static_cast<size_t>(node) < _arenas.size());
    _arenas[node].execute([&] {
          group.wait();
    });
  }

  void terminate() {
    for (tbb::task_arena& arena : _arenas) {
      arena.terminate();
    }

    for (ThreadPinningObserver& observer : _observers) {
      observer.observe(false);
    }

    if ( _global_observer ) {
      _global_observer->observe(false);
    }

    if ( _init ) {
      _init->terminate();
    }
  }

  void initialize(int num_threads) {
    _num_threads = num_threads;

    if (_init) {
      _init->initialize(num_threads);
    } else {
      _init = std::make_unique<tbb::task_scheduler_init>(num_threads);
    }

    HwTopology& topology = HwTopology::instance();
    int num_numa_nodes = topology.num_numa_nodes();
    DBG << "Initialize TBB with" << num_threads << "threads";
    _arenas.reserve(num_numa_nodes);
    // TODO(heuer): fix copy constructor of observer
    _observers.reserve(num_numa_nodes);

    _cpus = topology.get_all_cpus();
    // Sort cpus in the following order
    // 1.) Non-hyperthread first
    // 2.) Increasing order of numa node
    // 3.) Increasing order of cpu id
    // ...
    std::sort(_cpus.begin(), _cpus.end(), [&](const int& lhs, const int& rhs) {
                int node_lhs = topology.numa_node_of_cpu(lhs);
                int node_rhs = topology.numa_node_of_cpu(rhs);
                bool is_hyperthread_lhs = topology.is_hyperthread(lhs);
                bool is_hyperthread_rhs = topology.is_hyperthread(rhs);
                return std::tie(is_hyperthread_lhs, node_lhs, lhs)
                        < std::tie(is_hyperthread_rhs, node_rhs, rhs);
    });
    // ... this ensures that we first pop nodes in hyperthreading
    while (static_cast<int>(_cpus.size()) > _num_threads) {
      _cpus.pop_back();
    }
    _global_observer = std::make_unique<ThreadPinningObserver>(_cpus);

    _numa_node_to_cpu_id.clear();
    _numa_node_to_cpu_id.resize(num_numa_nodes);
    for ( const int cpu_id : _cpus ) {
      int node = topology.numa_node_of_cpu(cpu_id);
      ASSERT(node < static_cast<int>(_numa_node_to_cpu_id.size()));
      _numa_node_to_cpu_id[node].push_back(cpu_id);
    }
    while( !_numa_node_to_cpu_id.empty() && _numa_node_to_cpu_id.back().empty() ) {
      _numa_node_to_cpu_id.pop_back();
    }

    _arenas.clear();
    _groups.resize(1);
    _observers.clear();
    if ( is_numa_aware ) {
      for ( size_t node = 0; node < _numa_node_to_cpu_id.size(); ++node ) {
        initialize_tbb_numa_arena(node, _numa_node_to_cpu_id[node]);
      }
    } else {
      initialize_tbb_numa_arena(0, _cpus);
    }
  }

 private:
  explicit TBBNumaArena(const int num_threads) :
    _init(nullptr),
    _global_observer(nullptr) {
    initialize(num_threads);
  }

  void initialize_tbb_numa_arena(const int node,
                                 const std::vector<int>& cpus_on_numa_node) {
    unused(node);
    _groups.reserve(1024);
    int num_cpus_on_numa_node = cpus_on_numa_node.size();
    if (num_cpus_on_numa_node > 0) {
      DBG << "Initialize TBB task arena on numa node" << node
          << "with" << num_cpus_on_numa_node << "threads";
      #ifndef KAHYPAR_TRAVIS_BUILD
      _arenas.emplace_back(num_cpus_on_numa_node, 0);
      #else
      _arenas.emplace_back(num_cpus_on_numa_node, 1 /* reserve for master */);
      #endif
      _groups[GLOBAL_TASK_GROUP].emplace_back();
      _arenas.back().initialize();
      _observers.emplace_back(_arenas.back(), node, cpus_on_numa_node);
    }
  }

  int _num_threads;
  std::unique_ptr<tbb::task_scheduler_init> _init;
  std::unique_ptr<ThreadPinningObserver> _global_observer;
  std::vector<tbb::task_arena> _arenas;
  std::shared_timed_mutex _task_group_read_write_mutex;
  std::vector<NumaTaskGroups> _groups;
  std::vector<ThreadPinningObserver> _observers;
  std::vector<int> _cpus;
  std::vector<std::vector<int>> _numa_node_to_cpu_id;
};

template <typename HwTopology, bool is_numa_aware>
typename TBBNumaArena<HwTopology, is_numa_aware>::TaskGroupID TBBNumaArena<HwTopology, is_numa_aware>::GLOBAL_TASK_GROUP = 0;

}  // namespace parallel
}  // namespace mt_kahypar
