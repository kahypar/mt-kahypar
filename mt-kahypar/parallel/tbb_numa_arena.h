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

  using ThreadPinningObserver = mt_kahypar::parallel::ThreadPinningObserver<HwTopology>;

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

  hwloc_cpuset_t used_cpuset() const {
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    for ( const auto& numa_node : _numa_node_to_cpu_id ) {
      for ( const int cpu_id : numa_node ) {
        hwloc_bitmap_set(cpuset, cpu_id);
      }
    }
    return cpuset;
  }

  RecursionTaskGroups create_tbb_task_groups_for_recursion() {
    return std::make_pair(0, 0);
  }


  void terminate() {

    if ( _global_observer ) {
      _global_observer->observe(false);
    }

    if ( _init ) {
      _init->terminate();
    }
  }

  void initialize(int num_threads) {
    _num_threads = num_threads;

    if ( _init ) {
      _init->initialize(num_threads);
    } else {
      _init = std::make_unique<tbb::task_scheduler_init>(num_threads);
    }

    if ( _global_observer ) {
      assert(num_threads <= _cpus.size());
      _global_observer->observe(true);
    } else {
      HwTopology& topology = HwTopology::instance();
      int num_numa_nodes = topology.num_numa_nodes();
      DBG << "Initialize TBB with" << num_threads << "threads";

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
    }

  }

 private:
  explicit TBBNumaArena(const int num_threads) :
    _init(nullptr),
    _global_observer(nullptr)
  {
    initialize(num_threads);
  }

  int _num_threads;
  std::unique_ptr<tbb::task_scheduler_init> _init;
  std::unique_ptr<ThreadPinningObserver> _global_observer;
  std::shared_timed_mutex _task_group_read_write_mutex;
  std::vector<int> _cpus;
  std::vector<std::vector<int>> _numa_node_to_cpu_id;
};

template <typename HwTopology, bool is_numa_aware>
typename TBBNumaArena<HwTopology, is_numa_aware>::TaskGroupID TBBNumaArena<HwTopology, is_numa_aware>::GLOBAL_TASK_GROUP = 0;

}  // namespace parallel
}  // namespace mt_kahypar
