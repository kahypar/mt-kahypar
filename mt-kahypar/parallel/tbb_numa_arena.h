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
#include <memory>

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
template <typename HwTopology>
class TBBNumaArena {
  static constexpr bool debug = false;

 private:
  using Self = TBBNumaArena<HwTopology>;
  using ThreadPinningObserver = mt_kahypar::parallel::ThreadPinningObserver<HwTopology>;
  using SplittedTBBNumaArena = std::pair<std::unique_ptr<Self>, std::unique_ptr<Self>>;

 public:
  TBBNumaArena(const TBBNumaArena&) = delete;
  TBBNumaArena & operator= (const TBBNumaArena &) = delete;

  TBBNumaArena(TBBNumaArena&&) = delete;
  TBBNumaArena & operator= (TBBNumaArena &&) = delete;

  static TBBNumaArena& instance(const size_t num_threads = 1) {
    static TBBNumaArena instance(num_threads);
    return instance;
  }

  int total_number_of_threads() const {
    return _num_threads;
  }

  int number_of_threads_on_numa_node(const int node) const {
    ASSERT(node < (int)_arenas.size());
    return _arenas[node].max_concurrency();
  }

  int num_used_numa_nodes() const {
    return _arenas.size();
  }

  tbb::task_arena& numa_task_arena(const int node) {
    ASSERT(node < (int)_arenas.size());
    return _arenas[node];
  }

  tbb::task_group& numa_task_group(const int node) {
    ASSERT(_arenas.size() <= _groups.size());
    ASSERT(node < (int)_arenas.size());
    return _groups[node];
  }

  SplittedTBBNumaArena split_tbb_numa_arena(const size_t num_threads_0, const size_t num_threads_1) {
    ASSERT(num_threads_0 + num_threads_1 <= static_cast<size_t>(total_number_of_threads()),
      V(num_threads_0) << V(num_threads_1) << V(total_number_of_threads()));
    const size_t num_numa_nodes = num_used_numa_nodes();
    // At least one open slot should be available on each numa node
    std::vector<std::vector<int>> cpu_to_numa_node_0(num_numa_nodes);
    std::vector<std::vector<int>> cpu_to_numa_node_1(num_numa_nodes);

    auto num_threads_on_numa_node = [&](const size_t node, const size_t num_threads) {
      ASSERT(node < _cpus_to_numa_node.size());
      return static_cast<double>(num_threads) * (
        static_cast<double>(_cpus_to_numa_node[node].size()) /
        static_cast<double>(num_threads_0 + num_threads_1) );
    };


    size_t threads_left_0 = num_threads_0;
    size_t threads_left_1 = num_threads_1;
    for ( size_t node = 0; node < num_numa_nodes; ++node ) {
      size_t num_cpus_on_numa_node_0 = 0;
      size_t num_cpus_on_numa_node_1 = 0;
      if ( threads_left_0 >= threads_left_1 ) {
        num_cpus_on_numa_node_0 = std::max( std::min(
          static_cast<size_t>(std::ceil(num_threads_on_numa_node(node, num_threads_0))),
          threads_left_0 ), 1UL );
        num_cpus_on_numa_node_1 = std::max( std::min(
          static_cast<size_t>(std::floor(num_threads_on_numa_node(node, num_threads_1))),
          threads_left_1 ), 1UL );
      } else {
        num_cpus_on_numa_node_0 = std::max( std::min(
          static_cast<size_t>(std::floor(num_threads_on_numa_node(node, num_threads_0))),
          threads_left_0 ), 1UL );
        num_cpus_on_numa_node_1 = std::max( std::min(
          static_cast<size_t>(std::ceil(num_threads_on_numa_node(node, num_threads_1))),
          threads_left_1 ), 1UL );
      }
      threads_left_0 -= (num_cpus_on_numa_node_0 <= threads_left_0 ? num_cpus_on_numa_node_0 : threads_left_0);
      threads_left_1 -= (num_cpus_on_numa_node_1 <= threads_left_1 ? num_cpus_on_numa_node_1 : threads_left_1);

      std::vector<int>& cpus_on_numa_node = _cpus_to_numa_node[node];
      ASSERT(num_cpus_on_numa_node_0 <= cpus_on_numa_node.size());
      ASSERT(num_cpus_on_numa_node_1 <= cpus_on_numa_node.size());
      for ( size_t i = 0; i < num_cpus_on_numa_node_0; ++i ) {
        cpu_to_numa_node_0[node].push_back(cpus_on_numa_node[i]);
      }
      for ( size_t i = cpus_on_numa_node.size() - num_cpus_on_numa_node_1;
            i < cpus_on_numa_node.size(); ++i ) {
        cpu_to_numa_node_1[node].push_back(cpus_on_numa_node[i]);
      }
    }

    return std::make_pair(
      std::unique_ptr<Self>(new Self(cpu_to_numa_node_0)),
      std::unique_ptr<Self>(new Self(cpu_to_numa_node_1)) );
  }

  template <typename F>
  void execute_sequential_on_all_numa_nodes(F&& func) {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      numa_task_arena(node).execute([&] {
            numa_task_group(node).run([&, node] {
              func(node);
            });
          });
      wait(node, numa_task_group(node));
    }
  }

  template <typename F>
  void execute_parallel_on_all_numa_nodes(F&& func) {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      numa_task_arena(node).execute([&] {
            numa_task_group(node).run([&, node] {
              func(node);
            });
          });
    }
    wait();
  }

  void wait() {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      _arenas[node].execute([&, node] {
            _groups[node].wait();
          });
    }
  }

  void wait(const int node, tbb::task_group& group) {
    ASSERT(node < (int)_arenas.size());
    _arenas[node].execute([&] {
          group.wait();
        });
  }

  void terminate() {
    for (tbb::task_arena& arena : _arenas) {
      arena.terminate();
    }

    for (ThreadPinningObserver& observer : _observer) {
      observer.observe(false);
    }

    if ( _global_observer ) {
      _global_observer->observe(false);
    }

    if ( _init ) {
      _init->terminate();
    }
  }

 private:
  explicit TBBNumaArena(const int num_threads) :
    _num_threads(num_threads),
    _init(std::make_unique<tbb::task_scheduler_init>(num_threads)),
    _global_observer(nullptr),
    _arenas(),
    _groups(HwTopology::instance().num_numa_nodes()),
    _observer(),
    _cpus_to_numa_node() {
    HwTopology& topology = HwTopology::instance();
    int num_numa_nodes = topology.num_numa_nodes();
    DBG << "Initialize TBB with" << num_threads << "threads";
    _arenas.reserve(num_numa_nodes);
    // TODO(heuer): fix copy constructor of observer
    _observer.reserve(num_numa_nodes);

    std::vector<int> cpus = topology.get_all_cpus();
    // Sort cpus in the following order
    // 1.) Non-hyperthread first
    // 2.) Increasing order of numa node
    // 3.) Increasing order of cpu id
    // ...
    std::sort(cpus.begin(), cpus.end(),
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
    while (static_cast<int>(cpus.size()) > _num_threads) {
      cpus.pop_back();
    }
    _global_observer = std::make_unique<ThreadPinningObserver>(cpus);

    _cpus_to_numa_node.resize(num_numa_nodes);
    for ( const int cpu_id : cpus ) {
      int node = topology.numa_node_of_cpu(cpu_id);
      ASSERT(node < static_cast<int>(_cpus_to_numa_node.size()));
      _cpus_to_numa_node[node].push_back(cpu_id);
    }
    while( !_cpus_to_numa_node.empty() && _cpus_to_numa_node.back().empty() ) {
      _cpus_to_numa_node.pop_back();
    }

    initialize_tbb_numa_arenas();
  }

  explicit TBBNumaArena(const std::vector<std::vector<int>>& cpus_to_numa_node)
    : _num_threads(0),
      _init(nullptr),
      _global_observer(nullptr),
      _arenas(),
      _groups(cpus_to_numa_node.size()),
      _observer(),
      _cpus_to_numa_node(cpus_to_numa_node) {
    ASSERT(_cpus_to_numa_node.size() <= static_cast<size_t>(HwTopology::instance().num_numa_nodes()));
    int num_numa_nodes = _cpus_to_numa_node.size();
    _arenas.reserve(num_numa_nodes);
    _observer.reserve(num_numa_nodes);

    for ( size_t node = 0; node < _cpus_to_numa_node.size(); ++node ) {
      _num_threads += _cpus_to_numa_node[node].size();
    }

    initialize_tbb_numa_arenas();
  }

  void initialize_tbb_numa_arenas() {
    for (size_t node = 0; node < _cpus_to_numa_node.size(); ++node) {
      int num_cpus_on_numa_node = _cpus_to_numa_node[node].size();
      ASSERT(num_cpus_on_numa_node <= HwTopology::instance().num_cpus_on_numa_node(node));
      if (num_cpus_on_numa_node > 0) {
        DBG << "Initialize TBB task arena on numa node" << node
            << "with" << num_cpus_on_numa_node << "threads";
        _arenas.emplace_back(num_cpus_on_numa_node, 0);
        _arenas.back().initialize();
        _observer.emplace_back(_arenas.back(), node, _cpus_to_numa_node[node]);
      }
    }
  }

  int _num_threads;
  std::unique_ptr<tbb::task_scheduler_init> _init;
  std::unique_ptr<ThreadPinningObserver> _global_observer;
  std::vector<tbb::task_arena> _arenas;
  std::vector<tbb::task_group> _groups;
  std::vector<ThreadPinningObserver> _observer;
  std::vector<std::vector<int>> _cpus_to_numa_node;
};
}  // namespace parallel
}  // namespace mt_kahypar
