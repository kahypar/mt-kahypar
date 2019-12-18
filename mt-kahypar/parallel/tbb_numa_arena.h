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
#include "mt-kahypar/parallel/global_thread_pinning_observer.h"
#include "mt-kahypar/parallel/numa_thread_pinning_observer.h"

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
  using GlobalThreadPinning = mt_kahypar::parallel::GlobalThreadPinning<HwTopology>;
  using GlobalThreadPinningObserver = mt_kahypar::parallel::GlobalThreadPinningObserver<HwTopology>;
  using NumaThreadPinningObserver = mt_kahypar::parallel::NumaThreadPinningObserver<HwTopology>;
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
    ASSERT(num_threads_0 + num_threads_1 == static_cast<size_t>(total_number_of_threads()));
    const size_t num_numa_nodes = num_used_numa_nodes();
    // At least one open slot should be available on each numa node
    std::vector<int> used_cpus_on_numa_nodes_0(num_numa_nodes, 0);
    std::vector<int> used_cpus_on_numa_nodes_1(num_numa_nodes, 0);

    auto splitted_number_of_threads_on_numa_node = [&](const size_t node, const size_t threads) {
      return static_cast<double>(threads) * ( static_cast<double>(number_of_threads_on_numa_node(node)) /
        static_cast<double>(num_threads_0 + num_threads_1) );
    };

    int threads_left_0 = num_threads_0;
    int threads_left_1 = num_threads_1;
    for ( size_t node = 0; node < num_numa_nodes; ++node ) {
      double threads_on_numa_node_0 = splitted_number_of_threads_on_numa_node(node, num_threads_0);
      double threads_on_numa_node_1 = splitted_number_of_threads_on_numa_node(node, num_threads_1);
      if ( threads_left_0 >= threads_left_1 ) {
        used_cpus_on_numa_nodes_0[node] = std::max( std::min(
          static_cast<size_t>(std::ceil(threads_on_numa_node_0)),
          static_cast<size_t>(threads_left_0) ), 1UL );
        used_cpus_on_numa_nodes_1[node] = std::max( std::min(
          static_cast<size_t>(std::floor(threads_on_numa_node_1)),
          static_cast<size_t>(threads_left_1) ), 1UL );
      } else {
        used_cpus_on_numa_nodes_0[node] = std::max( std::min(
          static_cast<size_t>(std::floor(threads_on_numa_node_0)),
          static_cast<size_t>(threads_left_0) ), 1UL );
        used_cpus_on_numa_nodes_1[node] = std::max( std::min(
          static_cast<size_t>(std::ceil(threads_on_numa_node_1)),
          static_cast<size_t>(threads_left_1) ), 1UL );
      }
      threads_left_0 -= (used_cpus_on_numa_nodes_0[node] <= threads_left_0 ?
        used_cpus_on_numa_nodes_0[node] : threads_left_0);
      threads_left_1 -= (used_cpus_on_numa_nodes_1[node] <= threads_left_1 ?
        used_cpus_on_numa_nodes_1[node] : threads_left_1);
    }

    return std::make_pair(
      std::unique_ptr<Self>(new Self(used_cpus_on_numa_nodes_0)),
      std::unique_ptr<Self>(new Self(used_cpus_on_numa_nodes_1)) );
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
    for (NumaThreadPinningObserver& observer : _observer) {
      observer.observe(false);
    }

    for (tbb::task_arena& arena : _arenas) {
      arena.terminate();
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
    _global_observer(std::make_unique<GlobalThreadPinningObserver>()),
    _arenas(),
    _groups(HwTopology::instance().num_numa_nodes()),
    _observer() {
    HwTopology& topology = HwTopology::instance();
    int threads_left = num_threads;
    int num_numa_nodes = topology.num_numa_nodes();
    DBG << "Initialize TBB with" << num_threads << "threads";
    _arenas.reserve(num_numa_nodes);
    // TODO(heuer): fix copy constructor of observer
    _observer.reserve(num_numa_nodes);

    std::vector<int> used_cpus_on_numa_node(num_numa_nodes);
    // First use cores (to prevent using hyperthreads)
    for (int node = 0; node < num_numa_nodes && threads_left > 0; ++node) {
      used_cpus_on_numa_node[node] = std::min(threads_left, topology.num_cores_on_numa_node(node));
      threads_left -= used_cpus_on_numa_node[node];
    }

    // If there are still thread to assign left we use hyperthreading
    for (int node = 0; node < num_numa_nodes && threads_left > 0; ++node) {
      int num_hyperthreads = std::min(threads_left, topology.num_cpus_on_numa_node(node) - used_cpus_on_numa_node[node]);
      used_cpus_on_numa_node[node] += num_hyperthreads;
      threads_left -= num_hyperthreads;
    }
    _num_threads -= threads_left;

    initialize_tbb_numa_arenas(used_cpus_on_numa_node, false);

    // Initialize Global Thread Pinning
    GlobalThreadPinning::instance(num_threads);
    for (int node = 0; node < num_numa_nodes; ++node) {
      // Seems that there is one extra worker threads when num_threads is equal to,
      // but only one will participate in task scheduling at a time
      int num_cpus = std::max(used_cpus_on_numa_node[node], 2);
      topology.use_only_num_cpus_on_numa_node(node, num_cpus);
    }
    _global_observer->observe(true);
  }

  explicit TBBNumaArena(const std::vector<int>& used_cpus_on_numa_node)
    : _num_threads(0),
      _init(nullptr),
      _global_observer(nullptr),
      _arenas(),
      _groups(used_cpus_on_numa_node.size()),
      _observer() {
    ASSERT(used_cpus_on_numa_node.size() <= static_cast<size_t>(HwTopology::instance().num_numa_nodes()));
    int num_numa_nodes = used_cpus_on_numa_node.size();
    _arenas.reserve(num_numa_nodes);
    _observer.reserve(num_numa_nodes);

    for ( size_t node = 0; node < used_cpus_on_numa_node.size(); ++node ) {
      _num_threads += used_cpus_on_numa_node[node];
    }

    initialize_tbb_numa_arenas(used_cpus_on_numa_node, true, false);
  }

  void initialize_tbb_numa_arenas(const std::vector<int>& used_cpus_on_numa_node,
                                  const bool oversubscription_allowed,
                                  const bool reserve_for_master = true) {
    for (size_t node = 0; node < used_cpus_on_numa_node.size(); ++node) {
      ASSERT(used_cpus_on_numa_node[node] <= HwTopology::instance().num_cpus_on_numa_node(node));
      if (used_cpus_on_numa_node[node] > 0) {
        int num_cpus = used_cpus_on_numa_node[node];
        DBG << "Initialize TBB task arena on numa node" << node
            << "with" << num_cpus << "threads";
        _arenas.emplace_back(num_cpus, num_cpus == 1 && reserve_for_master ? 1 : 0);
        _observer.emplace_back(_arenas.back(), node, oversubscription_allowed);
      }
    }
  }

  int _num_threads;
  std::unique_ptr<tbb::task_scheduler_init> _init;
  std::unique_ptr<GlobalThreadPinningObserver> _global_observer;
  std::vector<tbb::task_arena> _arenas;
  std::vector<tbb::task_group> _groups;
  std::vector<NumaThreadPinningObserver> _observer;
};
}  // namespace parallel
}  // namespace mt_kahypar
