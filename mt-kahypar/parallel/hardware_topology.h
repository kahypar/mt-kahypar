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

#include <hwloc.h>
#include <mutex>
#include <thread>
#include <vector>

#include "mt-kahypar/macros.h"

#include "mt-kahypar/parallel/global_thread_pinning.h"
#include "mt-kahypar/parallel/hwloc_topology.h"

namespace mt_kahypar {
namespace parallel {
/**
 * Class represents the hardware topology of the system.
 * Internally it uses hwloc library to find numa nodes and corresponding
 * cpus. Furthermore, it implements functionalities to pin logical threads
 * to cpus of a specific numa node.
 *
 * Template parameters can be replaced in order to mock hardware topology and
 * simulate a NUMA on a UMA system.
 */
template <typename HwTopology = HwlocTopology,
          typename Topology = hwloc_topology_t,
          typename Node = hwloc_obj_t>
class HardwareTopology {
 private:
  static constexpr bool debug = false;

  using Self = HardwareTopology<HwTopology, Topology, Node>;
  using GlobalThreadPinning = mt_kahypar::parallel::GlobalThreadPinning<Self>;

  struct Cpu {
    int cpu_id;
    bool is_hyperthread;
  };

  class NumaNode {
   public:
    NumaNode(Node node) :
      _node_id(node->os_index),
      _cpuset(node->cpuset),
      _free_cpus(0),
      _cpus(),
      _mutex(),
      _num_cores(0) {
      for (const int cpu_id :  HwTopology::get_cpus_of_numa_node_without_hyperthreads(node)) {
        _cpus.emplace_back(Cpu { cpu_id, false });
        _num_cores++;
      }
      for (const int cpu_id :  HwTopology::get_cpus_of_numa_node_only_hyperthreads(node)) {
        _cpus.emplace_back(Cpu { cpu_id, true });
      }
      _free_cpus = _cpus.size();
    }

    NumaNode(const NumaNode&) = delete;
    NumaNode & operator= (const NumaNode &) = delete;

    NumaNode(NumaNode&& other) :
      _node_id(other._node_id),
      _cpuset(std::move(other._cpuset)),
      _free_cpus(other._free_cpus),
      _cpus(std::move(other._cpus)),
      _mutex(),
      _num_cores(other._num_cores) { }

    int get_id() const {
      return _node_id;
    }

    hwloc_cpuset_t get_cpuset() const {
      return _cpuset;
    }

    std::vector<int> cpus() {
      std::lock_guard<std::mutex> lock(_mutex);
      std::vector<int> cpus;
      for (const Cpu& cpu : _cpus) {
        cpus.push_back(cpu.cpu_id);
      }
      return cpus;
    }

    size_t num_cores_on_numa_node() const {
      return _num_cores;
    }

    size_t num_cpus_on_numa_node() const {
      return _cpus.size();
    }

    bool is_hyperthread(const int cpu_id) {
      std::lock_guard<std::mutex> lock(_mutex);
      size_t pos = 0;
      for ( ; pos < _cpus.size(); ++pos) {
        if (_cpus[pos].cpu_id == cpu_id) {
          break;
        }
      }
      ASSERT(pos < _cpus.size(), "CPU" << cpu_id << "not found on numa node" << _node_id);
      return _cpus[pos].is_hyperthread;
    }

    void use_only_num_cpus(const size_t num_cpus) {
      std::lock_guard<std::mutex> lock(_mutex);
      // Sort cpus such that hyperthreads are at the end of the list
      // such that they are removed first
      std::sort(_cpus.begin(), _cpus.end(), [&](const Cpu& lhs, const Cpu& rhs) {
            return lhs.is_hyperthread < rhs.is_hyperthread ||
            (lhs.is_hyperthread == rhs.is_hyperthread && lhs.cpu_id < rhs.cpu_id);
          });
      while (_cpus.size() > num_cpus) {
        _cpus.pop_back();
      }
      _free_cpus = _cpus.size();
    }

    int pin_thread_to_cpu() {
      std::lock_guard<std::mutex> lock(_mutex);
      int cpu_id = -1;
      if ( _free_cpus > 0 ) {
        Cpu& cpu = _cpus.front();
        cpu_id = cpu.cpu_id;
        std::swap(_cpus[0], _cpus[--_free_cpus]);
      }
      return cpu_id;
    }

    void unpin_thread_from_cpu(int cpu_id) {
      std::lock_guard<std::mutex> lock(_mutex);
      size_t pos = _free_cpus;
      // Find corresponding cpu
      while (pos < _cpus.size()) {
        if (_cpus[pos].cpu_id == cpu_id) {
          break;
        }
        ++pos;
      }
      ASSERT(pos < _cpus.size(), "CPU" << cpu_id << "not found on numa node"
        << _node_id << "( TID =" << std::this_thread::get_id() << ")");
      std::swap(_cpus[_free_cpus++], _cpus[pos]);
    }

   private:
    int _node_id;
    hwloc_cpuset_t _cpuset;
    size_t _free_cpus;
    std::vector<Cpu> _cpus;
    std::mutex _mutex;
    size_t _num_cores;
  };

 public:
  HardwareTopology(const HardwareTopology&) = delete;
  HardwareTopology & operator= (const HardwareTopology &) = delete;

  HardwareTopology(HardwareTopology&&) = delete;
  HardwareTopology & operator= (HardwareTopology &&) = delete;

  ~HardwareTopology() {
    HwTopology::destroy_topology(_topology);
  }

  static HardwareTopology& instance() {
    static HardwareTopology instance;
    return instance;
  }

  size_t num_numa_nodes() const {
    return _numa_nodes.size();
  }

  size_t num_cpus() const {
    return _num_cpus;
  }

  int numa_node_of_cpu(const int cpu_id) const {
    ASSERT(cpu_id < (int)_cpu_to_numa_node.size());
    ASSERT(_cpu_to_numa_node[cpu_id] != std::numeric_limits<int>::max());
    return _cpu_to_numa_node[cpu_id];
  }

  bool is_hyperthread(const int cpu_id) {
    int node = numa_node_of_cpu(cpu_id);
    return _numa_nodes[node].is_hyperthread(cpu_id);
  }

  // ! Restrict the number of cpus on numa node
  void use_only_num_cpus_on_numa_node(const int node, const size_t num_cpus) {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    _numa_nodes[node].use_only_num_cpus(num_cpus);
  }

  // ! Number of Cores on NUMA node
  int num_cores_on_numa_node(const int node) const {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].num_cores_on_numa_node();
  }

  // ! Number of CPUs on NUMA node
  int num_cpus_on_numa_node(const int node) const {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].num_cpus_on_numa_node();
  }

  // ! CPU bitmap of NUMA node
  hwloc_cpuset_t get_cpuset_of_numa_node(int node) const {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].get_cpuset();
  }

  // ! List of CPUs of NUMA node
  std::vector<int> get_cpus_of_numa_node(int node) {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].cpus();
  }

  // ! List of all available CPUs
  std::vector<int> get_all_cpus() {
    std::vector<int> cpus;
    for ( size_t node = 0; node < num_numa_nodes(); ++node ) {
      for ( const int cpu_id : _numa_nodes[node].cpus() ) {
        cpus.push_back(cpu_id);
      }
    }
    return cpus;
  }

  // ! Pins a thread to a NUMA node
  bool pin_thread_to_numa_node(const int node) {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    int cpu_id = _numa_nodes[node].pin_thread_to_cpu();
    if ( cpu_id != -1 ) {
      GlobalThreadPinning::instance().pin_thread_to_numa_node(node, cpu_id);
    }
    return cpu_id != -1;
  }

  // ! Unpin a thread from a NUMA node
  void unpin_thread_from_numa_node(const int node) {
    ASSERT(node < (int)_numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    int cpu_id = sched_getcpu();
    _numa_nodes[node].unpin_thread_from_cpu(cpu_id);
    GlobalThreadPinning::instance().unpin_thread_from_numa_node(node);
  }

 private:
  HardwareTopology() :
    _num_cpus(0),
    _topology(),
    _numa_nodes(),
    _cpu_to_numa_node(std::thread::hardware_concurrency(),
      std::numeric_limits<int>::max()) {
    HwTopology::initialize(_topology);
    init_numa_nodes();
  }

  void init_numa_nodes() {
    Node node = HwTopology::get_first_numa_node(_topology);
    while (node != nullptr) {
      _numa_nodes.emplace_back(node);
      node = node->next_cousin;
      for (const int cpu_id : _numa_nodes.back().cpus()) {
        ASSERT(cpu_id < (int)_cpu_to_numa_node.size());
        _cpu_to_numa_node[cpu_id] = _numa_nodes.back().get_id();
        ++_num_cpus;
      }
    }
  }

  static std::mutex _mutex;

  size_t _num_cpus;
  Topology _topology;
  std::vector<NumaNode> _numa_nodes;
  std::vector<int> _cpu_to_numa_node;
};

template <typename HwTopology, typename Topology, typename Node>
std::mutex HardwareTopology<HwTopology, Topology, Node>::_mutex;
}  // namespace parallel
}  // namespace mt_kahypar
