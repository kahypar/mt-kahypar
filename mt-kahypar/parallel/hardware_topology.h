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
#include <vector>
#include <mutex>
#include <thread>

#include "kahypar/macros.h"

#include "mt-kahypar/parallel/hwloc_topology.h"

namespace kahypar {
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
template < typename HwTopology = HwlocTopology,
           typename Topology = hwloc_topology_t,
           typename Node = hwloc_obj_t >
class HardwareTopology {

 private:
  static constexpr bool debug = false;

  struct Cpu {
    int cpu_id;
    size_t num_assigned_threads;
  };

  class NumaNode {

   public:
    NumaNode(Node node) :
      _node_id(node->os_index),
      _cpuset(node->cpuset),
      _cpus(),
      _mutex() { 
		  // assign CPUs from this node's CPU-set
      int cpu_id;
      hwloc_bitmap_foreach_begin(cpu_id, node->cpuset) {
        _cpus.emplace_back(Cpu { cpu_id, 0 });
      }
      hwloc_bitmap_foreach_end();
    }

    NumaNode(const NumaNode&) = delete;
    NumaNode& operator= (const NumaNode&) = delete;

    NumaNode(NumaNode&& other) :
     _node_id(other._node_id),
     _cpuset(std::move(other._cpuset)),
     _cpus(std::move(other._cpus)),
     _mutex() { }

    NumaNode& operator= (NumaNode&&) = default;

    int get_id() const {
      return _node_id;
    }

    hwloc_cpuset_t get_cpuset() const {
      return _cpuset;
    }

    // ! List of CPUs of NUMA node (only testing)
    std::vector<int> cpus() {
      std::lock_guard<std::mutex> lock(_mutex);
      std::vector<int> cpus;
      for ( const Cpu& cpu : _cpus ) {
        cpus.push_back(cpu.cpu_id);
      }
      return cpus;
    }

    size_t num_cpus_on_numa_node() const {
      return _cpus.size();
    }

    int pin_thread_to_cpu() {
      std::lock_guard<std::mutex> lock(_mutex);
      Cpu& cpu = _cpus.front();
      int cpu_id = cpu.cpu_id;
      cpu.num_assigned_threads++;
      size_t pos = 0;
      // Keep cpus sorted in increasing order of their number of assigned logical threads,
      // such that a thread is always assigned to the cpu with least number of logical threads.
      while ( pos < _cpus.size() - 1 && 
              _cpus[pos].num_assigned_threads > _cpus[pos + 1].num_assigned_threads ) {
        std::swap(_cpus[pos], _cpus[pos + 1]);
        ++pos;
      }
      return cpu_id;
    }

    void unpin_thread_from_cpu(int cpu_id) {
      std::lock_guard<std::mutex> lock(_mutex);
      size_t pos = 0;
      // Find corresponding cpu
      while ( pos < _cpus.size() ) {
        if ( _cpus[pos].cpu_id == cpu_id ) {
          break;
        }
        ++pos;
      }
      ASSERT(pos != _cpus.size(), "CPU" << cpu_id << "not found on numa node" << _node_id);
      ASSERT(_cpus[pos].num_assigned_threads > 0, "No thread assigned to cpu" << cpu_id);
      _cpus[pos].num_assigned_threads--;
      // Keep cpus sorted in increasing order of their number of assigned logical threads,
      // such that a thread is always assigned to the cpu with least number of logical threads.
      while ( pos > 0 && _cpus[pos - 1].num_assigned_threads > _cpus[pos].num_assigned_threads) {
        std::swap(_cpus[pos - 1], _cpus[pos]);
        --pos;
      }
    }

   private:
    int _node_id;
    hwloc_cpuset_t _cpuset;
    std::vector<Cpu> _cpus;
    std::mutex _mutex;
  };

 public:
  HardwareTopology(const HardwareTopology&) = delete;
  HardwareTopology& operator= (const HardwareTopology&) = delete;

  HardwareTopology(HardwareTopology&&) = delete;
  HardwareTopology& operator= (HardwareTopology&&) = delete;

  ~HardwareTopology() {
    HwTopology::destroy_topology(_topology);
  }

  static HardwareTopology& instance() {
    if ( _instance == nullptr ) {
      std::lock_guard<std::mutex> _lock(_mutex);
      if ( _instance == nullptr ) {
        _instance = new HardwareTopology();
      }
    }
    return *_instance;
  }

  size_t num_numa_nodes() const {
    return _numa_nodes.size();
  }

  // ! Number of CPUs on NUMA node
  int num_cpus_on_numa_node(const int node) const {
    ASSERT(node < (int) _numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].num_cpus_on_numa_node();
  }

  // ! CPU bitmap of NUMA node
  hwloc_cpuset_t get_cpuset_of_numa_node(int node) const {
    ASSERT(node < (int) _numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].get_cpuset();
  }

  // ! List of CPUs of NUMA node (only testing)
  std::vector<int> get_cpus_of_numa_node(int node) {
    ASSERT(node < (int) _numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    return _numa_nodes[node].cpus();
  }

  // ! Pins a thread to a NUMA node
  void pin_thread_to_numa_node(const int node) {
    ASSERT(node < (int) _numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
		const size_t size = CPU_ALLOC_SIZE( _num_cpus );
    int cpu_id = _numa_nodes[node].pin_thread_to_cpu();
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu_id, &mask);
    const int err = sched_setaffinity(0, size, &mask);

		if ( err ) {
			LOG << "Failed to set thread affinity on numa node" << node;
			exit( EXIT_FAILURE );
		}
    DBG << "Assigned thread with PID" << std::this_thread::get_id()
        << "to cpu" << cpu_id << "on numa node" << node;
  }

  // ! Unpin a thread from a NUMA node
  void unpin_thread_from_numa_node(const int node)  {
    ASSERT(node < (int) _numa_nodes.size());
    ASSERT(_numa_nodes[node].get_id() == node);
    int cpu_id = sched_getcpu();
    _numa_nodes[node].unpin_thread_from_cpu(cpu_id);
    DBG << "Free thread with PID" << std::this_thread::get_id()
        << "on cpu" << cpu_id << "from numa node" << node;
  }

 private:
   HardwareTopology() :
    _num_cpus(std::thread::hardware_concurrency()),
    _topology(),
    _numa_nodes() { 
    HwTopology::initialize(_topology);
    init_numa_nodes();
  }

  void init_numa_nodes() {
    Node node = HwTopology::get_first_numa_node(_topology);
    while ( node != nullptr ) {
      _numa_nodes.emplace_back(node);
      node = node->next_cousin;
    }
  }

  static std::mutex _mutex;
  static HardwareTopology* _instance;

  const size_t _num_cpus;
  Topology _topology;
  std::vector<NumaNode> _numa_nodes;
};

template < typename HwTopology, typename Topology, typename Node >
HardwareTopology<HwTopology, Topology, Node>* HardwareTopology<HwTopology, Topology, Node>::_instance { nullptr };
template < typename HwTopology, typename Topology, typename Node >
std::mutex HardwareTopology<HwTopology, Topology, Node>::_mutex;

} // namespace parallel
} // namespace kahypar