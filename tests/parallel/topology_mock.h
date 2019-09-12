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

#include <vector>
#include <algorithm>
#include <numeric>
#include <thread>
#include <hwloc.h>

#include "kahypar/macros.h"


namespace mt_kahypar {
namespace parallel {

using numa_to_cpu_t = std::vector<std::vector<int>>;

struct Cpu {
  int cpu_id;
  bool is_hyperthread;
};

struct Node {

  explicit Node(const int node) :
    os_index(node),
    cpuset(),
    next_cousin(),
    cpus() { 
    cpuset = hwloc_bitmap_alloc();
  }

  Node(const Node& other) :
    os_index(other.os_index),
    cpuset(hwloc_bitmap_dup(other.cpuset)),
    next_cousin(other.next_cousin),
    cpus(other.cpus) { }

  Node& operator= (const Node& other) {
    os_index = other.os_index;
    cpuset = hwloc_bitmap_dup(other.cpuset);
    next_cousin = other.next_cousin;
    cpus = other.cpus;
    return *this;
  };

  ~Node() {
    hwloc_bitmap_free(cpuset);
  }

  void init_cpuset(const std::vector<int>& cores,
                   const std::vector<int>& hyperthreads) {
    for ( const int cpu : cores ) {
      hwloc_bitmap_set(cpuset, cpu);
      cpus.emplace_back(Cpu { cpu, false });
    }

    for ( const int cpu : hyperthreads ) {
      hwloc_bitmap_set(cpuset, cpu);
      cpus.emplace_back(Cpu { cpu, true });
    }
  }

  int os_index;
  hwloc_cpuset_t cpuset;
  Node* next_cousin;
  std::vector<Cpu> cpus;
};

class Topology {
 public:
  explicit Topology(const numa_to_cpu_t& numa_to_cpu) :
    _nodes() { 
    init(numa_to_cpu);
  }

  Node* get_node(const int node) {
    ASSERT(node < (int) _nodes.size());
    return &_nodes[node];
  }

 private:
  void init(const numa_to_cpu_t& numa_to_cpu) {
    ASSERT(numa_to_cpu.size() % 2 == 0);
    int numa_nodes = numa_to_cpu.size() / 2;
    for ( int node = 0; node < numa_nodes; ++node ) {
      _nodes.emplace_back(node);
      _nodes.back().init_cpuset(numa_to_cpu[node], numa_to_cpu[node + numa_nodes]);
    }
    for ( int node = 0; node < numa_nodes - 1; ++node ) {
      _nodes[node].next_cousin = &_nodes[node + 1];
    }
  }

  std::vector<Node> _nodes;
};

static numa_to_cpu_t split_physical_cpus_into_numa_nodes(const int num_numa_nodes) {
  numa_to_cpu_t numa_to_cpu;

  std::vector<int> cpus(std::thread::hardware_concurrency());
  std::iota(cpus.begin(), cpus.end(), 0);
  ASSERT(num_numa_nodes <= (int) cpus.size());
  int cpus_per_node = cpus.size() / num_numa_nodes + (cpus.size() % num_numa_nodes > 0);
  std::reverse(cpus.begin(), cpus.end());
  for ( int node = 0; node < num_numa_nodes; ++node ) {
    int current_cpus = 0;
    numa_to_cpu.emplace_back();
    while ( current_cpus < cpus_per_node && !cpus.empty() ) {
      numa_to_cpu.back().push_back(cpus.back());
      cpus.pop_back();
      current_cpus++;
    }
  }

  return numa_to_cpu;
}

using node_t = Node*;
using topology_t = Topology*;

template< int NUM_NUMA_NODES >
class TopologyMock {

 public:
  static void initialize(topology_t& topology) {
    // We split into 2 * NUM_NUMA_NODES in order to simulate hyperthreading
    topology = new Topology(split_physical_cpus_into_numa_nodes(2 * NUM_NUMA_NODES));
  }

  static node_t get_first_numa_node(topology_t topology) {
    return topology->get_node(0);
  }

  static std::vector<int> get_cpus_of_numa_node_without_hyperthreads(node_t node) {
    std::vector<int> cpus;
    for ( const Cpu& cpu : node->cpus ) {
      if ( !cpu.is_hyperthread ) {
        cpus.push_back(cpu.cpu_id);
      }
    }
    return cpus;
  }

  static std::vector<int> get_cpus_of_numa_node_only_hyperthreads(node_t node) {
    std::vector<int> cpus;
    for ( const Cpu& cpu : node->cpus ) {
      if ( cpu.is_hyperthread ) {
        cpus.push_back(cpu.cpu_id);
      }
    }
    return cpus;
  }

  static void destroy_topology(topology_t) { }

 private:
  TopologyMock() { }
};

} // namespace parallel
} // namespace mt_kahypar