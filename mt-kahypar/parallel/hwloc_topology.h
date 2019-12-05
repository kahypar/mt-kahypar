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

namespace mt_kahypar {
namespace parallel {
/**
 * Static class responsible for initializing, destroying and
 * calling hwloc library. Calls to hwloc library are outsourced
 * to this class such that hardware topology can be mocked.
 */
class HwlocTopology {
 public:
  static void initialize(hwloc_topology_t& topology) {
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);
  }

  static hwloc_obj_t get_first_numa_node(hwloc_topology_t topology) {
    int numa_depth = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_NUMANODE);
    return hwloc_get_obj_by_depth(topology, numa_depth, 0);
  }

  static std::vector<int> get_cpus_of_numa_node_without_hyperthreads(hwloc_obj_t node) {
    std::vector<int> cpus;

    auto add_cpu_of_core = [&](hwloc_obj_t node) {
                             ASSERT(node->type == HWLOC_OBJ_CORE);
                             std::vector<int> core_cpus;
                             int cpu_id;
                             hwloc_bitmap_foreach_begin(cpu_id, node->cpuset) {
                               core_cpus.emplace_back(cpu_id);
                             }
                             hwloc_bitmap_foreach_end();
                             // Assume that core consists of two processing units (hyperthreads)
                             ASSERT(!core_cpus.empty());
                             cpus.push_back(core_cpus[0]);
                           };
    enumerate_all_core_units(node, add_cpu_of_core);

    return cpus;
  }

  static std::vector<int> get_cpus_of_numa_node_only_hyperthreads(hwloc_obj_t node) {
    std::vector<int> cpus;

    auto add_cpu_of_core = [&](hwloc_obj_t node) {
                             ASSERT(node->type == HWLOC_OBJ_CORE);
                             std::vector<int> core_cpus;
                             int cpu_id;
                             hwloc_bitmap_foreach_begin(cpu_id, node->cpuset) {
                               core_cpus.emplace_back(cpu_id);
                             }
                             hwloc_bitmap_foreach_end();
                             // Assume that core consists of two processing units (hyperthreads)
                             if ( core_cpus.size() >= 2 ) {
                              cpus.push_back(core_cpus[1]);
                             }
                           };
    enumerate_all_core_units(node, add_cpu_of_core);

    return cpus;
  }

  static void destroy_topology(hwloc_topology_t topology) {
    hwloc_topology_destroy(topology);
  }

 private:
  template <class F>
  static void enumerate_all_core_units(hwloc_obj_t node, F& func) {
    if (node->type == HWLOC_OBJ_CORE) {
      func(node);
      return;
    }

    for (size_t i = 0; i < node->arity; ++i) {
      enumerate_all_core_units(node->children[i], func);
    }
  }

  HwlocTopology() { }
};
}  // namespace parallel
}  // namespace mt_kahypar
