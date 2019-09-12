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

  static void destroy_topology(hwloc_topology_t topology) {
    hwloc_topology_destroy(topology);
  }

 private:
  HwlocTopology() { }
};

} // namespace parallel
} // namespace mt_kahypar