/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/parallel/thread_management.h"

#ifndef KAHYPAR_DISABLE_HWLOC
  #include <hwloc.h>
#endif

#include "mt-kahypar/parallel/tbb_initializer.h"
#ifndef KAHYPAR_DISABLE_HWLOC
  #include "mt-kahypar/parallel/hardware_topology.h"
#endif

namespace mt_kahypar {

#ifndef KAHYPAR_DISABLE_HWLOC
  using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
  using TBBInitializer = mt_kahypar::parallel::TBBInitializer<HardwareTopology, false>;
#else
  using TBBInitializer = mt_kahypar::parallel::SimpleTBBInitializer;
#endif

static_assert(parallel::provides_hardware_information == TBBInitializer::provides_numa_information);

namespace parallel {

void initialize_tbb(const size_t num_threads) {
  mt_kahypar::TBBInitializer::initialize(num_threads);
}

void terminate_tbb() {
  mt_kahypar::TBBInitializer::instance().terminate();
}

int total_number_of_threads() {
  return mt_kahypar::TBBInitializer::instance().total_number_of_threads();
}

#ifndef KAHYPAR_DISABLE_HWLOC
size_t num_hardware_cpus() {
  return HardwareTopology<>::instance().num_cpus();
}

int num_used_numa_nodes() {
  return mt_kahypar::TBBInitializer::instance().num_used_numa_nodes();
}

void activate_interleaved_membind_policy() {
  hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
  parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
  hwloc_bitmap_free(cpuset);
}
#endif

}  // namespace parallel
}  // namespace mt_kahypar
