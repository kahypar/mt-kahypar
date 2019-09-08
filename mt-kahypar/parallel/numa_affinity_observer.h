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

#undef __TBB_ARENA_OBSERVER
#define __TBB_ARENA_OBSERVER true
#include "tbb/task_scheduler_observer.h"
#undef __TBB_ARENA_OBSERVER

#include "kahypar/macros.h"

namespace kahypar {
namespace parallel {

template< typename HwTopology >
class NumaAffinityObserver : public tbb::task_scheduler_observer {
  using Base = tbb::task_scheduler_observer;

 public:
  explicit NumaAffinityObserver(tbb::task_arena& arena,
                                int numa_node) :
    Base(arena),
    _arena(arena),
    _topology(HwTopology::instance()),
    _numa_node(numa_node) {
    observe(true);
  }

  NumaAffinityObserver(const NumaAffinityObserver&) = delete;
  NumaAffinityObserver& operator= (const NumaAffinityObserver&) = delete;

  NumaAffinityObserver(NumaAffinityObserver&& other) = default;
  NumaAffinityObserver& operator= (NumaAffinityObserver&&) = delete;

  ~NumaAffinityObserver() {
    observe(false);
  }

  void on_scheduler_entry(bool) override {
    _topology.set_affinity_to_numa_node(_numa_node);
  }

  void on_scheduler_exit(bool) override {
    _topology.free_thread_from_numa_node(_numa_node);
  }

 private:
  tbb::task_arena& _arena;
  HwTopology& _topology;
  int _numa_node;
};

} // namespace parallel
} // namespace kahypar