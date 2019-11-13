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

#include <algorithm>
#include <thread>

#include "tbb/task_scheduler_observer.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/global_thread_pinning.h"

namespace mt_kahypar {
namespace parallel {


template< typename HwTopology >
class GlobalThreadPinningObserver : public tbb::task_scheduler_observer {
  using Base = tbb::task_scheduler_observer;

  static constexpr bool debug = false;

  using GlobalThreadPinning = mt_kahypar::parallel::GlobalThreadPinning<HwTopology>;

 public:
  explicit GlobalThreadPinningObserver() :
    Base(true) { }

  GlobalThreadPinningObserver(const GlobalThreadPinningObserver&) = delete;
  GlobalThreadPinningObserver& operator= (const GlobalThreadPinningObserver&) = delete;

  GlobalThreadPinningObserver& operator= (GlobalThreadPinningObserver&&) = delete;

  ~GlobalThreadPinningObserver() {
    observe(false);
  }

  void on_scheduler_entry(bool) override {
    GlobalThreadPinning::instance().pin_thread();
  }

  void on_scheduler_exit(bool) override {
    GlobalThreadPinning::instance().unregister_thread();
  }
};

} // namespace parallel
} // namespace mt_kahypar