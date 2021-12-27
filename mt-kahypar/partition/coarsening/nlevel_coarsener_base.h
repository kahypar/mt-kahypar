/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"
#include <mt-kahypar/partition/coarsening/coarsening_commons.h>

namespace mt_kahypar {

class NLevelCoarsenerBase {
 private:

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using ParallelHyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>>;

 public:
  NLevelCoarsenerBase(Hypergraph& hypergraph,
                      const Context& context,
                      UncoarseningData& uncoarseningData) :
    _hg(hypergraph),
    _context(context),
    _uncoarseningData(uncoarseningData) { }

  NLevelCoarsenerBase(const NLevelCoarsenerBase&) = delete;
  NLevelCoarsenerBase(NLevelCoarsenerBase&&) = delete;
  NLevelCoarsenerBase & operator= (const NLevelCoarsenerBase &) = delete;
  NLevelCoarsenerBase & operator= (NLevelCoarsenerBase &&) = delete;

  virtual ~NLevelCoarsenerBase() = default;

 protected:

  Hypergraph& compactifiedHypergraph() {
    ASSERT(_uncoarseningData.is_finalized);
    return *_uncoarseningData.compactified_hg;
  }

  PartitionedHypergraph& compactifiedPartitionedHypergraph() {
    ASSERT(_uncoarseningData.is_finalized);
    return *_uncoarseningData.compactified_phg;
  }

  void removeSinglePinAndParallelNets(const HighResClockTimepoint& round_start) {
    utils::Timer::instance().start_timer("remove_single_pin_and_parallel_nets", "Remove Single Pin and Parallel Nets");
    _uncoarseningData.removed_hyperedges_batches.emplace_back(_hg.removeSinglePinAndParallelHyperedges());
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _uncoarseningData.round_coarsening_times.push_back(elapsed_time);
    utils::Timer::instance().stop_timer("remove_single_pin_and_parallel_nets");
  }

 protected:
  // ! Original hypergraph
  Hypergraph& _hg;

  const Context& _context;

  UncoarseningData& _uncoarseningData;
};
}  // namespace mt_kahypar
