/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

class NLevelCoarsenerBase {
 private:

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using ParallelHyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>>;

 public:
  NLevelCoarsenerBase(Hypergraph& hypergraph,
                      const Context& context,
                      const bool top_level) :
    _is_finalized(false),
    _hg(hypergraph),
    _context(context),
    _top_level(top_level) { }

  NLevelCoarsenerBase(const NLevelCoarsenerBase&) = delete;
  NLevelCoarsenerBase(NLevelCoarsenerBase&&) = delete;
  NLevelCoarsenerBase & operator= (const NLevelCoarsenerBase &) = delete;
  NLevelCoarsenerBase & operator= (NLevelCoarsenerBase &&) = delete;

  virtual ~NLevelCoarsenerBase() = default;

 protected:

  Hypergraph& compactifiedHypergraph() {
    ASSERT(_is_finalized);
    return *_compactified_hg;
  }

  PartitionedHypergraph& compactifiedPartitionedHypergraph() {
    ASSERT(_is_finalized);
    return *_compactified_phg;
  }

  void finalize();

  void removeSinglePinAndParallelNets(const HighResClockTimepoint& round_start) {
    utils::Timer::instance().start_timer("remove_single_pin_and_parallel_nets", "Remove Single Pin and Parallel Nets");
    _removed_hyperedges_batches->emplace_back(_hg.removeSinglePinAndParallelHyperedges());
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _round_coarsening_times->push_back(elapsed_time);
    utils::Timer::instance().stop_timer("remove_single_pin_and_parallel_nets");
  }

 protected:
  // ! True, if coarsening terminates and finalize function was called
  bool _is_finalized;

  // ! Original hypergraph
  Hypergraph& _hg;

  const Context& _context;
  const bool _top_level;

  // ! Original partitioned hypergraph
  std::shared_ptr<PartitionedHypergraph> _phg;
  // ! Once coarsening terminates we generate a compactified hypergraph
  // ! containing only enabled vertices and hyperedges within a consecutive
  // ! ID range, which is then used for initial partitioning
  std::shared_ptr<Hypergraph> _compactified_hg;
  // ! Compactified partitioned hypergraph
  std::shared_ptr<PartitionedHypergraph> _compactified_phg;
  // ! Mapping from vertex IDs of the original hypergraph to the IDs
  // ! in the compactified hypergraph
  std::shared_ptr<vec<HypernodeID>> _compactified_hn_mapping;

  // ! Represents the n-level hierarchy
  // ! A batch is vector of uncontractions/mementos that can be uncontracted in parallel
  // ! without conflicts. All batches of a specific version of the hypergraph are assembled
  // ! in a batch vector. Each time we perform single-pin and parallel net detection we create
  // ! a new version (simply increment a counter) of the hypergraph. Once a batch vector is
  // ! completly processed single-pin and parallel nets have to be restored.
  std::shared_ptr<VersionedBatchVector> _hierarchy;
  // ! Removed single-pin and parallel nets.
  // ! All hyperedges that are contained in one vector must be restored once
  // ! we completly processed a vector of batches.
  std::shared_ptr<ParallelHyperedgeVector> _removed_hyperedges_batches;
  // ! Contains timings how long a coarsening pass takes for each round
  std::shared_ptr<vec<double>> _round_coarsening_times;
};
}  // namespace mt_kahypar
