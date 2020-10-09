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
                      const TaskGroupID task_group_id,
                      const bool top_level) :
    _is_finalized(false),
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _top_level(top_level),
    _phg(),
    _compactified_hg(),
    _compactified_phg(),
    _compactified_hn_mapping(),
    _hierarchy(),
    _removed_hyperedges_batches(),
    _round_coarsening_times() { }

  NLevelCoarsenerBase(const NLevelCoarsenerBase&) = delete;
  NLevelCoarsenerBase(NLevelCoarsenerBase&&) = delete;
  NLevelCoarsenerBase & operator= (const NLevelCoarsenerBase &) = delete;
  NLevelCoarsenerBase & operator= (NLevelCoarsenerBase &&) = delete;

  virtual ~NLevelCoarsenerBase() = default;

 protected:

  Hypergraph& compactifiedHypergraph() {
    ASSERT(_is_finalized);
    return _compactified_hg;
  }

  PartitionedHypergraph& compactifiedPartitionedHypergraph() {
    ASSERT(_is_finalized);
    return _compactified_phg;
  }

  void finalize();

  void removeSinglePinAndParallelNets(const HighResClockTimepoint& round_start) {
    utils::Timer::instance().start_timer("remove_single_pin_and_parallel_nets", "Remove Single Pin and Parallel Nets");
    _removed_hyperedges_batches.emplace_back(_hg.removeSinglePinAndParallelHyperedges());
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _round_coarsening_times.push_back(elapsed_time);
    utils::Timer::instance().stop_timer("remove_single_pin_and_parallel_nets");
  }

  PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm);

 protected:
  kahypar::Metrics computeMetrics(PartitionedHypergraph& phg) {
    HyperedgeWeight cut = 0;
    HyperedgeWeight km1 = 0;
    tbb::parallel_invoke([&] {
      cut = metrics::hyperedgeCut(phg);
    }, [&] {
      km1 = metrics::km1(phg);
    });
    return { cut, km1,  metrics::imbalance(phg, _context) };
  }

  kahypar::Metrics initialize(PartitionedHypergraph& current_hg);

  void localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
                       const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                       std::unique_ptr<IRefiner>& label_propagation,
                       std::unique_ptr<IRefiner>& fm,
                       kahypar::Metrics& current_metrics,
                       const bool force_measure_timings);

  void globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                    std::unique_ptr<IRefiner>& fm,
                    kahypar::Metrics& current_metrics,
                    const double time_limit);

  // ! True, if coarsening terminates and finalize function was called
  bool _is_finalized;

  // ! Original hypergraph
  Hypergraph& _hg;

  const Context& _context;
  const TaskGroupID _task_group_id;
  const bool _top_level;

  // ! Original partitioned hypergraph
  PartitionedHypergraph _phg;
  // ! Once coarsening terminates we generate a compactified hypergraph
  // ! containing only enabled vertices and hyperedges within a consecutive
  // ! ID range, which is then used for initial partitioning
  Hypergraph _compactified_hg;
  // ! Compactified partitioned hypergraph
  PartitionedHypergraph _compactified_phg;
  // ! Mapping from vertex IDs of the original hypergraph to the IDs
  // ! in the compactified hypergraph
  parallel::scalable_vector<HypernodeID> _compactified_hn_mapping;

  // ! Represents the n-level hierarchy
  // ! A batch is vector of uncontractions/mementos that can be uncontracted in parallel
  // ! without conflicts. All batches of a specific version of the hypergraph are assembled
  // ! in a batch vector. Each time we perform single-pin and parallel net detection we create
  // ! a new version (simply increment a counter) of the hypergraph. Once a batch vector is
  // ! completly processed single-pin and parallel nets have to be restored.
  VersionedBatchVector _hierarchy;
  // ! Removed single-pin and parallel nets.
  // ! All hyperedges that are contained in one vector must be restored once
  // ! we completly processed a vector of batches.
  ParallelHyperedgeVector _removed_hyperedges_batches;
  // ! Contains timings how long a coarsening pass takes for each round
  parallel::scalable_vector<double> _round_coarsening_times;
};
}  // namespace mt_kahypar
