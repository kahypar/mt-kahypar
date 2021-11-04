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

#include <mt-kahypar/datastructures/async/array_lock_manager.h>
#include <tbb/parallel_do.h>
#include <mt-kahypar/partition/refinement/label_propagation/async_lp_refiner.h>
#include <mt-kahypar/parallel/atomic_wrapper.h>
#include <mt-kahypar/utils/progress_bar.h>
#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"

namespace mt_kahypar {

class NLevelCoarsenerBase {
 private:

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = true;

  using ParallelHyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>>;
  using SeedDeduplicator = kahypar::ds::FastResetFlagArray<HypernodeID>;
  using RegionComparator = ds::NodeRegionComparator<Hypergraph>;
  using TreeGroupPool = ds::ConcurrentQueueGroupPool<ds::UncontractionGroupTree, RegionComparator>;
  using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;

 public:
  NLevelCoarsenerBase(Hypergraph& hypergraph,
                      const Context& context,
                      const bool top_level) :
    _is_finalized(false),
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _phg(),
    _compactified_hg(),
    _compactified_phg(),
    _compactified_hn_mapping(),
    _hierarchy(),
    _removed_hyperedges_batches(),
    _round_coarsening_times(),
    _group_pools_for_versions(),
    _lock_manager_for_async() { }

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

  PartitionedHypergraph&& doAsynchronousUncoarsen(std::unique_ptr<IRefiner>& global_fm);

  void
  uncoarsenAsyncTask(TreeGroupPool *pool, metrics::ThreadSafeMetrics &current_metrics,
                     IAsyncRefiner *async_lp_refiner, IAsyncRefiner *async_fm_refiner,
                     HypernodeID &uncontraction_counter,
                     utils::ProgressBar &uncontraction_progress,
                     AsyncNodeTracker &async_node_tracker,
                     RegionComparator &node_region_comparator,
                     SeedDeduplicator &seed_deduplicator, const size_t task_id,
                     const bool alwaysInsertIntoPQ, size_t &local_calls_to_localized_refine,
                     size_t &local_iterations_in_localized_refine,
                     utils::PerTaskTimerForAsync &task_local_timer);

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void uncontractGroupAsyncSubtask(const ds::ContractionGroup &group,
                                                                      const ds::ContractionGroupID groupID);


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

  metrics::ThreadSafeMetrics computeMetricsForAsync(PartitionedHypergraph& phg) {
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
  metrics::ThreadSafeMetrics initializeForAsync(PartitionedHypergraph& current_hg);

  void localizedRefine(PartitionedHypergraph &partitioned_hypergraph,
                       const parallel::scalable_vector <HypernodeID> &refinement_nodes,
                       std::unique_ptr<IRefiner> &label_propagation, std::unique_ptr<IRefiner> &fm,
                       kahypar::Metrics &current_metrics, const bool force_measure_timings,
                       size_t &num_refinement_iterations);

  void localizedRefineForAsync(PartitionedHypergraph &partitioned_hypergraph,
                               const IteratorRange<IAsyncRefiner::NodeIteratorT> &refinement_nodes,
                               IAsyncRefiner *async_lp, IAsyncRefiner *async_fm,
                               ds::ContractionGroupID group_id,
                               metrics::ThreadSafeMetrics &current_metrics, size_t &num_iterations);

  void globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                    std::unique_ptr<IRefiner>& fm,
                    kahypar::Metrics& current_metrics,
                    const double time_limit);

  // ! True, if coarsening terminates and finalize function was called
  bool _is_finalized;

  // ! Original hypergraph
  Hypergraph& _hg;

  const Context& _context;
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

  // ! Contains the contraction group pools that are used in asynchronous (or sequential) uncoarsening. There is one pool
  // ! per hypergraph version which is used to manage the uncontraction hierarchy.
  VersionedPoolVector _group_pools_for_versions;

  // ! A lock manager for locks on hypernodes used in asynchronous n-level uncoarsening
  std::unique_ptr<ds::GroupLockManager> _lock_manager_for_async;

  static constexpr HypernodeID ASYNC_UPDATE_PROGRESS_BAR_THRESHOLD = 25000;

};
}  // namespace mt_kahypar
