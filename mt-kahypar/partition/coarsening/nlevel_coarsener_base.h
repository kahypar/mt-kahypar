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

#include <mt-kahypar/datastructures/asynch/array_lock_manager.h>
#include <tbb/parallel_do.h>
#include <mt-kahypar/partition/refinement/label_propagation/asynch_lp_refiner.h>
#include <mt-kahypar/partition/refinement/thread_local_asynch_refiners.h>
#include <mt-kahypar/parallel/atomic_wrapper.h>
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

  PartitionedHypergraph&& doSequentialUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                                std::unique_ptr<IRefiner>& fm);


    class PoolUncoarseningParallelBody {
    public:
        PoolUncoarseningParallelBody(const ds::UncontractionGroupTree *hierarchy, NLevelCoarsenerBase &base,
                                     metrics::ThreadSafeMetrics &current_metrics, bool &force_measure_timings,
                                     CAtomic<size_t> &total_uncontractions) :
                                              _hierarchy(hierarchy),
                                              _base(base),
                                              _current_metrics(current_metrics),
                                              _force_measure_timings(force_measure_timings),
                                              _total_uncontractions(total_uncontractions) {}

        void operator()(ds::ContractionGroupID groupID, tbb::parallel_do_feeder<ds::ContractionGroupID>& feeder) const;
    private:

        const ds::UncontractionGroupTree* _hierarchy;
        NLevelCoarsenerBase& _base;
        metrics::ThreadSafeMetrics& _current_metrics;
        bool& _force_measure_timings;
        CAtomic<size_t>& _total_uncontractions;
    };

  PartitionedHypergraph&& doAsynchronousUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
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

  metrics::ThreadSafeMetrics computeMetricsForAsynch(PartitionedHypergraph& phg) {
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
  metrics::ThreadSafeMetrics initializeForAsynch(PartitionedHypergraph& current_hg);

  void localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
                       const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                       std::unique_ptr<IRefiner>& label_propagation,
                       std::unique_ptr<IRefiner>& fm,
                       kahypar::Metrics& current_metrics,
                       const bool force_measure_timings);

  void localizedRefineForAsynch(PartitionedHypergraph &partitioned_hypergraph,
                                const parallel::scalable_vector <HypernodeID> &refinement_nodes,
                                IAsynchRefiner *asynch_lp, ds::ContractionGroupID group_id,
                                metrics::ThreadSafeMetrics &current_metrics, const bool force_measure_timings);

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

  // ! Contains the contraction group pools that are used in asynchronous (or sequential) uncoarsening. There is one pool
  // ! per hypergraph version which is used to manage the uncontraction hierarchy.
  ds::VersionedPoolVector _group_pools_for_versions;

  // ! A lock manager for locks on hypernodes used in asynchronous n-level uncoarsening
  std::unique_ptr<ds::GroupLockManager> _lock_manager_for_async;

  // ! Static references to all thread-local refiners created in asynchronous uncoarsening. As the thread-local refiners
  // ! are implicitly static and their lifetime lasts until program exit, it is necessary to manually delete them whenever one
  // ! call to doAsynchronousUncoarsening() finishes by calling reset() on the unique_ptr's. On the next call to
  // ! doAsynchronousUncoarsening() the threads will realize they each need a new Refiner and build one into the
  // ! previously declared unique_ptr. The memory overhead for this in the case of different threads being used for
  // ! different calls to doAsynchronousUncoarsening() should not be large as used Refiners get destroyed, i.e. the only
  // ! overhead is from this vector itself which holds as many pointers as there are threads used in asynchronous uncoarsening.
  static tbb::concurrent_vector<std::unique_ptr<IAsynchRefiner>*> _thread_local_refiner_ptrs;

  static void destroy_thread_local_refiners() {
      for (std::unique_ptr<IAsynchRefiner>* refiner_ptr : _thread_local_refiner_ptrs) {
          refiner_ptr->reset();
      }
  }
  static void register_thread_local_refiner_ptr(std::unique_ptr<IAsynchRefiner>* refiner_ptr) {
      if (std::find(_thread_local_refiner_ptrs.begin(), _thread_local_refiner_ptrs.end(), refiner_ptr) == _thread_local_refiner_ptrs.end()) {
          _thread_local_refiner_ptrs.push_back(refiner_ptr);
      }
  }

};
}  // namespace mt_kahypar
