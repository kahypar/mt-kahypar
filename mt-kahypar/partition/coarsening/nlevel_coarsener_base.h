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
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/stats.h"

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
    _removed_hyperedges_batches() { }

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

  void finalize() {
    // Create compactified hypergraph containing only enabled vertices and hyperedges
    // with consecutive IDs => Less complexity in initial partitioning.
    utils::Timer::instance().start_timer("compactify_hypergraph", "Compactify Hypergraph");
    auto compactification = HypergraphFactory::compactify(_task_group_id, _hg);
    _compactified_hg = std::move(compactification.first);
    _compactified_hn_mapping = std::move(compactification.second);
    _compactified_phg = PartitionedHypergraph(_context.partition.k, _task_group_id, _compactified_hg);
    utils::Timer::instance().stop_timer("compactify_hypergraph");

    // Create n-level batch uncontraction hierarchy
    utils::Timer::instance().start_timer("create_batch_uncontraction_hierarchy", "Create n-Level Hierarchy");
    _hierarchy = _hg.createBatchUncontractionHierarchy(_context.refinement.max_batch_size);
    ASSERT(_removed_hyperedges_batches.size() == _hierarchy.size() - 1);
    utils::Timer::instance().stop_timer("create_batch_uncontraction_hierarchy");

    _is_finalized = true;
  }

  void removeSinglePinAndParallelNets() {
    utils::Timer::instance().start_timer("remove_single_pin_and_parallel_nets", "Remove Single Pin and Parallel Nets");
    _removed_hyperedges_batches.emplace_back(_hg.removeSinglePinAndParallelHyperedges());
    utils::Timer::instance().stop_timer("remove_single_pin_and_parallel_nets");
  }

  PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm) {
    ASSERT(_is_finalized);
    kahypar::Metrics current_metrics = initialize(_compactified_phg);

    // Project partition from compactified hypergraph to original hypergraph
    utils::Timer::instance().start_timer("initialize_partition", "Initialize Partition");
    _phg = PartitionedHypergraph(_context.partition.k, _task_group_id, _hg);
    _phg.doParallelForAllNodes([&](const HypernodeID hn) {
      ASSERT(static_cast<size_t>(hn) < _compactified_hn_mapping.size());
      const HypernodeID compactified_hn = _compactified_hn_mapping[hn];
      const PartitionID block_id = _compactified_phg.partID(compactified_hn);
      ASSERT(block_id != kInvalidPartition && block_id < _context.partition.k);
      _phg.setOnlyNodePart(hn, block_id);
    });
    _phg.initializePartition(_task_group_id);
    if ( _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
      _phg.initializeGainInformation();
    }

    ASSERT(metrics::objective(_compactified_phg, _context.partition.objective) ==
            metrics::objective(_phg, _context.partition.objective),
            V(metrics::objective(_compactified_phg, _context.partition.objective)) <<
            V(metrics::objective(_phg, _context.partition.objective)));
    ASSERT(metrics::imbalance(_compactified_phg, _context) ==
            metrics::imbalance(_phg, _context),
            V(metrics::imbalance(_compactified_phg, _context)) <<
            V(metrics::imbalance(_phg, _context)));
    utils::Timer::instance().stop_timer("initialize_partition");

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
      _context.partition.objective == kahypar::Objective::km1 ? current_metrics.km1 : current_metrics.cut,
      _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += _compactified_hg.initialNumNodes();

    // Initialize Refiner
    if ( label_propagation ) {
      label_propagation->initialize(_phg);
    }
    if ( fm ) {
      fm->initialize(_phg);
    }

    // Perform batch uncontractions
    bool is_timer_disabled = false;
    bool force_measure_timings = _context.partition.measure_detailed_uncontraction_timings && _top_level;
    if ( utils::Timer::instance().isEnabled() ) {
      utils::Timer::instance().disable();
      is_timer_disabled = true;
    }

    size_t num_batches = 0;
    size_t total_batches_size = 0;
    while ( !_hierarchy.empty() ) {
      BatchVector& batches = _hierarchy.back();

      // Uncontract all batches of a specific version of the hypergraph
      while ( !batches.empty() ) {
        const Batch& batch = batches.back();
        if ( batch.size() > 0 ) {
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                V(metrics::objective(_phg, _context.partition.objective)));
          utils::Timer::instance().start_timer("batch_uncontractions", "Batch Uncontractions", false, force_measure_timings);
          _phg.uncontract(batch);
          utils::Timer::instance().stop_timer("batch_uncontractions", force_measure_timings);
          HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                V(metrics::objective(_phg, _context.partition.objective)));

          // Perform refinement
          refine(_phg, batch, label_propagation, fm, current_metrics, force_measure_timings);

          ++num_batches;
          total_batches_size += batch.size();
          // Update Progress Bar
          uncontraction_progress.setObjective(current_metrics.getMetric(
            _context.partition.mode, _context.partition.objective));
          uncontraction_progress += batch.size();
        }
        batches.pop_back();
      }

      // Restore single-pin and parallel nets to continue with the next vector of batches
      if ( !_removed_hyperedges_batches.empty() ) {
        utils::Timer::instance().start_timer("restore_single_pin_and_parallel_nets", "Restore Single Pin and Parallel Nets", false, force_measure_timings);
        _phg.restoreSinglePinAndParallelNets(_removed_hyperedges_batches.back());
        _removed_hyperedges_batches.pop_back();
        utils::Timer::instance().stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
      }
      _hierarchy.pop_back();
    }

    if ( is_timer_disabled ) {
      utils::Timer::instance().enable();
    }

    // If we finish batch uncontractions and partition is imbalanced, we try to rebalance it
    if ( _top_level && !metrics::isBalanced(_phg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
            << metrics::imbalance(_phg, _context) << ") ->"
            << "Rebalancer is activated" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(_phg, _context);
      }

      utils::Timer::instance().start_timer("rebalance", "Rebalance");
      if ( _context.partition.objective == kahypar::Objective::km1 ) {
        Km1Rebalancer rebalancer(_phg, _context);
        rebalancer.rebalance(current_metrics);
      } else if ( _context.partition.objective == kahypar::Objective::cut ) {
        CutRebalancer rebalancer(_phg, _context);
        rebalancer.rebalance(current_metrics);
      }
      utils::Timer::instance().stop_timer("rebalance");

      const HyperedgeWeight quality_after = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        const HyperedgeWeight quality_delta = quality_after - quality_before;
        if ( quality_delta > 0 ) {
          LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
              << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
        } else {
          LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
              << "(Current Imbalance:" << metrics::imbalance(_phg, _context) << ")" << END;
        }
      }
    }

    double avg_batch_size = static_cast<double>(total_batches_size) / num_batches;
    utils::Stats::instance().add_stat("num_batches", static_cast<int64_t>(num_batches));
    utils::Stats::instance().add_stat("avg_batch_size", avg_batch_size);
    DBG << V(num_batches) << V(avg_batch_size);

    ASSERT(metrics::objective(_phg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_phg, _context.partition.objective)));

    return std::move(_phg);
  }

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

  kahypar::Metrics initialize(PartitionedHypergraph& current_hg) {
    kahypar::Metrics current_metrics = computeMetrics(current_hg);
    int64_t num_nodes = current_hg.initialNumNodes();
    int64_t num_edges = current_hg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);
    return current_metrics;
  }

  void refine(PartitionedHypergraph& partitioned_hypergraph,
              const Batch& batch,
              std::unique_ptr<IRefiner>& label_propagation,
              std::unique_ptr<IRefiner>& fm,
              kahypar::Metrics& current_metrics,
              const bool force_measure_timings) {
    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph, "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation &&
           _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("label_propagation", "Label Propagation", false, force_measure_timings);
        improvement_found |= label_propagation->refine(partitioned_hypergraph, batch, current_metrics, std::numeric_limits<double>::max());
        utils::Timer::instance().stop_timer("label_propagation", force_measure_timings);
      }

      if ( fm &&
           _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("fm", "FM", false, force_measure_timings);
        improvement_found |= fm->refine(partitioned_hypergraph, batch, current_metrics, std::numeric_limits<double>::max());
        utils::Timer::instance().stop_timer("fm", force_measure_timings);
      }

      if ( _top_level ) {
        ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      if ( !_context.refinement.refine_until_no_improvement ) {
        break;
      }
    }

    if ( _top_level) {
      DBG << "--------------------------------------------------\n";
    }
  }

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
};
}  // namespace mt_kahypar
