/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/refinement/flows/scheduler.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {

  PartitionedHypergraph&& NLevelUncoarsener::doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                                         std::unique_ptr<IRefiner>& fm) {
    ASSERT(_uncoarseningData.is_finalized);
    Metrics current_metrics = initialize(*_uncoarseningData.compactified_phg);
    if (_context.type == ContextType::main) {
      _context.initial_km1 = current_metrics.km1;
    }

    // Project partition from compactified hypergraph to original hypergraph
    _timer.start_timer("initialize_partition", "Initialize Partition");
    *_uncoarseningData.partitioned_hg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
    _uncoarseningData.partitioned_hg->doParallelForAllNodes([&](const HypernodeID hn) {
      ASSERT(static_cast<size_t>(hn) < _uncoarseningData.compactified_hn_mapping.size());
      const HypernodeID compactified_hn = _uncoarseningData.compactified_hn_mapping[hn];
      const PartitionID block_id = _uncoarseningData.compactified_phg->partID(compactified_hn);
      ASSERT(block_id != kInvalidPartition && block_id < _context.partition.k);
      _uncoarseningData.partitioned_hg->setOnlyNodePart(hn, block_id);
    });
    _uncoarseningData.partitioned_hg->initializePartition();

    if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache
        || _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache_on_demand ) {
      _uncoarseningData.partitioned_hg->allocateGainTableIfNecessary();
      if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache ) {
        _uncoarseningData.partitioned_hg->initializeGainCache();
      }
    }


    ASSERT(metrics::objective(*_uncoarseningData.compactified_phg, _context.partition.objective) ==
           metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective),
           V(metrics::objective(*_uncoarseningData.compactified_phg, _context.partition.objective)) <<
           V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));
    ASSERT(metrics::imbalance(*_uncoarseningData.compactified_phg, _context) ==
           metrics::imbalance(*_uncoarseningData.partitioned_hg, _context),
           V(metrics::imbalance(*_uncoarseningData.compactified_phg, _context)) <<
           V(metrics::imbalance(*_uncoarseningData.partitioned_hg, _context)));
    _timer.stop_timer("initialize_partition");

    // Initialize Flow Refinement Scheduler
    std::unique_ptr<IRefiner> flows(nullptr);
    if ( _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
      flows = std::make_unique<FlowRefinementScheduler>(_hg, _context);
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == Objective::km1 ? current_metrics.km1 : current_metrics.cut,
                                              _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += _uncoarseningData.compactified_hg->initialNumNodes();

    // Initialize Refiner
    if ( label_propagation ) {
      label_propagation->initialize(*_uncoarseningData.partitioned_hg);
    }
    if ( fm ) {
      fm->initialize(*_uncoarseningData.partitioned_hg);
    }

    // Perform batch uncontractions
    bool is_timer_disabled = false;
    bool force_measure_timings = _context.partition.measure_detailed_uncontraction_timings && _context.type == ContextType::main;
    if ( _timer.isEnabled() ) {
      _timer.disable();
      is_timer_disabled = true;
    }

    ASSERT(_uncoarseningData.round_coarsening_times.size() == _uncoarseningData.removed_hyperedges_batches.size());
    _uncoarseningData.round_coarsening_times.push_back(_uncoarseningData.round_coarsening_times.size() > 0 ?
                                                       _uncoarseningData.round_coarsening_times.back() : std::numeric_limits<double>::max()); // Sentinel

    size_t num_batches = 0;
    size_t total_batches_size = 0;
    const size_t minimum_required_number_of_border_vertices = std::max(_context.refinement.max_batch_size,
                                                                       _context.shared_memory.num_threads * _context.refinement.min_border_vertices_per_thread);
    ds::StreamingVector<HypernodeID> tmp_refinement_nodes;
    kahypar::ds::FastResetFlagArray<> border_vertices_of_batch(_uncoarseningData.partitioned_hg->initialNumNodes());
    auto do_localized_refinement = [&]() {
      parallel::scalable_vector<HypernodeID> refinement_nodes = tmp_refinement_nodes.copy_parallel();
      tmp_refinement_nodes.clear_parallel();
      border_vertices_of_batch.reset();
      localizedRefine(*_uncoarseningData.partitioned_hg, refinement_nodes, label_propagation,
                      fm, current_metrics, force_measure_timings);
    };

    while ( !_hierarchy.empty() ) {
      BatchVector& batches = _hierarchy.back();

      // Uncontract all batches of a specific version of the hypergraph
      while ( !batches.empty() ) {
        const Batch& batch = batches.back();
        if ( batch.size() > 0 ) {
          HEAVY_REFINEMENT_ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
                                  current_metrics.getMetric(Mode::direct, _context.partition.objective),
                                  V(current_metrics.getMetric(Mode::direct, _context.partition.objective)) <<
                                  V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));
          _timer.start_timer("batch_uncontractions", "Batch Uncontractions", false, force_measure_timings);
          _uncoarseningData.partitioned_hg->uncontract(batch);
          _timer.stop_timer("batch_uncontractions", force_measure_timings);
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(_uncoarseningData.partitioned_hg->checkTrackedPartitionInformation());
          HEAVY_REFINEMENT_ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
                                  current_metrics.getMetric(Mode::direct, _context.partition.objective),
                                  V(current_metrics.getMetric(Mode::direct, _context.partition.objective)) <<
                                  V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));

          _timer.start_timer("collect_border_vertices", "Collect Border Vertices", false, force_measure_timings);
          tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
            const Memento& memento = batch[i];
            if ( !border_vertices_of_batch[memento.u] && _uncoarseningData.partitioned_hg->isBorderNode(memento.u) ) {
              border_vertices_of_batch.set(memento.u, true);
              tmp_refinement_nodes.stream(memento.u);
            }
            if ( !border_vertices_of_batch[memento.v] && _uncoarseningData.partitioned_hg->isBorderNode(memento.v) ) {
              border_vertices_of_batch.set(memento.v, true);
              tmp_refinement_nodes.stream(memento.v);
            }
          });
          _timer.stop_timer("collect_border_vertices", force_measure_timings);

          if ( tmp_refinement_nodes.size() >= minimum_required_number_of_border_vertices ) {
            // Perform localized refinement if we uncontract more
            // than the minimum required number of border vertices
            do_localized_refinement();
          }

          ++num_batches;
          total_batches_size += batch.size();
          // Update Progress Bar
          uncontraction_progress.setObjective(current_metrics.getMetric(
              _context.partition.mode, _context.partition.objective));
          uncontraction_progress += batch.size();
        }
        batches.pop_back();
      }

      if ( tmp_refinement_nodes.size() > 0 ) {
        // Perform localized refinement on remaining border vertices
        do_localized_refinement();
      }

      // Restore single-pin and parallel nets to continue with the next vector of batches
      if ( !_uncoarseningData.removed_hyperedges_batches.empty() ) {
        _timer.start_timer("restore_single_pin_and_parallel_nets", "Restore Single Pin and Parallel Nets", false, force_measure_timings);
        _uncoarseningData.partitioned_hg->restoreSinglePinAndParallelNets(_uncoarseningData.removed_hyperedges_batches.back());
        _uncoarseningData.removed_hyperedges_batches.pop_back();
        _timer.stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
        HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
        HEAVY_REFINEMENT_ASSERT(_uncoarseningData.partitioned_hg->checkTrackedPartitionInformation());

        // Perform refinement on all vertices
        const double time_limit = refinementTimeLimit(_context, _uncoarseningData.round_coarsening_times.back());
        globalRefine(*_uncoarseningData.partitioned_hg, fm, flows, current_metrics, time_limit);
        uncontraction_progress.setObjective(current_metrics.getMetric(
            _context.partition.mode, _context.partition.objective));
        _uncoarseningData.round_coarsening_times.pop_back();
      }
      _hierarchy.pop_back();
    }

    // Top-Level Refinement on all vertices
    const HyperedgeWeight objective_before = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    const double time_limit = refinementTimeLimit(_context, _uncoarseningData.round_coarsening_times.back());
    globalRefine(*_uncoarseningData.partitioned_hg, fm, flows, current_metrics, time_limit);
    _uncoarseningData.round_coarsening_times.pop_back();
    ASSERT(_uncoarseningData.round_coarsening_times.size() == 0);
    const HyperedgeWeight objective_after = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    if ( _context.partition.verbose_output && objective_after < objective_before ) {
      LOG << GREEN << "Top-Level Refinment improved objective from"
      << objective_before << "to" << objective_after << END;
    }

    if ( is_timer_disabled ) {
      _timer.enable();
    }

    // If we finish batch uncontractions and partition is imbalanced, we try to rebalance it
    if ( _context.type == ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        Mode::direct, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
        << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ") ->"
        << "Rebalancer is activated" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(*_uncoarseningData.partitioned_hg, _context);
      }

      _timer.start_timer("rebalance", "Rebalance");
      if ( _context.partition.objective == Objective::km1 ) {
        Km1Rebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
        rebalancer.rebalance(current_metrics);
      } else if ( _context.partition.objective == Objective::cut ) {
        CutRebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
        rebalancer.rebalance(current_metrics);
      }
      _timer.stop_timer("rebalance");

      const HyperedgeWeight quality_after = current_metrics.getMetric(
        Mode::direct, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        const HyperedgeWeight quality_delta = quality_after - quality_before;
        if ( quality_delta > 0 ) {
          LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
          << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
        } else {
          LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
          << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
        }
      }
    }

    double avg_batch_size = static_cast<double>(total_batches_size) / num_batches;
    utils::Utilities::instance().getStats(_context.utility_id).add_stat("num_batches", static_cast<int64_t>(num_batches));
    utils::Utilities::instance().getStats(_context.utility_id).add_stat("avg_batch_size", avg_batch_size);
    DBG << V(num_batches) << V(avg_batch_size);

    ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
           current_metrics.getMetric(Mode::direct, _context.partition.objective),
           V(current_metrics.getMetric(Mode::direct, _context.partition.objective)) <<
           V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));

    return std::move(*_uncoarseningData.partitioned_hg);
  }
  void NLevelUncoarsener::localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
                                          const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                                          std::unique_ptr<IRefiner>& label_propagation,
                                          std::unique_ptr<IRefiner>& fm,
                                          Metrics& current_metrics,
                                          const bool force_measure_timings) {
    if ( debug && _context.type == ContextType::main ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
      << ", imbalance = " << current_metrics.imbalance;
    }

    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation &&
          _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        _timer.start_timer("label_propagation", "Label Propagation", false, force_measure_timings);
        improvement_found |= label_propagation->refine(partitioned_hypergraph,
                                                       refinement_nodes, current_metrics, std::numeric_limits<double>::max());
        _timer.stop_timer("label_propagation", force_measure_timings);
      }

      if ( fm &&
          _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        _timer.start_timer("fm", "FM", false, force_measure_timings);
        improvement_found |= fm->refine(partitioned_hypergraph,
                                        refinement_nodes, current_metrics, std::numeric_limits<double>::max());
        _timer.stop_timer("fm", force_measure_timings);
      }

      if ( _context.type == ContextType::main ) {
        ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      if ( !_context.refinement.refine_until_no_improvement ) {
        break;
      }
    }

    if ( _context.type == ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
  }

  void NLevelUncoarsener::globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                                       std::unique_ptr<IRefiner>& fm,
                                       std::unique_ptr<IRefiner>& flows,
                                       Metrics& current_metrics,
                                       const double time_limit) {

    auto applyGlobalFMParameters = [&](const FMParameters& fm, const NLevelGlobalFMParameters global_fm){
      NLevelGlobalFMParameters tmp_global_fm;
      tmp_global_fm.num_seed_nodes = fm.num_seed_nodes;
      tmp_global_fm.obey_minimal_parallelism = fm.obey_minimal_parallelism;
      fm.num_seed_nodes = global_fm.num_seed_nodes;
      fm.obey_minimal_parallelism = global_fm.obey_minimal_parallelism;
      return tmp_global_fm;
    };

    if ( _context.refinement.global_fm.use_global_fm ) {
      if ( debug && _context.type == ContextType::main ) {
        io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
        DBG << "Start Refinement - km1 = " << current_metrics.km1
        << ", imbalance = " << current_metrics.imbalance;
      }

      // Enable Timings
      bool was_enabled = false;
      if ( !_timer.isEnabled() &&
           _context.type == ContextType::main ) {
        _timer.enable();
        was_enabled = true;
      }

      // Apply global FM parameters to FM context and temporary store old fm context
      _timer.start_timer("global_refinement", "Global Refinement");
      NLevelGlobalFMParameters tmp_global_fm = applyGlobalFMParameters(
        _context.refinement.fm, _context.refinement.global_fm);
      bool improvement_found = true;
      while( improvement_found ) {
        improvement_found = false;
        const HyperedgeWeight metric_before = current_metrics.getMetric(
          Mode::direct, _context.partition.objective);

        if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
          _timer.start_timer("fm", "FM");
          improvement_found |= fm->refine(partitioned_hypergraph, {}, current_metrics, time_limit);
          _timer.stop_timer("fm");
        }

        if ( flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
          _timer.start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
          flows->initialize(partitioned_hypergraph);
          _timer.stop_timer("initialize_flow_scheduler");

          _timer.start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
          improvement_found |= flows->refine(partitioned_hypergraph, {}, current_metrics, time_limit);
          _timer.stop_timer("flow_refinement_scheduler");
        }

        if ( _context.type == ContextType::main ) {
          ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
                 "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                 << "does not match the metric updated by the refiners" << V(current_metrics.km1));
        }

        const HyperedgeWeight metric_after = current_metrics.getMetric(
          Mode::direct, _context.partition.objective);
        const double relative_improvement = 1.0 -
          static_cast<double>(metric_after) / metric_before;
        if ( !_context.refinement.global_fm.refine_until_no_improvement ||
            relative_improvement <= _context.refinement.relative_improvement_threshold ) {
          break;
        }
      }
      // Reset FM context
      applyGlobalFMParameters(_context.refinement.fm, tmp_global_fm);
      _timer.stop_timer("global_refinement");

      if ( was_enabled ) {
        _timer.disable();
      }

      if ( _context.type == ContextType::main) {
        DBG << "--------------------------------------------------\n";
      }
    }
  }

}
