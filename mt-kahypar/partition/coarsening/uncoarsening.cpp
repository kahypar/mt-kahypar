#include "mt-kahypar/partition/coarsening/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener_base.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/partitioning_output.h"

namespace mt_kahypar {

  void MultilevelCoarsenerBase::finalize() {
    utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
    // Free memory of temporary contraction buffer and
    // release coarsening memory in memory pool
    currentHypergraph().freeTmpContractionBuffer();
    if (_top_level) {
      parallel::MemoryPool::instance().release_mem_group("Coarsening");
    }

    // Construct top level partitioned hypergraph (memory is taken from memory pool)
    _partitioned_hg = PartitionedHypergraph(
            _context.partition.k, _task_group_id, _hg);

    // Construct partitioned hypergraphs parallel
    tbb::task_group group;
    // Construct partitioned hypergraph for each coarsened hypergraph in the hierarchy
    for (size_t i = 0; i < _hierarchy.size(); ++i) {
      group.run([&, i] {
        _hierarchy[i].contractedPartitionedHypergraph() = PartitionedHypergraph(
                _context.partition.k, _task_group_id, _hierarchy[i].contractedHypergraph());
      });
    }
    group.wait();

    // Set the representative partitioned hypergraph for each hypergraph
    // in the hierarchy
    if (_hierarchy.size() > 0) {
      _hierarchy[0].setRepresentativeHypergraph(&_partitioned_hg);
      for (size_t i = 1; i < _hierarchy.size(); ++i) {
        _hierarchy[i].setRepresentativeHypergraph(&_hierarchy[i - 1].contractedPartitionedHypergraph());
      }
    }
    _is_finalized = true;
    utils::Timer::instance().stop_timer("finalize_multilevel_hierarchy");
  }

  void NLevelCoarsenerBase::finalize() {
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

  namespace {
    double refinementTimeLimit(const Context& context, const double time) {
      if ( context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
        const double time_limit_factor = std::max(1.0,  context.refinement.fm.time_limit_factor * context.partition.k);
        return std::max(5.0, time_limit_factor * time);
      } else {
        return std::numeric_limits<double>::max();
      }
    }
  } // namespace

  void MultilevelCoarsenerBase::performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities,
          const HighResClockTimepoint& round_start) {
    ASSERT(!_is_finalized);
    Hypergraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities, _task_group_id);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

  PartitionedHypergraph&& MultilevelCoarsenerBase::doUncoarsen(
          std::unique_ptr<IRefiner>& label_propagation,
          std::unique_ptr<IRefiner>& fm) {
    PartitionedHypergraph& coarsest_hg = currentPartitionedHypergraph();
    kahypar::Metrics current_metrics = initialize(coarsest_hg);

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == kahypar::Objective::km1
                                              ? current_metrics.km1 : current_metrics.cut,
                                              _context.partition.verbose_output &&
                                              _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += coarsest_hg.initialNumNodes();

    // Refine Coarsest Partitioned Hypergraph
    double time_limit = refinementTimeLimit(_context, _hierarchy.back().coarseningTime());
    refine(coarsest_hg, label_propagation, fm, current_metrics, time_limit);

    for (int i = _hierarchy.size() - 1; i >= 0; --i) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      PartitionedHypergraph& representative_hg = _hierarchy[i].representativeHypergraph();
      PartitionedHypergraph& contracted_hg = _hierarchy[i].contractedPartitionedHypergraph();
      representative_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = _hierarchy[i].mapToContractedHypergraph(hn);
        const PartitionID block = contracted_hg.partID(coarse_hn);
        ASSERT(block != kInvalidPartition && block < representative_hg.k());
        representative_hg.setOnlyNodePart(hn, block);
      });
      representative_hg.initializePartition(_task_group_id);

      ASSERT(metrics::objective(representative_hg, _context.partition.objective) ==
             metrics::objective(contracted_hg, _context.partition.objective),
             V(metrics::objective(representative_hg, _context.partition.objective)) <<
                                                                                    V(metrics::objective(
                                                                                            contracted_hg,
                                                                                            _context.partition.objective)));
      ASSERT(metrics::imbalance(representative_hg, _context) ==
             metrics::imbalance(contracted_hg, _context),
             V(metrics::imbalance(representative_hg, _context)) <<
                                                                V(metrics::imbalance(contracted_hg, _context)));
      utils::Timer::instance().stop_timer("projecting_partition");

      // Refinement
      time_limit = refinementTimeLimit(_context, _hierarchy[i].coarseningTime());
      refine(representative_hg, label_propagation, fm, current_metrics, time_limit);

      // Update Progress Bar
      uncontraction_progress.setObjective(
              current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective));
      uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if (_top_level && !metrics::isBalanced(_partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
              kahypar::Mode::direct_kway, _context.partition.objective);
      if (_context.partition.verbose_output) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
            << metrics::imbalance(_partitioned_hg, _context) << ") ->"
            << "Rebalancer is activated" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(_partitioned_hg, _context);
      }

      utils::Timer::instance().start_timer("rebalance", "Rebalance");
      if (_context.partition.objective == kahypar::Objective::km1) {
        Km1Rebalancer rebalancer(_partitioned_hg, _context);
        rebalancer.rebalance(current_metrics);
      } else if (_context.partition.objective == kahypar::Objective::cut) {
        CutRebalancer rebalancer(_partitioned_hg, _context);
        rebalancer.rebalance(current_metrics);
      }
      utils::Timer::instance().stop_timer("rebalance");

      const HyperedgeWeight quality_after = current_metrics.getMetric(
              kahypar::Mode::direct_kway, _context.partition.objective);
      if (_context.partition.verbose_output) {
        const HyperedgeWeight quality_delta = quality_after - quality_before;
        if (quality_delta > 0) {
          LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
              << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
        } else {
          LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
              << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
        }
      }
    }

    ASSERT(metrics::objective(_partitioned_hg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective))
           << V(metrics::objective(_partitioned_hg, _context.partition.objective)));
    return std::move(_partitioned_hg);
  }

  PartitionedHypergraph&& NLevelCoarsenerBase::doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
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

    if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache ) {
      _phg.initializeGainCache();
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

    ASSERT(_round_coarsening_times.size() == _removed_hyperedges_batches.size());
    _round_coarsening_times.push_back(_round_coarsening_times.size() > 0 ?
      _round_coarsening_times.back() : std::numeric_limits<double>::max()); // Sentinel

    size_t num_batches = 0;
    size_t total_batches_size = 0;
    const size_t minimum_required_number_of_border_vertices = std::max(_context.refinement.max_batch_size,
      _context.shared_memory.num_threads * _context.refinement.min_border_vertices_per_thread);
    ds::StreamingVector<HypernodeID> tmp_refinement_nodes;
    kahypar::ds::FastResetFlagArray<> border_vertices_of_batch(_phg.initialNumNodes());
    auto do_localized_refinement = [&]() {
      parallel::scalable_vector<HypernodeID> refinement_nodes = tmp_refinement_nodes.copy_parallel();
      tmp_refinement_nodes.clear_parallel();
      border_vertices_of_batch.reset();
      localizedRefine(_phg, refinement_nodes, label_propagation,
        fm, current_metrics, force_measure_timings);
    };

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
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());
          HEAVY_REFINEMENT_ASSERT(metrics::objective(_phg, _context.partition.objective) ==
                current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                V(metrics::objective(_phg, _context.partition.objective)));

          utils::Timer::instance().start_timer("collect_border_vertices", "Collect Border Vertices", false, force_measure_timings);
          tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
            const Memento& memento = batch[i];
            if ( !border_vertices_of_batch[memento.u] && _phg.isBorderNode(memento.u) ) {
              border_vertices_of_batch.set(memento.u, true);
              tmp_refinement_nodes.stream(memento.u);
            }
            if ( !border_vertices_of_batch[memento.v] && _phg.isBorderNode(memento.v) ) {
              border_vertices_of_batch.set(memento.v, true);
              tmp_refinement_nodes.stream(memento.v);
            }
          });
          utils::Timer::instance().stop_timer("collect_border_vertices", force_measure_timings);

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
      if ( !_removed_hyperedges_batches.empty() ) {
        utils::Timer::instance().start_timer("restore_single_pin_and_parallel_nets", "Restore Single Pin and Parallel Nets", false, force_measure_timings);
        _phg.restoreSinglePinAndParallelNets(_removed_hyperedges_batches.back());
        _removed_hyperedges_batches.pop_back();
        utils::Timer::instance().stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
        HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
        HEAVY_REFINEMENT_ASSERT(_phg.checkTrackedPartitionInformation());

        // Perform refinement on all vertices
        const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
        globalRefine(_phg, fm, current_metrics, time_limit);
        uncontraction_progress.setObjective(current_metrics.getMetric(
          _context.partition.mode, _context.partition.objective));
        _round_coarsening_times.pop_back();
      }
      _hierarchy.pop_back();
    }

    // Top-Level Refinement on all vertices
    const HyperedgeWeight objective_before = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    const double time_limit = refinementTimeLimit(_context, _round_coarsening_times.back());
    globalRefine(_phg, fm, current_metrics, time_limit);
    _round_coarsening_times.pop_back();
    ASSERT(_round_coarsening_times.size() == 0);
    const HyperedgeWeight objective_after = current_metrics.getMetric(
      _context.partition.mode, _context.partition.objective);
    if ( _context.partition.verbose_output && objective_after < objective_before ) {
      LOG << GREEN << "Top-Level Refinment improved objective from"
          << objective_before << "to" << objective_after << END;
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


  void MultilevelCoarsenerBase::refine(
          PartitionedHypergraph& partitioned_hypergraph,
          std::unique_ptr<IRefiner>& label_propagation,
          std::unique_ptr<IRefiner>& fm,
          kahypar::Metrics& current_metrics,
          const double time_limit) {

    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        label_propagation->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_lp_refiner");

        utils::Timer::instance().start_timer("label_propagation", "Label Propagation");
        improvement_found |= label_propagation->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("label_propagation");
      }

      if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        fm->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_fm_refiner");

        utils::Timer::instance().start_timer("fm", "FM");
        improvement_found |= fm->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("fm");
      }

      if ( _top_level ) {
        ASSERT(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)
               == metrics::objective(partitioned_hypergraph, _context.partition.objective),
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

  void NLevelCoarsenerBase::localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
                                            const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                                            std::unique_ptr<IRefiner>& label_propagation,
                                            std::unique_ptr<IRefiner>& fm,
                                            kahypar::Metrics& current_metrics,
                                            const bool force_measure_timings) {
    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation &&
           _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("label_propagation", "Label Propagation", false, force_measure_timings);
        improvement_found |= label_propagation->refine(partitioned_hypergraph,
          refinement_nodes, current_metrics, std::numeric_limits<double>::max());
        utils::Timer::instance().stop_timer("label_propagation", force_measure_timings);
      }

      if ( fm &&
           _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("fm", "FM", false, force_measure_timings);
        improvement_found |= fm->refine(partitioned_hypergraph,
          refinement_nodes, current_metrics, std::numeric_limits<double>::max());
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

  namespace {
    NLevelGlobalFMParameters applyGlobalFMParameters(const FMParameters& fm,
                                                     const NLevelGlobalFMParameters global_fm) {
      NLevelGlobalFMParameters tmp_global_fm;
      tmp_global_fm.num_seed_nodes = fm.num_seed_nodes;
      tmp_global_fm.obey_minimal_parallelism = fm.obey_minimal_parallelism;
      fm.num_seed_nodes = global_fm.num_seed_nodes;
      fm.obey_minimal_parallelism = global_fm.obey_minimal_parallelism;
      return tmp_global_fm;
    }
  } // namespace

  void NLevelCoarsenerBase::globalRefine(PartitionedHypergraph& partitioned_hypergraph,
                                         std::unique_ptr<IRefiner>& fm,
                                         kahypar::Metrics& current_metrics,
                                         const double time_limit) {
    if ( _context.refinement.global_fm.use_global_fm ) {
      if ( debug && _top_level ) {
        io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
        DBG << "Start Refinement - km1 = " << current_metrics.km1
            << ", imbalance = " << current_metrics.imbalance;
      }

      // Apply global FM parameters to FM context and temporary store old fm context
      NLevelGlobalFMParameters tmp_global_fm = applyGlobalFMParameters(
        _context.refinement.fm, _context.refinement.global_fm);
      bool improvement_found = true;
      while( improvement_found ) {
        improvement_found = false;

        if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
          utils::Timer::instance().start_timer("global_fm", "Global FM", false, _top_level);
          improvement_found |= fm->refine(partitioned_hypergraph, {}, current_metrics, time_limit);
          utils::Timer::instance().stop_timer("global_fm", _top_level);
        }

        if ( _top_level ) {
          ASSERT(current_metrics.km1 == metrics::km1(partitioned_hypergraph),
                "Actual metric" << V(metrics::km1(partitioned_hypergraph))
                                << "does not match the metric updated by the refiners" << V(current_metrics.km1));
        }

        if ( !_context.refinement.global_fm.refine_until_no_improvement ) {
          break;
        }
      }
      // Reset FM context
      applyGlobalFMParameters(_context.refinement.fm, tmp_global_fm);

      if ( _top_level) {
        DBG << "--------------------------------------------------\n";
      }
    }
  }

  kahypar::Metrics MultilevelCoarsenerBase::initialize(PartitionedHypergraph& phg) {
    kahypar::Metrics m = { 0, 0, 0.0 };
    tbb::parallel_invoke([&] {
      m.cut = metrics::hyperedgeCut(phg);
    }, [&] {
      m.km1 = metrics::km1(phg);
    });
    m.imbalance = metrics::imbalance(phg, _context);

    int64_t num_nodes = phg.initialNumNodes();
    int64_t num_edges = phg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", m.cut);
    utils::Stats::instance().add_stat("initial_km1", m.km1);
    utils::Stats::instance().add_stat("initial_imbalance", m.imbalance);
    return m;
  }

  kahypar::Metrics NLevelCoarsenerBase::initialize(PartitionedHypergraph& current_hg) {
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

}