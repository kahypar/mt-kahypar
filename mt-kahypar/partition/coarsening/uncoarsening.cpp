#include "multilevel_coarsener_base.h"

#include "mt-kahypar/parallel/memory_pool.h"
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

  void MultilevelCoarsenerBase::performMultilevelContraction(parallel::scalable_vector<HypernodeID>&& communities,
                                    const HighResClockTimepoint& round_start) {
    ASSERT(!_is_finalized);
    Hypergraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities, _task_group_id);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    _hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

  PartitionedHypergraph&& MultilevelCoarsenerBase::doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
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
    double time_limit = refinementTimeLimit(_hierarchy.back());
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
      time_limit = refinementTimeLimit(_hierarchy[i]);
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


  void MultilevelCoarsenerBase::refine(PartitionedHypergraph& partitioned_hypergraph,
              std::unique_ptr<IRefiner>& label_propagation,
              std::unique_ptr<IRefiner>& fm,
              kahypar::Metrics& current_metrics,
              const double time_limit) {

    if ( debug && _top_level ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
          << ", imbalance = " << current_metrics.imbalance;
    }

    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;

      if ( label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        label_propagation->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_lp_refiner");

        utils::Timer::instance().start_timer("label_propagation", "Label Propagation");
        improvement_found |= label_propagation->refine(partitioned_hypergraph, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("label_propagation");
      }

      if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        fm->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_fm_refiner");

        utils::Timer::instance().start_timer("fm", "FM");
        improvement_found |= fm->refine(partitioned_hypergraph, current_metrics, time_limit);
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

  double MultilevelCoarsenerBase::refinementTimeLimit(const Level& level) const {
    if ( _context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
      const double time_limit_factor = std::max(1.0,  _context.refinement.fm.time_limit_factor * _context.partition.k);
      return std::max(5.0, time_limit_factor * level.coarseningTime());
    } else {
      return std::numeric_limits<double>::max();
    }
  }


}