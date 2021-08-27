/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/partition/coarsening/i_uncoarsener.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/partition/coarsening/level.h"
namespace mt_kahypar {

  class NLevelUncoarsener : public IUncoarsener {

  private:
    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    using ParallelHyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>>;

  public:
    NLevelUncoarsener(Hypergraph& hypergraph,
                      const Context& context,
                      const bool top_level,
                      UncoarseningData& uncoarseningData) :
      _hg(hypergraph),
      _context(context),
      _top_level(top_level),
      _phg(std::move(uncoarseningData.partitioned_hypergraph)),
      _compactified_hg(std::move(uncoarseningData.compactified_hg)),
      _compactified_phg(std::move(uncoarseningData.compactified_phg)),
      _compactified_hn_mapping(std::move(uncoarseningData.compactified_hn_mapping)),
      _hierarchy(std::move(uncoarseningData.n_level_hierarchy)),
      _removed_hyperedges_batches(std::move(uncoarseningData.removed_hyperedges_batches)),
      _round_coarsening_times(std::move(uncoarseningData.round_coarsening_times)) { }

    NLevelUncoarsener(const NLevelUncoarsener&) = delete;
    NLevelUncoarsener(NLevelUncoarsener&&) = delete;
    NLevelUncoarsener & operator= (const NLevelUncoarsener &) = delete;
    NLevelUncoarsener & operator= (NLevelUncoarsener &&) = delete;

  private:
    Hypergraph& _hg;
    const Context& _context;
    const bool _top_level;
    std::shared_ptr<PartitionedHypergraph> _phg;
    std::shared_ptr<Hypergraph> _compactified_hg;
    std::shared_ptr<PartitionedHypergraph> _compactified_phg;
    std::shared_ptr<vec<HypernodeID>> _compactified_hn_mapping;
    std::shared_ptr<VersionedBatchVector> _hierarchy;
    std::shared_ptr<ParallelHyperedgeVector> _removed_hyperedges_batches;
    std::shared_ptr<vec<double>> _round_coarsening_times;

  protected:

    PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) {
      /*ASSERT(_is_finalized);*/
      kahypar::Metrics current_metrics = initialize(*_compactified_phg);
      if (_top_level) {
        _context.initial_km1 = current_metrics.km1;
      }

      // Project partition from compactified hypergraph to original hypergraph
      utils::Timer::instance().start_timer("initialize_partition", "Initialize Partition");
      *_phg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
      _phg->doParallelForAllNodes([&](const HypernodeID hn) {
        ASSERT(static_cast<size_t>(hn) < _compactified_hn_mapping->size());
        const HypernodeID compactified_hn = (*_compactified_hn_mapping)[hn];
        const PartitionID block_id = _compactified_phg->partID(compactified_hn);
        ASSERT(block_id != kInvalidPartition && block_id < _context.partition.k);
        _phg->setOnlyNodePart(hn, block_id);
      });
      _phg->initializePartition();

      if ( _context.refinement.fm.algorithm == FMAlgorithm::fm_gain_cache ) {
        _phg->initializeGainCache();
      }

      ASSERT(metrics::objective(*_compactified_phg, _context.partition.objective) ==
             metrics::objective(*_phg, _context.partition.objective),
             V(metrics::objective(*_compactified_phg, _context.partition.objective)) <<
             V(metrics::objective(*_phg, _context.partition.objective)));
      ASSERT(metrics::imbalance(*_compactified_phg, _context) ==
             metrics::imbalance(*_phg, _context),
             V(metrics::imbalance(*_compactified_phg, _context)) <<
             V(metrics::imbalance(*_phg, _context)));
      utils::Timer::instance().stop_timer("initialize_partition");

      utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                                _context.partition.objective == kahypar::Objective::km1 ? current_metrics.km1 : current_metrics.cut,
                                                _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug);
      uncontraction_progress += _compactified_hg->initialNumNodes();

      // Initialize Refiner
      if ( label_propagation ) {
        label_propagation->initialize(*_phg);
      }
      if ( fm ) {
        fm->initialize(*_phg);
      }

      // Perform batch uncontractions
      bool is_timer_disabled = false;
      bool force_measure_timings = _context.partition.measure_detailed_uncontraction_timings && _top_level;
      if ( utils::Timer::instance().isEnabled() ) {
        utils::Timer::instance().disable();
        is_timer_disabled = true;
      }

      ASSERT(_round_coarsening_times->size() == _removed_hyperedges_batches->size());
      _round_coarsening_times->push_back(_round_coarsening_times->size() > 0 ?
                                        _round_coarsening_times->back() : std::numeric_limits<double>::max()); // Sentinel

      size_t num_batches = 0;
      size_t total_batches_size = 0;
      const size_t minimum_required_number_of_border_vertices = std::max(_context.refinement.max_batch_size,
                                                                         _context.shared_memory.num_threads * _context.refinement.min_border_vertices_per_thread);
      ds::StreamingVector<HypernodeID> tmp_refinement_nodes;
      kahypar::ds::FastResetFlagArray<> border_vertices_of_batch(_phg->initialNumNodes());
      auto do_localized_refinement = [&]() {
        parallel::scalable_vector<HypernodeID> refinement_nodes = tmp_refinement_nodes.copy_parallel();
        tmp_refinement_nodes.clear_parallel();
        border_vertices_of_batch.reset();
        localizedRefine(*_phg, refinement_nodes, label_propagation,
                        fm, current_metrics, force_measure_timings);
      };

      while ( !_hierarchy->empty() ) {
        BatchVector& batches = _hierarchy->back();

        // Uncontract all batches of a specific version of the hypergraph
        while ( !batches.empty() ) {
          const Batch& batch = batches.back();
          if ( batch.size() > 0 ) {
            HEAVY_REFINEMENT_ASSERT(metrics::objective(*_phg, _context.partition.objective) ==
                                    current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                                    V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                                    V(metrics::objective(*_phg, _context.partition.objective)));
            utils::Timer::instance().start_timer("batch_uncontractions", "Batch Uncontractions", false, force_measure_timings);
            _phg->uncontract(batch);
            utils::Timer::instance().stop_timer("batch_uncontractions", force_measure_timings);
            HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
            HEAVY_REFINEMENT_ASSERT(_phg->checkTrackedPartitionInformation());
            HEAVY_REFINEMENT_ASSERT(metrics::objective(*_phg, _context.partition.objective) ==
                                    current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
                                    V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
                                    V(metrics::objective(*_phg, _context.partition.objective)));

            utils::Timer::instance().start_timer("collect_border_vertices", "Collect Border Vertices", false, force_measure_timings);
            tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
              const Memento& memento = batch[i];
              if ( !border_vertices_of_batch[memento.u] && _phg->isBorderNode(memento.u) ) {
                border_vertices_of_batch.set(memento.u, true);
                tmp_refinement_nodes.stream(memento.u);
              }
              if ( !border_vertices_of_batch[memento.v] && _phg->isBorderNode(memento.v) ) {
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
        if ( !_removed_hyperedges_batches->empty() ) {
          utils::Timer::instance().start_timer("restore_single_pin_and_parallel_nets", "Restore Single Pin and Parallel Nets", false, force_measure_timings);
          _phg->restoreSinglePinAndParallelNets(_removed_hyperedges_batches->back());
          _removed_hyperedges_batches->pop_back();
          utils::Timer::instance().stop_timer("restore_single_pin_and_parallel_nets", force_measure_timings);
          HEAVY_REFINEMENT_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
          HEAVY_REFINEMENT_ASSERT(_phg->checkTrackedPartitionInformation());

          // Perform refinement on all vertices
          const double time_limit = refinementTimeLimit(_context, _round_coarsening_times->back());
          globalRefine(*_phg, fm, current_metrics, time_limit);
          uncontraction_progress.setObjective(current_metrics.getMetric(
              _context.partition.mode, _context.partition.objective));
          _round_coarsening_times->pop_back();
        }
        _hierarchy->pop_back();
      }

      // Top-Level Refinement on all vertices
      const HyperedgeWeight objective_before = current_metrics.getMetric(
        _context.partition.mode, _context.partition.objective);
      const double time_limit = refinementTimeLimit(_context, _round_coarsening_times->back());
      globalRefine(*_phg, fm, current_metrics, time_limit);
      _round_coarsening_times->pop_back();
      ASSERT(_round_coarsening_times->size() == 0);
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
      if ( _top_level && !metrics::isBalanced(*_phg, _context)) {
        const HyperedgeWeight quality_before = current_metrics.getMetric(
          kahypar::Mode::direct_kway, _context.partition.objective);
        if ( _context.partition.verbose_output ) {
          LOG << RED << "Partition is imbalanced (Current Imbalance:"
          << metrics::imbalance(*_phg, _context) << ") ->"
          << "Rebalancer is activated" << END;

          LOG << "Part weights: (violations in red)";
          io::printPartWeightsAndSizes(*_phg, _context);
        }

        utils::Timer::instance().start_timer("rebalance", "Rebalance");
        if ( _context.partition.objective == kahypar::Objective::km1 ) {
          Km1Rebalancer rebalancer(*_phg, _context);
          rebalancer.rebalance(current_metrics);
        } else if ( _context.partition.objective == kahypar::Objective::cut ) {
          CutRebalancer rebalancer(*_phg, _context);
          rebalancer.rebalance(current_metrics);
        }
        utils::Timer::instance().stop_timer("rebalance");

        const HyperedgeWeight quality_after = current_metrics.getMetric(
          kahypar::Mode::direct_kway, _context.partition.objective);
        if ( _context.partition.verbose_output ) {
          const HyperedgeWeight quality_delta = quality_after - quality_before;
          if ( quality_delta > 0 ) {
            LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
            << "(Current Imbalance:" << metrics::imbalance(*_phg, _context) << ")" << END;
          } else {
            LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
            << "(Current Imbalance:" << metrics::imbalance(*_phg, _context) << ")" << END;
          }
        }
      }

      double avg_batch_size = static_cast<double>(total_batches_size) / num_batches;
      utils::Stats::instance().add_stat("num_batches", static_cast<int64_t>(num_batches));
      utils::Stats::instance().add_stat("avg_batch_size", avg_batch_size);
      DBG << V(num_batches) << V(avg_batch_size);

      ASSERT(metrics::objective(*_phg, _context.partition.objective) ==
             current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
             V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
             V(metrics::objective(*_phg, _context.partition.objective)));

      return std::move(*_phg);
    }

    void localizedRefine(PartitionedHypergraph& partitioned_hypergraph,
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

    NLevelGlobalFMParameters applyGlobalFMParameters(const FMParameters& fm,
                                                     const NLevelGlobalFMParameters global_fm) {
      NLevelGlobalFMParameters tmp_global_fm;
      tmp_global_fm.num_seed_nodes = fm.num_seed_nodes;
      tmp_global_fm.obey_minimal_parallelism = fm.obey_minimal_parallelism;
      fm.num_seed_nodes = global_fm.num_seed_nodes;
      fm.obey_minimal_parallelism = global_fm.obey_minimal_parallelism;
      return tmp_global_fm;
    }

    void globalRefine(PartitionedHypergraph& partitioned_hypergraph,
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

    double refinementTimeLimit(const Context& context, const double time) {
      if ( context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
        const double time_limit_factor = std::max(1.0,  context.refinement.fm.time_limit_factor * context.partition.k);
        return std::max(5.0, time_limit_factor * time);
      } else {
        return std::numeric_limits<double>::max();
      }
    }

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

  private:
  PartitionedHypergraph&& uncoarsenImpl(
      std::unique_ptr<IRefiner>& label_propagation,
      std::unique_ptr<IRefiner>& fm) override {
    return doUncoarsen(label_propagation, fm);
  }
  };
}
