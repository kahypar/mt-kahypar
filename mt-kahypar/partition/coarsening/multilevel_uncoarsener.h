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
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/partition/coarsening/level.h"
#include "mt-kahypar/partition/coarsening/i_uncoarsener.h"
namespace mt_kahypar {

  class MultilevelUncoarsener : public IUncoarsener {

  private:
  static constexpr bool debug = false;

  public:
    MultilevelUncoarsener(Hypergraph& hypergraph,
                        const Context& context,
                        const bool top_level,
                        UncoarseningData& uncoarseningData) :
      _hg(hypergraph),
      _partitioned_hg(std::move(uncoarseningData.partitioned_hypergraph)),
      _context(context),
      _top_level(top_level),
      _hierarchy(std::move(uncoarseningData.hierarchy)) { }

  MultilevelUncoarsener(const MultilevelUncoarsener&) = delete;
  MultilevelUncoarsener(MultilevelUncoarsener&&) = delete;
  MultilevelUncoarsener & operator= (const MultilevelUncoarsener &) = delete;
  MultilevelUncoarsener & operator= (MultilevelUncoarsener &&) = delete;

  private:
    Hypergraph& _hg;
    std::shared_ptr<PartitionedHypergraph> _partitioned_hg;
    const Context& _context;
    const bool _top_level;
    std::shared_ptr<vec<Level>> _hierarchy;

protected:
    PartitionedHypergraph&& doUncoarsen(
      std::unique_ptr<IRefiner>& label_propagation,
      std::unique_ptr<IRefiner>& fm) {
      PartitionedHypergraph& coarsest_hg = currentPartitionedHypergraph();
      kahypar::Metrics current_metrics = initialize(coarsest_hg);

      if (_top_level) {
        _context.initial_km1 = current_metrics.km1;
      }

      utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                                _context.partition.objective == kahypar::Objective::km1
                                                ? current_metrics.km1 : current_metrics.cut,
                                                _context.partition.verbose_output &&
                                                _context.partition.enable_progress_bar && !debug);
      uncontraction_progress += coarsest_hg.initialNumNodes();

      // Refine Coarsest Partitioned Hypergraph
      double time_limit = refinementTimeLimit(_context, _hierarchy->back().coarseningTime());
      refine(coarsest_hg, label_propagation, fm, current_metrics, time_limit);

      for (int i = _hierarchy->size() - 1; i >= 0; --i) {
        // Project partition to next level finer hypergraph
        utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
        PartitionedHypergraph& representative_hg = (*_hierarchy)[i].representativeHypergraph();
        PartitionedHypergraph& contracted_hg = (*_hierarchy)[i].contractedPartitionedHypergraph();
        representative_hg.doParallelForAllNodes([&](const HypernodeID hn) {
          const HypernodeID coarse_hn = (*_hierarchy)[i].mapToContractedHypergraph(hn);
          const PartitionID block = contracted_hg.partID(coarse_hn);
          ASSERT(block != kInvalidPartition && block < representative_hg.k());
          representative_hg.setOnlyNodePart(hn, block);
        });
        representative_hg.initializePartition();

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
        time_limit = refinementTimeLimit(_context, (*_hierarchy)[i].coarseningTime());
        refine(representative_hg, label_propagation, fm, current_metrics, time_limit);

        // Update Progress Bar
        uncontraction_progress.setObjective(
          current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective));
        uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
      }

      // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
      if (_top_level && !metrics::isBalanced(*_partitioned_hg, _context)) {
        const HyperedgeWeight quality_before = current_metrics.getMetric(
          kahypar::Mode::direct_kway, _context.partition.objective);
        if (_context.partition.verbose_output) {
          LOG << RED << "Partition is imbalanced (Current Imbalance:"
          << metrics::imbalance(*_partitioned_hg, _context) << ")" << END;

          LOG << "Part weights: (violations in red)";
          io::printPartWeightsAndSizes(*_partitioned_hg, _context);
        }

        if (_context.partition.deterministic) {
          if (_context.partition.verbose_output) {
            LOG << RED << "Skip rebalancing since deterministic mode is activated" << END;
          }
        } else {
          if (_context.partition.verbose_output) {
            LOG << RED << "Start rebalancing!" << END;
          }
          utils::Timer::instance().start_timer("rebalance", "Rebalance");
          if (_context.partition.objective == kahypar::Objective::km1) {
            Km1Rebalancer rebalancer(*_partitioned_hg, _context);
            rebalancer.rebalance(current_metrics);
          } else if (_context.partition.objective == kahypar::Objective::cut) {
            CutRebalancer rebalancer(*_partitioned_hg, _context);
            rebalancer.rebalance(current_metrics);
          }
          utils::Timer::instance().stop_timer("rebalance");

          const HyperedgeWeight quality_after = current_metrics.getMetric(
            kahypar::Mode::direct_kway, _context.partition.objective);
          if (_context.partition.verbose_output) {
            const HyperedgeWeight quality_delta = quality_after - quality_before;
            if (quality_delta > 0) {
              LOG << RED << "Rebalancer decreased solution quality by" << quality_delta
              << "(Current Imbalance:" << metrics::imbalance(*_partitioned_hg, _context) << ")" << END;
            } else {
              LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
              << "(Current Imbalance:" << metrics::imbalance(*_partitioned_hg, _context) << ")" << END;
            }
          }
        }

        ASSERT(metrics::objective(*_partitioned_hg, _context.partition.objective) ==
               current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
               V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective))
               << V(metrics::objective(*_partitioned_hg, _context.partition.objective)));
      }
      return std::move(*_partitioned_hg);
    }

  PartitionedHypergraph& currentPartitionedHypergraph() {
    /*ASSERT(_is_finalized);*/
    if ( _hierarchy->empty() ) {
      return *_partitioned_hg;
    } else {
      return _hierarchy->back().contractedPartitionedHypergraph();
    }
  }
  kahypar::Metrics initialize(PartitionedHypergraph& phg) {
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
  double refinementTimeLimit(const Context& context, const double time) {
    if ( context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
      const double time_limit_factor = std::max(1.0,  context.refinement.fm.time_limit_factor * context.partition.k);
      return std::max(5.0, time_limit_factor * time);
    } else {
      return std::numeric_limits<double>::max();
    }
  }
  void refine(
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

private:
  PartitionedHypergraph&& uncoarsenImpl(
      std::unique_ptr<IRefiner>& label_propagation,
      std::unique_ptr<IRefiner>& fm) override {
    return doUncoarsen(label_propagation, fm);
  }
  };

}
