/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <mt-kahypar/partition/coarsening/multilevel_uncoarsener.h>
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {

  PartitionedHypergraph&& MultilevelUncoarsener::doUncoarsen(
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm) {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    kahypar::Metrics current_metrics = initialize(partitioned_hg);

    if (_context.type == kahypar::ContextType::main) {
      _context.initial_km1 = current_metrics.km1;
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == kahypar::Objective::km1
                                              ? current_metrics.km1 : current_metrics.cut,
                                              _context.partition.verbose_output &&
                                              _context.partition.enable_progress_bar && !debug);
    uncontraction_progress += partitioned_hg.initialNumNodes();

    // Refine Coarsest Partitioned Hypergraph
    double time_limit = refinementTimeLimit(_context, _uncoarseningData.hierarchy.back().coarseningTime());
    refine(partitioned_hg, label_propagation, fm, current_metrics, time_limit);

    ds::Array<PartitionID> part_ids(_hg.initialNumNodes(), kInvalidPartition);

    if (!_uncoarseningData.hierarchy.empty()) {
      // restore partition of initial partitioning on larger phg
      ds::Array<PartitionID> ip_part_ids;
      partitioned_hg.extractPartIDs(ip_part_ids);

      partitioned_hg = PartitionedHypergraph(
        _context.partition.k, _hg, parallel_tag_t());
      partitioned_hg.setHypergraph(
        _uncoarseningData.hierarchy.back().contractedHypergraph());
      tbb::parallel_for(0UL, ip_part_ids.size(), [&](const HypernodeID hn) {
        part_ids[hn] = ip_part_ids[hn];
      });
    }

    for (int i = _uncoarseningData.hierarchy.size() - 1; i >= 0; --i) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      const size_t num_nodes = partitioned_hg.initialNumNodes();
      if (i == 0) {
        partitioned_hg.setHypergraph(_hg);
      } else {
        partitioned_hg.setHypergraph((_uncoarseningData.hierarchy)[i-1].contractedHypergraph());
      }
      // extract part_ids to reset partition
      if (static_cast<size_t>(i) < _uncoarseningData.hierarchy.size() - 1) {
        partitioned_hg.extractPartIDs(part_ids);
        partitioned_hg.resetData();
      }

      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = (_uncoarseningData.hierarchy)[i].mapToContractedHypergraph(hn);
        const PartitionID block = part_ids[coarse_hn];
        ASSERT(block != kInvalidPartition && block < partitioned_hg.k());
        partitioned_hg.setOnlyNodePart(hn, block);
      });
      partitioned_hg.initializePartition();
      utils::Timer::instance().stop_timer("projecting_partition");

      // Refinement
      time_limit = refinementTimeLimit(_context, (_uncoarseningData.hierarchy)[i].coarseningTime());
      refine(partitioned_hg, label_propagation, fm, current_metrics, time_limit);

      // Update Progress Bar
      uncontraction_progress.setObjective(
        current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective));
      uncontraction_progress += partitioned_hg.initialNumNodes() - num_nodes;
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if (_context.type == kahypar::ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if (_context.partition.verbose_output) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
        << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(*_uncoarseningData.partitioned_hg, _context);
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
          Km1Rebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        } else if (_context.partition.objective == kahypar::Objective::cut) {
          CutRebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        }
        utils::Timer::instance().stop_timer("rebalance");

        const HyperedgeWeight quality_after = current_metrics.getMetric(
          kahypar::Mode::direct_kway, _context.partition.objective);
        if (_context.partition.verbose_output) {
          const HyperedgeWeight quality_delta = quality_after - quality_before;
          if (quality_delta > 0) {
            LOG << RED << "Rebalancer decreased solution quality by" << quality_delta
            << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
          } else {
            LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
            << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
          }
        }
      }

      ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
             current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
             V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective))
             << V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));
    }
    return std::move(*_uncoarseningData.partitioned_hg);
  }

  void MultilevelUncoarsener::refine(
    PartitionedHypergraph& partitioned_hypergraph,
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm,
    kahypar::Metrics& current_metrics,
    const double time_limit) {

    if ( debug && _context.type == kahypar::ContextType::main ) {
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

      if ( _context.type == kahypar::ContextType::main ) {
        ASSERT(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)
               == metrics::objective(partitioned_hypergraph, _context.partition.objective),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      if ( !_context.refinement.refine_until_no_improvement ) {
        break;
      }
    }

    if ( _context.type == kahypar::ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
  }

}
