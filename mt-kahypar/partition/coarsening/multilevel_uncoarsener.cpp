/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#include <mt-kahypar/partition/coarsening/multilevel_uncoarsener.h>
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flows/scheduler.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {

  PartitionedHypergraph&& MultilevelUncoarsener::doUncoarsen(
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm) {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    Metrics current_metrics = initialize(partitioned_hg);

    if (_context.type == ContextType::main) {
      _context.initial_km1 = current_metrics.km1;
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == Objective::km1
                                              ? current_metrics.km1 : current_metrics.cut,
                                              _context.partition.verbose_output &&
                                              _context.partition.enable_progress_bar && !debug);

    // Initialize Flow Refinement Scheduler
    std::unique_ptr<IRefiner> flows(nullptr);
    if ( _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
      flows = std::make_unique<FlowRefinementScheduler>(_hg, _context);
    }

    // Refine Coarsest Partitioned Hypergraph
    double time_limit = refinementTimeLimit(_context, _uncoarseningData.hierarchy.back().coarseningTime());
    refine(partitioned_hg, label_propagation, fm, flows, current_metrics, time_limit);
    uncontraction_progress.setObjective(
      current_metrics.getMetric(Mode::direct, _context.partition.objective));
    uncontraction_progress += partitioned_hg.initialNumNodes();

    ds::Array<PartIdType> part_ids(_hg.initialNumNodes(), PartIdType(kInvalidPartition));

    for (int i = _uncoarseningData.hierarchy.size() - 1; i >= 0; --i) {
      // Project partition to next level finer hypergraph
      _timer.start_timer("projecting_partition", "Projecting Partition");
      const size_t num_nodes = partitioned_hg.initialNumNodes();
      if (i == 0) {
        partitioned_hg.setHypergraph(_hg);
      } else {
        partitioned_hg.setHypergraph((_uncoarseningData.hierarchy)[i-1].contractedHypergraph());
      }
      // extract part_ids to reset partition
      partitioned_hg.extractPartIDs(part_ids);
      partitioned_hg.resetData();

      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = (_uncoarseningData.hierarchy)[i].mapToContractedHypergraph(hn);
        const PartitionID block = part_ids[coarse_hn];
        ASSERT(block != kInvalidPartition && block < partitioned_hg.k());
        partitioned_hg.setOnlyNodePart(hn, block);
      });
      partitioned_hg.initializePartition();
      _timer.stop_timer("projecting_partition");

      // Refinement
      time_limit = refinementTimeLimit(_context, (_uncoarseningData.hierarchy)[i].coarseningTime());
      refine(partitioned_hg, label_propagation, fm, flows, current_metrics, time_limit);

      // Update Progress Bar
      uncontraction_progress.setObjective(
        current_metrics.getMetric(Mode::direct, _context.partition.objective));
      uncontraction_progress += partitioned_hg.initialNumNodes() - num_nodes;
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if (_context.type == ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        Mode::direct, _context.partition.objective);
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
        _timer.start_timer("rebalance", "Rebalance");
        if (_context.partition.objective == Objective::km1) {
          Km1Rebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        } else if (_context.partition.objective == Objective::cut) {
          CutRebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(current_metrics);
        }
        _timer.stop_timer("rebalance");

        const HyperedgeWeight quality_after = current_metrics.getMetric(
          Mode::direct, _context.partition.objective);
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
    }

    ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
            current_metrics.getMetric(Mode::direct, _context.partition.objective),
            V(current_metrics.getMetric(Mode::direct, _context.partition.objective))
            << V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));
    return std::move(*_uncoarseningData.partitioned_hg);
  }

  void MultilevelUncoarsener::refine(
    PartitionedHypergraph& partitioned_hypergraph,
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm,
    std::unique_ptr<IRefiner>& flows,
    Metrics& current_metrics,
    const double time_limit) {

    if ( debug && _context.type == ContextType::main ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << current_metrics.km1
      << ", imbalance = " << current_metrics.imbalance;
    }

    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;
      const HyperedgeWeight metric_before = current_metrics.getMetric(
        Mode::direct, _context.partition.objective);

      if ( label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        label_propagation->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_lp_refiner");

        _timer.start_timer("label_propagation", "Label Propagation");
        improvement_found |= label_propagation->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        _timer.stop_timer("label_propagation");
      }

      if ( fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        fm->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_fm_refiner");

        _timer.start_timer("fm", "FM");
        improvement_found |= fm->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        _timer.stop_timer("fm");
      }

      if ( flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
        flows->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_flow_scheduler");

        _timer.start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
        improvement_found |= flows->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        _timer.stop_timer("flow_refinement_scheduler");
      }

      if ( _context.type == ContextType::main ) {
        ASSERT(current_metrics.getMetric(Mode::direct, _context.partition.objective)
               == metrics::objective(partitioned_hypergraph, _context.partition.objective),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
               << "does not match the metric updated by the refiners" << V(current_metrics.km1));
      }

      const HyperedgeWeight metric_after = current_metrics.getMetric(
        Mode::direct, _context.partition.objective);
      const double relative_improvement = 1.0 -
        static_cast<double>(metric_after) / metric_before;
      if ( !_context.refinement.refine_until_no_improvement ||
           relative_improvement <= _context.refinement.relative_improvement_threshold ) {
        break;
      }
    }

    if ( _context.type == ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
  }

}
