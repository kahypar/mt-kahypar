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
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {

  void MultilevelUncoarsener::initializeImpl() {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    _current_metrics = initializeMetrics(partitioned_hg);
    initializeRefinementAlgorithms();

    if (_context.type == ContextType::main) {
      _context.initial_km1 = _current_metrics.km1;
    }

    // Enable progress bar if verbose output is enabled
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug ) {
      _progress.enable();
      _progress.setObjective(_context.partition.objective == Objective::km1
        ? _current_metrics.km1 : _current_metrics.cut);
    }

    _current_level = _uncoarseningData.hierarchy.size();
    _num_levels = _current_level;
  }

  bool MultilevelUncoarsener::isTopLevelImpl() const {
    return _current_level < 0;
  }

  void MultilevelUncoarsener::projectToNextLevelAndRefineImpl() {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    if ( _current_level == _num_levels ) {
      // We always start with a refinement pass on the smallest hypergraph.
      // The next calls to this function will then project the partition to the next level
      // and perform refinement until we reach the input hypergraph.
      const double time_limit = refinementTimeLimit(_context, _uncoarseningData.hierarchy.back().coarseningTime());
      refine(partitioned_hg, time_limit);
      _progress.setObjective(
        _current_metrics.getMetric(Mode::direct, _context.partition.objective));
      _progress += partitioned_hg.initialNumNodes();
    } else {
      ASSERT(_current_level >= 0);
      // Project partition to the hypergraph on the next level of the hierarchy
      _timer.start_timer("projecting_partition", "Projecting Partition");
      const size_t num_nodes_on_previous_level = partitioned_hg.initialNumNodes();
      if (_current_level == 0) {
        partitioned_hg.setHypergraph(_hg);
      } else {
        partitioned_hg.setHypergraph((_uncoarseningData.hierarchy)[_current_level-1].contractedHypergraph());
      }
      // Hypergraph stores partition from previous level.
      // We now extract the block ids and reset the partition to reuse the
      // data structure for the current level.
      partitioned_hg.extractPartIDs(_block_ids);
      partitioned_hg.resetData();

      // Assign nodes of current level to their corresponding representative of the previous level
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = (_uncoarseningData.hierarchy)[_current_level].mapToContractedHypergraph(hn);
        const PartitionID block = _block_ids[coarse_hn];
        ASSERT(block != kInvalidPartition && block < partitioned_hg.k());
        partitioned_hg.setOnlyNodePart(hn, block);
      });
      partitioned_hg.initializePartition();
      _timer.stop_timer("projecting_partition");

      // Improve partition
      const double time_limit = refinementTimeLimit(_context, (_uncoarseningData.hierarchy)[_current_level].coarseningTime());
      refine(partitioned_hg, time_limit);

      // Update Progress Bar
      _progress.setObjective(
        _current_metrics.getMetric(Mode::direct, _context.partition.objective));
      _progress += partitioned_hg.initialNumNodes() - num_nodes_on_previous_level;
    }

    ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
            _current_metrics.getMetric(Mode::direct, _context.partition.objective),
            V(_current_metrics.getMetric(Mode::direct, _context.partition.objective))
            << V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));

    --_current_level;
  }

  void MultilevelUncoarsener::rebalancingImpl() {
    // If we reach the top-level hypergraph and the partition is still imbalanced,
    // we use a rebalancing algorithm to restore balance.
    if (_context.type == ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
      const HyperedgeWeight quality_before = _current_metrics.getMetric(
        Mode::direct, _context.partition.objective);
      if (_context.partition.verbose_output) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
        << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(*_uncoarseningData.partitioned_hg, _context);
      }

      if ( !_context.partition.deterministic ) {
        if (_context.partition.verbose_output) {
          LOG << RED << "Start rebalancing!" << END;
        }

        // Preform rebalancing
        _timer.start_timer("rebalance", "Rebalance");
        if (_context.partition.objective == Objective::km1) {
          Km1Rebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(_current_metrics);
        } else if (_context.partition.objective == Objective::cut) {
          CutRebalancer rebalancer(*_uncoarseningData.partitioned_hg, _context);
          rebalancer.rebalance(_current_metrics);
        }
        _timer.stop_timer("rebalance");

        const HyperedgeWeight quality_after = _current_metrics.getMetric(
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
      } else {
        if (_context.partition.verbose_output) {
          LOG << RED << "Skip rebalancing since deterministic mode is activated" << END;
        }
      }


      ASSERT(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective) ==
              _current_metrics.getMetric(Mode::direct, _context.partition.objective),
              V(_current_metrics.getMetric(Mode::direct, _context.partition.objective))
              << V(metrics::objective(*_uncoarseningData.partitioned_hg, _context.partition.objective)));
    }
  }

  HypernodeID MultilevelUncoarsener::currentNumberOfNodesImpl() const {
    return _uncoarseningData.partitioned_hg->initialNumNodes();
  }

  PartitionedHypergraph&& MultilevelUncoarsener::movePartitionedHypergraphImpl() {
    ASSERT(isTopLevelImpl());
    return std::move(*_uncoarseningData.partitioned_hg);
  }

  void MultilevelUncoarsener::refine(
    PartitionedHypergraph& partitioned_hypergraph,
    const double time_limit) {

    if ( debug && _context.type == ContextType::main ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Refinement Hypergraph", false);
      DBG << "Start Refinement - km1 = " << _current_metrics.km1
      << ", imbalance = " << _current_metrics.imbalance;
    }

    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    while( improvement_found ) {
      improvement_found = false;
      const HyperedgeWeight metric_before = _current_metrics.getMetric(
        Mode::direct, _context.partition.objective);

      if ( _label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        _label_propagation->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_lp_refiner");

        _timer.start_timer("label_propagation", "Label Propagation");
        improvement_found |= _label_propagation->refine(partitioned_hypergraph, dummy, _current_metrics, time_limit);
        _timer.stop_timer("label_propagation");
      }

      if ( _fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        _fm->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_fm_refiner");

        _timer.start_timer("fm", "FM");
        improvement_found |= _fm->refine(partitioned_hypergraph, dummy, _current_metrics, time_limit);
        _timer.stop_timer("fm");
      }

      if ( _flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
        _flows->initialize(partitioned_hypergraph);
        _timer.stop_timer("initialize_flow_scheduler");

        _timer.start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
        improvement_found |= _flows->refine(partitioned_hypergraph, dummy, _current_metrics, time_limit);
        _timer.stop_timer("flow_refinement_scheduler");
      }

      if ( _context.type == ContextType::main ) {
        ASSERT(_current_metrics.getMetric(Mode::direct, _context.partition.objective)
               == metrics::objective(partitioned_hypergraph, _context.partition.objective),
               "Actual metric" << V(metrics::km1(partitioned_hypergraph))
               << "does not match the metric updated by the refiners" << V(_current_metrics.km1));
      }

      const HyperedgeWeight metric_after = _current_metrics.getMetric(
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
