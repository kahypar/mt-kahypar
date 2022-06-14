/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <mt-kahypar/partition/coarsening/multilevel_uncoarsener.h>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flows/scheduler.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"
#include "mt-kahypar/utils/progress_bar.h"
#include <tbb/parallel_sort.h>
#include <tbb/enumerable_thread_specific.h>

namespace mt_kahypar {
  using ds::Array;
  using ds::SeparatedNodes;

  PartitionedHypergraph&& MultilevelUncoarsener::doUncoarsen(
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm) {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    Metrics current_metrics = initialize(partitioned_hg);

    if (_context.type == kahypar::ContextType::main) {
      _context.initial_km1 = current_metrics.km1;
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
                                              _context.partition.objective == kahypar::Objective::km1
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
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
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
        if (coarse_hn != kInvalidHypernode) {
          const PartitionID block = part_ids[coarse_hn];
          ASSERT(block != kInvalidPartition && block < partitioned_hg.k());
          partitioned_hg.setOnlyNodePart(hn, block);
        }
      });
      partitioned_hg.initializePartition();
      utils::Timer::instance().stop_timer("projecting_partition");

      if (partitioned_hg.hasSeparatedNodes()) {
        utils::Timer::instance().start_timer("assign_separated_nodes", "Assign Separated Nodes");
        SeparatedNodes& separated_nodes = partitioned_hg.separatedNodes();
        tbb::enumerable_thread_specific<Array<HyperedgeWeight>> edge_weights(partitioned_hg.k());
        tbb::enumerable_thread_specific<Array<HypernodeWeight>> block_weights(partitioned_hg.k());

        const HypernodeID first_separated = separated_nodes.currentBatchIndex();
        const HypernodeID last_separated = separated_nodes.numNodes();
        Array<HyperedgeWeight> max_gains(last_separated - first_separated);
        tbb::parallel_for(first_separated, last_separated, [&](const HypernodeID s_node) {
          const HypernodeID node = separated_nodes.originalHypernodeID(s_node);
          Array<HyperedgeWeight>& local_edge_weights = edge_weights.local();
          local_edge_weights.assign(partitioned_hg.k(), 0);

          for (HyperedgeID e: partitioned_hg.incidentEdges(node)) {
            const PartitionID target_part = partitioned_hg.partID(partitioned_hg.edgeTarget(e));
            if (target_part != kInvalidPartition) {
              local_edge_weights[target_part] += partitioned_hg.edgeWeight(e);
            }
          }

          HyperedgeWeight max_gain = 0;
          for (PartitionID part = 0; part < partitioned_hg.k(); ++part) {
            if (local_edge_weights[part] >= max_gain) {
              max_gain = local_edge_weights[part];
            }
          }
          max_gains[s_node - first_separated] = max_gain;
        });

        Array<HypernodeID> id_order(last_separated - first_separated);
        tbb::parallel_for(first_separated, last_separated, [&](const HypernodeID s_node) {
          id_order[s_node - first_separated] = s_node;
        });

        tbb::parallel_sort(id_order.begin(), id_order.end(),
          [&](const HypernodeID& left, const HypernodeID& right) {
            return max_gains[left - first_separated] > max_gains[right - first_separated];
          }
        );

        for (size_t i = 0; i < id_order.size(); ++i) {
          const HypernodeID s_node = id_order[i];
          const HypernodeID node = separated_nodes.originalHypernodeID(s_node);

          if (separated_nodes.partID(s_node) != kInvalidPartition) {
            // Note: Because the block weight already respects separated nodes, we use
            // setOnlyNodePart instead of setNodePart.
            // Later, the nodes are removed via popBatch().
            partitioned_hg.setOnlyNodePart(node, separated_nodes.partID(s_node));
            continue;
          }

          Array<HyperedgeWeight>& local_edge_weights = edge_weights.local();
          Array<HypernodeWeight>& local_block_weights = block_weights.local();
          local_edge_weights.assign(partitioned_hg.k(), 0);
          local_block_weights.assign(partitioned_hg.k(), 0);
          for (PartitionID part = 0; part < partitioned_hg.k(); ++part) {
            if (local_block_weights[part] == 0) {
              local_block_weights[part] = partitioned_hg.partWeight(part);
            }
          }

          for (HyperedgeID e: partitioned_hg.incidentEdges(node)) {
            const PartitionID target_part = partitioned_hg.partID(partitioned_hg.edgeTarget(e));
            if (target_part != kInvalidPartition) {
              local_edge_weights[target_part] += partitioned_hg.edgeWeight(e);
            }
          }

          // greedily assign separated nodes
          while (true) {
            PartitionID max_part = kInvalidPartition;
            HyperedgeWeight max_gain = 0;
            for (PartitionID part = 0; part < partitioned_hg.k(); ++part) {
              const HypernodeWeight max_part_weight = _context.partition.max_part_weights[part];
              if (local_edge_weights[part] >= max_gain &&
                  local_block_weights[part] + separated_nodes.nodeWeight(s_node) <= max_part_weight) {
                max_part = part;
                max_gain = local_edge_weights[part];
              }
            }

            if (max_part == kInvalidPartition) {
              partitioned_hg.setNodePart(node, 0);
              break;
            }
            const HypernodeWeight max_allowed_weight = _context.partition.max_part_weights[max_part];
            auto [success, new_part_weight] = partitioned_hg.trySetNodePart(node, max_part, max_allowed_weight);
            if (success) {
              break;
            } else {
              local_block_weights[max_part] = new_part_weight;
            }
          }
        }
        separated_nodes.popBatch();

        tbb::parallel_invoke([&] {
          current_metrics.cut = metrics::hyperedgeCut(partitioned_hg);
        }, [&] {
          current_metrics.km1 = metrics::km1(partitioned_hg);
        }, [&] {
          current_metrics.imbalance = metrics::imbalance(partitioned_hg, _context);
        });

        ASSERT([&] {
          for (HypernodeID node: partitioned_hg.nodes()) {
            if (partitioned_hg.partID(node) == kInvalidPartition) {
              return false;
            }
          }
          return true;
        }() );

        utils::Timer::instance().stop_timer("assign_separated_nodes");
      }

      // Refinement
      time_limit = refinementTimeLimit(_context, (_uncoarseningData.hierarchy)[i].coarseningTime());
      refine(partitioned_hg, label_propagation, fm, flows, current_metrics, time_limit);

      // Update Progress Bar
      uncontraction_progress.setObjective(
        current_metrics.getMetric(Mode::direct, _context.partition.objective));
      uncontraction_progress += partitioned_hg.initialNumNodes() - num_nodes;
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if (_context.type == kahypar::ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
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

    if ( debug && _context.type == kahypar::ContextType::main ) {
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

      if ( flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
        utils::Timer::instance().start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
        flows->initialize(partitioned_hypergraph);
        utils::Timer::instance().stop_timer("initialize_flow_scheduler");

        utils::Timer::instance().start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
        improvement_found |= flows->refine(partitioned_hypergraph, dummy, current_metrics, time_limit);
        utils::Timer::instance().stop_timer("flow_refinement_scheduler");
      }

      if ( _context.type == kahypar::ContextType::main ) {
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

    if ( _context.type == kahypar::ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
  }

}
