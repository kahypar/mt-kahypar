/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <tbb/parallel_for_each.h>

#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/utils/timer.h>

#include <atomic>
#include <external_tools/kahypar/kahypar/partition/metrics.h>

#include "localized_kway_fm_core.h"
#include "global_rollback.h"


namespace mt_kahypar {
namespace refinement {

// TODO try variant in which, a bunch of searches are stored in a PQ, findMoves(..) yields frequently, and then the most promising search is scheduled next

class MultiTryKWayFM {
public:
  MultiTryKWayFM(const Context& context, TaskGroupID taskGroupID, size_t numNodes, size_t numHyperedges) :
          context(context),
          taskGroupID(taskGroupID),
          sharedData(numNodes, context.partition.k, context.shared_memory.num_threads),
          globalRollback(numNodes, numHyperedges, context.partition.k),
          ets_fm(context, numNodes, sharedData.vertexPQHandles.data())
  { }

  Gain refine(PartitionedHypergraph& phg, kahypar::Metrics& metrics) {
    Gain improvement = refine(phg);
    metrics.km1 -= improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    assert(metrics.km1 == metrics::km1(phg));
    return improvement;
  }

  Gain refine(PartitionedHypergraph& phg) {
    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("fm", "FM");
    timer.start_timer("fm_unnecessary_init", "FM Init");

    phg.initializeGainInformation();                // initialization only as long as LP refiner does not use these datastructures
    globalRollback.setRemainingOriginalPins(phg);   // initialization only as long as LP refiner does not use these datastructures

    timer.stop_timer("fm_unnecessary_init");
    //sharedData.partition_weight_budgets.initialize(phg, context.partition.max_part_weights);          // only for version with budgets

    Gain overall_improvement = 0;
    for (size_t round = 0; round < context.refinement.fm.rounds; ++round) {                    // global multi try rounds
      timer.start_timer("collect_border_nodes", "Collect Border Nodes");

      initialize(phg);

      timer.stop_timer("collect_border_nodes");
      timer.start_timer("find_moves", "Find Moves");

      vec<HypernodeWeight> initialPartWeights(size_t(sharedData.numParts));
      for (PartitionID i = 0; i < sharedData.numParts; ++i) initialPartWeights[i] = phg.partWeight(i);

      if (context.refinement.fm.multitry) {
        auto task = [&](const int, const int socket_local_task_id, const int task_id) {
          unused(socket_local_task_id); unused(task_id);
          LocalizedKWayFM& fm = ets_fm.local();
          while(fm.findMoves(phg, sharedData)) { /* keep running */ }
        };
        TBBNumaArena::instance().run_max_concurrency_tasks_on_all_sockets(taskGroupID, task);
        //task(0,0,0);
      } else {
        // Try boundary FM
        LOG << "start FM" << V(sharedData.moveTracker.numPerformedMoves());
        vec<HypernodeID> test_refinement_nodes;
        for (HypernodeID u = 0; u < phg.initialNumNodes(); ++u)
          if (context.refinement.fm.all_nodes || phg.isBorderNode(u))
            test_refinement_nodes.push_back(u);
        LocalizedKWayFM& fm = ets_fm.local();
        fm.findMoves(phg, sharedData, test_refinement_nodes);
      }

      FMStats stats;
      for (auto& fm : ets_fm) {
        fm.stats.merge(stats);
      }
      LOG << "Overall stats" << stats.serialize() << V(sharedData.moveTracker.numPerformedMoves());

      sharedData.refinementNodes.clear();  // calling clear is necessary since tryPop will reduce the size to -(num calling threads)

      timer.stop_timer("find_moves");
      timer.start_timer("rollback", "Rollback to Best Solution");

      HyperedgeWeight improvement = globalRollback.revertToBestPrefix(phg, sharedData, initialPartWeights, context.partition.max_part_weights[0]);
      overall_improvement += improvement;

      timer.stop_timer("rollback");

      if (improvement <= 0) {
        break;
      }
    }

    timer.stop_timer("fm");

    // sharedData.partition_weight_budgets.updatePartWeights(phg, context.partition.max_part_weights);  // only for version with budgets
    return overall_improvement;
  }

  void initialize(PartitionedHypergraph& phg) {
    // insert border nodes into work queues
    sharedData.refinementNodes.clear();
    tbb::parallel_for(HypernodeID(0), phg.initialNumNodes(), [&](const HypernodeID u) {
    //for (NodeID u = 0; u < phg.initialNumNodes(); ++u)
      if (phg.isBorderNode(u)) {
        sharedData.refinementNodes.push_back(u);
      }
    });

    // requesting new searches activates all nodes by raising the deactivated node marker
    // also clears the array tracking search IDs in case of overflow
    sharedData.nodeTracker.requestNewSearches(static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));

    // shuffle work queues if requested
    if (context.refinement.fm.shuffle) {
      sharedData.refinementNodes.shuffleQueue();
    }

  }
//protected:

  const Context& context;
  const TaskGroupID taskGroupID;
  FMSharedData sharedData;
  GlobalRollback globalRollback;
  tbb::enumerable_thread_specific<LocalizedKWayFM> ets_fm;
};

}
}