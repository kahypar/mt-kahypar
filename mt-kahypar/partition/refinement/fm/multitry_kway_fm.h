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

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/global_rollback.h"


namespace mt_kahypar {

class MultiTryKWayFM final : public IRefiner {

  static constexpr bool debug = false;

public:
  MultiTryKWayFM(Hypergraph& hypergraph, const Context& context, const TaskGroupID taskGroupID) :
          context(context),
          taskGroupID(taskGroupID),
          sharedData(hypergraph.initialNumNodes(), context),
          globalRollback(hypergraph.initialNumNodes(), hypergraph.initialNumEdges(), context.partition.k),
          ets_fm(context, hypergraph.initialNumNodes(), sharedData.vertexPQHandles.data()) { }

  bool refineImpl(PartitionedHypergraph& phg,
                  kahypar::Metrics& metrics) override final {
    Gain improvement = refine(phg);
    metrics.km1 -= improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    assert(metrics.km1 == metrics::km1(phg));
    return improvement > 0;
  }

  Gain refine(PartitionedHypergraph& phg) {
    if (!is_initialized) throw std::runtime_error("Call initialize on fm before calling refine");

    utils::Timer& timer = utils::Timer::instance();
    Gain overall_improvement = 0;
    for (size_t round = 0; round < context.refinement.fm.multitry_rounds; ++round) {                    // global multi try rounds
      timer.start_timer("collect_border_nodes", "Collect Border Nodes");

      roundInitialization(phg);

      timer.stop_timer("collect_border_nodes");
      timer.start_timer("find_moves", "Find Moves");

      vec<HypernodeWeight> initialPartWeights(size_t(sharedData.numParts));
      for (PartitionID i = 0; i < sharedData.numParts; ++i) initialPartWeights[i] = phg.partWeight(i);

      if (context.refinement.fm.algorithm == FMAlgorithm::fm_multitry) {
        auto task = [&](const int , const int , const int ) {
          LocalizedKWayFM& fm = ets_fm.local();
          while(fm.findMoves(phg, sharedData)) { /* keep running */ }
        };
        TBBNumaArena::instance().execute_task_on_each_thread(taskGroupID, task);
        //task(0,0,0);
      } else if (context.refinement.fm.algorithm == FMAlgorithm::fm_boundary){
        // Try boundary FM
        vec<HypernodeID> test_refinement_nodes;
        for (HypernodeID u = 0; u < phg.initialNumNodes(); ++u)
          if (context.refinement.fm.init_boundary_fm_with_all_nodes || phg.isBorderNode(u))
            test_refinement_nodes.push_back(u);
        LocalizedKWayFM& fm = ets_fm.local();
        fm.findMoves(phg, sharedData, test_refinement_nodes);
      }

      FMStats stats;
      for (auto& fm : ets_fm) {
        fm.stats.merge(stats);
      }
      DBG << stats.serialize();

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

    is_initialized = false;
    return overall_improvement;
  }

  void initializeImpl(PartitionedHypergraph& phg) override final {
    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("init_gain_info", "Initialize Gain Information");
    // initialization only as long as LP refiner does not use these datastructures TODO consolidate at some point
    phg.initializeGainInformation();
    timer.stop_timer("init_gain_info");
    timer.start_timer("set_remaining_original_pins", "Set remaining original pins");
    // initialization only as long as LP refiner does not use these datastructures
    globalRollback.setRemainingOriginalPins(phg);
    timer.stop_timer("set_remaining_original_pins");
    is_initialized = true;
  }

  void roundInitialization(PartitionedHypergraph& phg) {
    // insert border nodes into work queues
    sharedData.refinementNodes.clear();
    ds::StreamingVector<HypernodeID> tmpRefinementNodes;
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if (phg.isBorderNode(hn)) {
        tmpRefinementNodes.stream(hn);
      }
    });
    sharedData.refinementNodes.size.store(tmpRefinementNodes.size());
    tmpRefinementNodes.copy_parallel(sharedData.refinementNodes.elements);

    // requesting new searches activates all nodes by raising the deactivated node marker
    // also clears the array tracking search IDs in case of overflow
    sharedData.nodeTracker.requestNewSearches(static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));

    // shuffle work queues if requested
    if (context.refinement.fm.shuffle) {
      sharedData.refinementNodes.shuffle();
    }
  }

  bool is_initialized = false;
  const Context& context;
  const TaskGroupID taskGroupID;
  FMSharedData sharedData;
  GlobalRollback globalRollback;
  tbb::enumerable_thread_specific<LocalizedKWayFM> ets_fm;
};

}