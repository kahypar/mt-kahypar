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

#include <mt-kahypar/parallel/numa_work_queue.h>
#include <mt-kahypar/partition/context.h>

#include <atomic>

#include "localized_kway_fm_core.h"

namespace mt_kahypar {
namespace refinement {

class MultiTryKWayFM {
public:
  MultiTryKWayFM(const Context& context, TaskGroupID taskGroupID, size_t numNodes) :
          context(context),
          taskGroupID(taskGroupID),
          globalMoveTracker(std::make_shared<GlobalMoveTracker>(numNodes)),
          ets_fm(context, globalMoveTracker)
  { }


  bool refine(PartitionedHypergraph& phg) {

    bool overall_improved = false;
    bool round_improved = true;

    // global multi try rounds
    for (size_t round = 0; round < context.refinement.fm.multitry_rounds && round_improved; ++round) {
      initialize(phg);

      auto task = [&](const int socket, const int socket_local_task_id, const int task_id) {
        HypernodeID u = std::numeric_limits<HypernodeID>::max();
        LocalizedKWayFM& fm = ets_fm.local();
        while (refinementNodes.tryPop(u, socket) /* && u not marked */ ) {
          fm.findMoves(phg, u);
        }
      };

      TBBNumaArena::instance().run_max_concurrency_tasks_on_all_sockets(taskGroupID, task);

      HyperedgeWeight improvement = globalRollbackToBestPrefix(phg);

      round_improved = true;//improvement > 0;
      overall_improved |= round_improved;
    }

    return overall_improved;
  }

  void initialize(PartitionedHypergraph& phg) {
    assert(refinementNodes.empty());

    phg.setRemainingOriginalPins();

    // insert border nodes into work queues
    tbb::parallel_for(HypernodeID(0), phg.initialNumNodes(), [&](const HypernodeID u) {
      if (phg.isBorderNode(u)) {
        refinementNodes.push(u, common::get_numa_node_of_vertex(u));
      }
    });

    // shuffle work queues if requested
    if (context.refinement.fm.shuffle) {
      refinementNodes.shuffleQueues();
    }

  }

  // TODO move to its own class
  HyperedgeWeight globalRollbackToBestPrefix(PartitionedHypergraph& phg) {

    GlobalMoveTracker move_tracker = *globalMoveTracker;
    
    // gain recalculation
    const MoveID numMoves = move_tracker.numPerformedMoves();
    tbb::parallel_for(0U, numMoves, [&](MoveID localMoveID) {
      MoveID moveID = move_tracker.firstMoveID + localMoveID;
      Gain gain = 0;
      Move& m = move_tracker.globalMoveOrder[localMoveID];
      const HypernodeID u = m.node;
      for (HyperedgeID e : phg.incidentEdges(u)) {
          if (phg.remainingPinsFromBeginningOfMovePhase(e, m.from) == 0  && phg.lastMoveOut(e, m.from) == moveID && phg.firstMoveIn(e, m.from) > moveID) {
            gain += phg.edgeWeight(e);
          }

          const MoveID lastOutTo = phg.lastMoveOut(e, m.to);
          const bool lastOutCameEarlierInThisRound = lastOutTo >= move_tracker.firstMoveID && lastOutTo < moveID;
          if (phg.remainingPinsFromBeginningOfMovePhase(e, m.to) == 0 && phg.firstMoveIn(e, m.to) == moveID && lastOutCameEarlierInThisRound) {
            gain -= phg.edgeWeight(e);
          }
      }
      m.gain = gain;
    });

    // one-pass prefix sum style to find the index with the best gain
    size_t num_tasks = context.shared_memory.num_threads;
    const std::vector<size_t> range_boundaries = parallel::Chunking::getChunkBorders(numMoves, num_tasks);
    std::vector<size_t> best_index_in_subrange(num_tasks, 0);
    std::vector<HyperedgeWeight> sum_at_beginning_of_subrange(num_tasks + 1, 0), best_sum_in_subrange(num_tasks, 0);

    // index i represents that all moves up to but excluding i are performed

    tbb::task_group tg;
    for (size_t task_id = 0; task_id < num_tasks; ++task_id) {
      tg.run([&, task_id]() {
        HyperedgeWeight sum = 0, best_sum = 0;
        size_t best_index = range_boundaries[task_id];
        for (size_t i = range_boundaries[task_id]; i < range_boundaries[task_id + 1]; ++i) {
          sum += move_tracker.globalMoveOrder[i].gain;
          if (sum > best_sum) {
            best_sum = sum;
            best_index = i + 1;
          }
        }

        sum_at_beginning_of_subrange[task_id + 1] = sum;
        best_sum_in_subrange[task_id] = best_sum;
        best_index_in_subrange [task_id] = best_index;
      });
    }
    tg.wait();

    size_t best_index = 0;
    HyperedgeWeight best_sum = 0, sum = 0;
    for (size_t task_id = 0; task_id < num_tasks; ++task_id) {
      sum += sum_at_beginning_of_subrange[task_id];

      if (sum + best_sum_in_subrange[task_id] > best_sum) {
        best_index = best_index_in_subrange[task_id];
        best_sum = sum + best_sum_in_subrange[task_id];
      }
    }

    // revert rejected moves
    tbb::parallel_for(best_index, size_t(numMoves), [&](const size_t moveID) {
      const Move& m = move_tracker.globalMoveOrder[moveID];
      phg.changeNodePart(m.node, m.to, m.from);
    });

    // reset stored move IDs either by raising the firstMoveID or by a hard reset if necessary
    if (move_tracker.reset()) {
      phg.resetStoredMoveIDs();
    }

    return best_sum;
  }

protected:

  const Context& context;
  const TaskGroupID taskGroupID;
  std::shared_ptr<GlobalMoveTracker> globalMoveTracker;
  NumaWorkQueue<HypernodeID> refinementNodes;
  tbb::enumerable_thread_specific<LocalizedKWayFM> ets_fm;
};

}
}