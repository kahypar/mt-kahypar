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

#include <external_tools/kahypar/kahypar/partition/metrics.h>

#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/global_rollback.h"
#include "mt-kahypar/utils/memory_tree.h"


namespace mt_kahypar {

class MultiTryKWayFM final : public IRefiner {

  static constexpr bool debug = false;

public:
  MultiTryKWayFM(const Hypergraph& hypergraph,
                 const Context& c,
                 const TaskGroupID taskGroupID) :
    initial_num_nodes(hypergraph.initialNumNodes()),
    original_context(c),
    context(c),
    taskGroupID(taskGroupID),
    sharedData(hypergraph.initialNumNodes(), context),
    globalRollback(hypergraph, context, context.partition.k),
    ets_fm([&] {
      return constructLocalizedKWayFMSearch();
    }) {
    if (context.refinement.fm.obey_minimal_parallelism) {
      sharedData.finishedTasksLimit = std::min(8UL, context.shared_memory.num_threads);
    }
  }

  bool refineImpl(PartitionedHypergraph& phg,
                  kahypar::Metrics& metrics,
                  const double time_limit) override final {
    Gain improvement = refine(phg, metrics, time_limit);
    metrics.km1 -= improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    ASSERT(metrics.km1 == metrics::km1(phg), V(metrics.km1) << V(metrics::km1(phg)));
    return improvement > 0;
  }

  Gain refine(PartitionedHypergraph& phg,
              const kahypar::Metrics& metrics,
              const double time_limit) {
    if (!is_initialized) throw std::runtime_error("Call initialize on fm before calling refine");
    utils::Timer& timer = utils::Timer::instance();
    Gain overall_improvement = 0;
    size_t consecutive_rounds_with_too_little_improvement = 0;
    context = original_context;
    enable_light_fm = false;
    double current_time_limit = time_limit;
    HighResClockTimepoint fm_start = std::chrono::high_resolution_clock::now();
    for (size_t round = 0; round < context.refinement.fm.multitry_rounds; ++round) { // global multi try rounds
      timer.start_timer("collect_border_nodes", "Collect Border Nodes");

      roundInitialization(phg);
      size_t numBorderNodes = sharedData.refinementNodes.unsafe_size(); unused(numBorderNodes);

      timer.stop_timer("collect_border_nodes");
      timer.start_timer("find_moves", "Find Moves");

      vec<HypernodeWeight> initialPartWeights(size_t(sharedData.numParts));
      for (PartitionID i = 0; i < sharedData.numParts; ++i) initialPartWeights[i] = phg.partWeight(i);

      if (context.refinement.fm.algorithm == FMAlgorithm::fm_multitry) {
        sharedData.finishedTasks.store(0, std::memory_order_relaxed);
        auto task = [&](const int , const int task_id, const int ) {
          LocalizedKWayFM& fm = ets_fm.local();
          while(sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
                && fm.findMovesLocalized(phg, sharedData, static_cast<size_t>(task_id))) {
            /* keep running */
          }
          sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
        };
        TBBNumaArena::instance().execute_task_on_each_thread(taskGroupID, task);
      } else if (context.refinement.fm.algorithm == FMAlgorithm::fm_boundary){
        LocalizedKWayFM& fm = ets_fm.local();
        fm.findMovesUsingFullBoundary(phg, sharedData);
      }

      FMStats stats;
      for (auto& fm : ets_fm) {
        fm.stats.merge(stats);
      }
      peak_reinsertions = std::max(peak_reinsertions, stats.task_queue_reinsertions);

      timer.stop_timer("find_moves");
      timer.start_timer("rollback", "Rollback to Best Solution");

      HyperedgeWeight improvement = globalRollback.revertToBestPrefix(phg, sharedData, initialPartWeights);
      const double roundImprovementFraction = improvementFraction(improvement, metrics.km1 - overall_improvement);
      overall_improvement += improvement;

      timer.stop_timer("rollback");

      if (roundImprovementFraction < context.refinement.fm.min_improvement) {
        consecutive_rounds_with_too_little_improvement++;
      } else {
        consecutive_rounds_with_too_little_improvement = 0;
      }

      HighResClockTimepoint fm_timestamp = std::chrono::high_resolution_clock::now();
      const double elapsed_time = std::chrono::duration<double>(fm_timestamp - fm_start).count();
      if (debug && context.type == kahypar::ContextType::main) {
        LOG << V(round) << V(improvement) << V(metrics::km1(phg)) << V(metrics::imbalance(phg, context))
            << V(numBorderNodes) << V(roundImprovementFraction) << V(elapsed_time) << V(current_time_limit)
            << stats.serialize();
      }

      // We enforce during FM a time limit, which is calculated based on k and the coarsening time of
      // the current hypergraph. Especially for instances with low density, we observed that FM time
      // dominates. The root cause for this is that for such instances (with large hyperedges) maintaining
      // the gain cache inside the partitioned hypergraph becomes expensive.
      // After the FM reaches the time limit the first time, we switch to a light FM variant. Here, we do not
      // release vertices that are within a local PQ at the end of a localized search. If the time limit is
      // reached a second time we immediatly abort FM refinement.
      if ( elapsed_time > current_time_limit ) {
        if ( !enable_light_fm ) {
          DBG << RED << "Multitry FM reached time limit => switch to Light FM Configuration" << END;
          context.refinement.fm.release_nodes = false;
          current_time_limit *= 2;
          enable_light_fm = true;
        } else {
          DBG << RED << "Light version of Multitry FM reached time limit => ABORT" << END;
          break;
        }
      }

      if (improvement <= 0 || consecutive_rounds_with_too_little_improvement >= 2) {
        break;
      }
    }

    if (context.partition.show_memory_consumption && context.partition.verbose_output
        && context.type == kahypar::ContextType::main
        && phg.initialNumNodes() == sharedData.moveTracker.moveOrder.size() /* top level */) {
      printMemoryConsumption();
    }

    is_initialized = false;
    return overall_improvement;
  }

  void initializeImpl(PartitionedHypergraph& phg) override final {
    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("init_gain_info", "Initialize Gain Information");
    phg.initializeGainInformation();
    timer.stop_timer("init_gain_info");
    timer.start_timer("set_remaining_original_pins", "Set remaining original pins");
    globalRollback.setRemainingOriginalPins(phg);
    timer.stop_timer("set_remaining_original_pins");

    is_initialized = true;
  }

  void roundInitialization(PartitionedHypergraph& phg) {
    // clear border nodes
    sharedData.refinementNodes.clear();

    // iterate over all nodes and insert border nodes into task queue
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
      [&](const tbb::blocked_range<HypernodeID>& r) {
      const int task_id = tbb::this_task_arena::current_thread_index();
      ASSERT(task_id >= 0 && task_id < TBBNumaArena::instance().total_number_of_threads());
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
          sharedData.refinementNodes.safe_push(u, task_id);
        }
      }
    });

    // shuffle task queue if requested
    if (context.refinement.fm.shuffle) {
      sharedData.refinementNodes.shuffle();
    }

    // requesting new searches activates all nodes by raising the deactivated node marker
    // also clears the array tracking search IDs in case of overflow
    sharedData.nodeTracker.requestNewSearches(static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));

    sharedData.fruitlessSeed.reset();
  }

  LocalizedKWayFM constructLocalizedKWayFMSearch() {
    return LocalizedKWayFM(context, initial_num_nodes, sharedData.vertexPQHandles.data());
  }

  bool is_initialized = false;
  bool enable_light_fm = false;
  const HypernodeID initial_num_nodes;
  const Context& original_context;
  Context context;
  const TaskGroupID taskGroupID;
  FMSharedData sharedData;
  GlobalRollback globalRollback;
  tbb::enumerable_thread_specific<LocalizedKWayFM> ets_fm;
  size_t peak_reinsertions = 0;

  double improvementFraction(Gain gain, HyperedgeWeight old_km1) {
    if (old_km1 == 0)
      return 0;
    else
      return static_cast<double>(gain) / static_cast<double>(old_km1);
  }

  bool shouldStopSearch(const vec<double>& improvement_fractions, double threshold, size_t n) const {
    if (improvement_fractions.size() < n || context.type != kahypar::ContextType::main) {
      return false;
    } else {
      bool all_below = true;
      for (size_t i = improvement_fractions.size() - n; i < improvement_fractions.size(); ++i) {
        all_below &= (improvement_fractions[i] < threshold);
      }
      return all_below;
    }
  }

  void printMemoryConsumption() {
    utils::MemoryTreeNode fm_memory("Multitry k-Way FM", utils::OutputType::MEGABYTE);

    // Localized k-Way FM Memory Consumption
    for (const LocalizedKWayFM& fm : ets_fm) {
      fm.memoryConsumption(&fm_memory);
    }

    // Global Rollback
    globalRollback.memoryConsumption(&fm_memory);

    // Shared Data
    sharedData.memoryConsumption(&fm_memory);

    fm_memory.finalize();
    LOG << BOLD << "\n FM Memory Consumption" << END;
    LOG << fm_memory;
  }
};

}
