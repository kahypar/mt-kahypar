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

#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"

#include "mt-kahypar/utils/timer.h"
#include "kahypar/partition/metrics.h"
#include "mt-kahypar/utils/memory_tree.h"

namespace mt_kahypar {

  template<typename FMStrategy>
  bool MultiTryKWayFM<FMStrategy>::refineImpl(
              PartitionedHypergraph& phg,
              const vec<HypernodeID>& refinement_nodes,
              kahypar::Metrics& metrics,
              const double time_limit) {

    if (!is_initialized) throw std::runtime_error("Call initialize on fm before calling refine");

    Gain overall_improvement = 0;
    size_t consecutive_rounds_with_too_little_improvement = 0;
    enable_light_fm = false;
    sharedData.release_nodes = context.refinement.fm.release_nodes;
    sharedData.perform_moves_global = context.refinement.fm.perform_moves_global;
    double current_time_limit = time_limit;
    tbb::task_group tg;
    vec<HypernodeWeight> initialPartWeights(size_t(sharedData.numParts));
    HighResClockTimepoint fm_start = std::chrono::high_resolution_clock::now();
    utils::Timer& timer = utils::Timer::instance();

    for (size_t round = 0; round < context.refinement.fm.multitry_rounds; ++round) { // global multi try rounds
      for (PartitionID i = 0; i < sharedData.numParts; ++i) {
        initialPartWeights[i] = phg.partWeight(i);
      }

      timer.start_timer("collect_border_nodes", "Collect Border Nodes");
      roundInitialization(phg, refinement_nodes);
      timer.stop_timer("collect_border_nodes");

      size_t num_border_nodes = sharedData.refinementNodes.unsafe_size();
      if (num_border_nodes == 0) {
        break;
      }
      size_t num_seeds = context.refinement.fm.num_seed_nodes;
      if (context.type == kahypar::ContextType::main
          && !refinement_nodes.empty()  /* n-level */
          && num_border_nodes < 20 * context.shared_memory.num_threads) {
        num_seeds = num_border_nodes / (4 * context.shared_memory.num_threads);
        num_seeds = std::min(num_seeds, context.refinement.fm.num_seed_nodes);
        num_seeds = std::max(num_seeds, 1UL);
      }

      timer.start_timer("find_moves", "Find Moves");
      sharedData.finishedTasks.store(0, std::memory_order_relaxed);

      auto task = [&](const size_t task_id) {
        auto& fm = ets_fm.local();
        while(sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
              && fm.findMoves(phg, task_id, num_seeds)) { /* keep running*/ }
        sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
      };
      size_t num_tasks = std::min(num_border_nodes, context.shared_memory.num_threads);
      ASSERT(static_cast<int>(num_tasks) <= TBBNumaArena::instance().total_number_of_threads());
      for (size_t i = 0; i < num_tasks; ++i) {
        tg.run(std::bind(task, i));
      }
      tg.wait();
      timer.stop_timer("find_moves");

      timer.start_timer("rollback", "Rollback to Best Solution");
      HyperedgeWeight improvement = globalRollback.revertToBestPrefix
              <FMStrategy::maintain_gain_cache_between_rounds>(phg, sharedData, initialPartWeights);
      timer.stop_timer("rollback");

      const double roundImprovementFraction = improvementFraction(improvement, metrics.km1 - overall_improvement);
      overall_improvement += improvement;
      if (roundImprovementFraction < context.refinement.fm.min_improvement) {
        consecutive_rounds_with_too_little_improvement++;
      } else {
        consecutive_rounds_with_too_little_improvement = 0;
      }

      HighResClockTimepoint fm_timestamp = std::chrono::high_resolution_clock::now();
      const double elapsed_time = std::chrono::duration<double>(fm_timestamp - fm_start).count();
      if (debug && context.type == kahypar::ContextType::main) {
        FMStats stats;
        for (auto& fm : ets_fm) {
          fm.stats.merge(stats);
        }
        LOG << V(round) << V(improvement) << V(metrics::km1(phg)) << V(metrics::imbalance(phg, context))
            << V(num_border_nodes) << V(roundImprovementFraction) << V(elapsed_time) << V(current_time_limit)
            << stats.serialize();
      }

      // Enforce a time limit (based on k and coarsening time).
      // Switch to more "light-weight" FM after reaching it the first time. Abort after second time.
      if ( elapsed_time > current_time_limit ) {
        if ( !enable_light_fm ) {
          DBG << RED << "Multitry FM reached time limit => switch to Light FM Configuration" << END;
          sharedData.release_nodes = false;
          sharedData.perform_moves_global = true;
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

    #ifndef USE_STRONG_PARTITIONER
    is_initialized = false;
    #endif

    metrics.km1 -= overall_improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    ASSERT(metrics.km1 == metrics::km1(phg), V(metrics.km1) << V(metrics::km1(phg)));
    return overall_improvement > 0;
  }

  template<typename FMStrategy>
  void MultiTryKWayFM<FMStrategy>::roundInitialization(PartitionedHypergraph& phg,
                                                       const vec<HypernodeID>& refinement_nodes) {
    // clear border nodes
    sharedData.refinementNodes.clear();

    if ( refinement_nodes.empty() ) {
      // log(n) level case
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
    } else {
      // n-level case
      tbb::parallel_for(0UL, refinement_nodes.size(), [&](const size_t i) {
        const HypernodeID u = refinement_nodes[i];
        const int task_id = tbb::this_task_arena::current_thread_index();
        ASSERT(task_id >= 0 && task_id < TBBNumaArena::instance().total_number_of_threads());
        if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
          sharedData.refinementNodes.safe_push(u, task_id);
        }
      });
    }

    // shuffle task queue if requested
    if (context.refinement.fm.shuffle) {
      sharedData.refinementNodes.shuffle();
    }

    // requesting new searches activates all nodes by raising the deactivated node marker
    // also clears the array tracking search IDs in case of overflow
    sharedData.nodeTracker.requestNewSearches(static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));
  }


  template<typename FMStrategy>
  void MultiTryKWayFM<FMStrategy>::initializeImpl(PartitionedHypergraph& phg) {
    if ( !phg.isGainCacheInitialized() && FMStrategy::maintain_gain_cache_between_rounds ) {
      utils::Timer& timer = utils::Timer::instance();
      timer.start_timer("init_gain_info", "Initialize Gain Information");
      phg.initializeGainCache();
      timer.stop_timer("init_gain_info");
    }

    is_initialized = true;
  }

  template<typename FMStrategy>
  void MultiTryKWayFM<FMStrategy>::printMemoryConsumption() {
    utils::MemoryTreeNode fm_memory("Multitry k-Way FM", utils::OutputType::MEGABYTE);

    for (const auto& fm : ets_fm) {
      fm.memoryConsumption(&fm_memory);
    }
    sharedData.memoryConsumption(&fm_memory);
    fm_memory.finalize();

    LOG << BOLD << "\n FM Memory Consumption" << END;
    LOG << fm_memory;
  }

} // namespace mt_kahypar

#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"

namespace mt_kahypar {
  template class MultiTryKWayFM<GainCacheStrategy>;
  template class MultiTryKWayFM<GainDeltaStrategy>;
  template class MultiTryKWayFM<RecomputeGainStrategy>;
  template class MultiTryKWayFM<GainCacheOnDemandStrategy>;
}
