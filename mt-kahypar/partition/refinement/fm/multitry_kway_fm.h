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


namespace mt_kahypar {

class MultiTryKWayFM final : public IRefiner {

  static constexpr bool debug = true;

public:
  MultiTryKWayFM(const Hypergraph& hypergraph, const Context& context, const TaskGroupID taskGroupID) :
          context(context),
          taskGroupID(taskGroupID),
          sharedData(hypergraph.initialNumNodes(), context),
          globalRollback(hypergraph, context, context.partition.k),
          ets_fm(context, hypergraph.initialNumNodes(), sharedData.vertexPQHandles.data())
  {
    if (context.refinement.fm.obey_minimal_parallelism) {
      sharedData.finishedTasksLimit = std::min(8UL, context.shared_memory.num_threads);
    }
  }

  bool refineImpl(PartitionedHypergraph& phg,
                  kahypar::Metrics& metrics) override final {
    if (shouldStopSearchAltogether()) {
      return false;
    }
    Gain improvement = refine(phg, metrics);
    levelImprovementFractions.push_back( improvementFraction(improvement, metrics.km1) );
    metrics.km1 -= improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    ASSERT(metrics.km1 == metrics::km1(phg), V(metrics.km1) << V(metrics::km1(phg)));
    return improvement > 0;
  }

  Gain refine(PartitionedHypergraph& phg, const kahypar::Metrics& metrics) {
    if (!is_initialized) throw std::runtime_error("Call initialize on fm before calling refine");
    utils::Timer& timer = utils::Timer::instance();
    Gain overall_improvement = 0;
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
      roundImprovementFractions.push_back( improvementFraction(improvement, metrics.km1 - overall_improvement) );
      overall_improvement += improvement;

      timer.stop_timer("rollback");

      if (debug && context.type == kahypar::ContextType::main) {
        LOG << V(round) << V(improvement) << V(metrics::km1(phg)) << V(metrics::imbalance(phg, context))
            << V(numBorderNodes) << V(roundImprovementFractions.back()) << stats.serialize();
      }

      if (improvement <= 0 || shouldStopSearchOnThisLevel()) {
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
    if (shouldStopSearchAltogether()) {
      return;
    }
    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("init_gain_info", "Initialize Gain Information");
    phg.initializeGainInformation();
    timer.stop_timer("init_gain_info");
    timer.start_timer("set_remaining_original_pins", "Set remaining original pins");
    globalRollback.setRemainingOriginalPins(phg);
    timer.stop_timer("set_remaining_original_pins");

    // clear gain tracking for the next level
    roundImprovementFractions.clear();

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
    sharedData.nodeTracker.requestNewSearches(
      static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));

    sharedData.fruitlessSeed.reset();
  }

  bool is_initialized = false;
  const Context& context;
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
    if (context.type == kahypar::ContextType::main)
	    LOG << V(threshold) << V(n) << V(improvement_fractions.size());
    if (improvement_fractions.size() < n || context.type != kahypar::ContextType::main) {
      return false;
    } else {
      bool all_below = true;
      for (size_t i = improvement_fractions.size() - n; i < improvement_fractions.size(); ++i) {
        all_below &= (improvement_fractions[i] < threshold);
      }
      LOG << V(threshold) << V(all_below);
      return all_below;
    }
  }

  bool shouldStopSearchAltogether() const {
    static constexpr size_t levels_to_consider = 3;
    double level_improvement_fraction_threshold = context.refinement.fm.min_improvement;
    return shouldStopSearch(levelImprovementFractions, level_improvement_fraction_threshold, levels_to_consider);
  }

  bool shouldStopSearchOnThisLevel() const {
    static constexpr size_t rounds_to_consider = 2;
    double round_improvement_fraction_threshold = context.refinement.fm.min_improvement;
    return shouldStopSearch(roundImprovementFractions, round_improvement_fraction_threshold, rounds_to_consider);
  }



  vec<double> roundImprovementFractions;
  vec<double> levelImprovementFractions;

  void printMemoryConsumption() {
    std::unordered_map<std::string, size_t> r;
    r["global rollback"] = globalRollback.memory_consumption();
    for (const LocalizedKWayFM& fm : ets_fm) {
      auto local_mem = fm.memory_consumption();
      for (const auto& it : local_mem) {
        r[it.first] += it.second;
      }
    }
    r["tbb concurrent task queue >="] = peak_reinsertions * sizeof(HypernodeID);
    for (const auto& it : sharedData.memory_consumption()) {
      r[it.first] += it.second;
    }
    LOG << "---------------------";
    LOG << "FM Memory Consumption";
    LOG << "---------------------";

    size_t left_width = 35; size_t right_width = 10;
    for (const auto& it : r) {
      std::string right_word = std::to_string(it.second / (1024*1024));
      size_t pad = (left_width - it.first.length()) + (right_width - right_word.length());
      std::cout << it.first;
      for (size_t i = 0; i < pad; ++i) std::cout << " ";
      std::cout << right_word << "MB\n";
    }
    LOG << "---------------------";
  }
};

}
