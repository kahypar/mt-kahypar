/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

  template<typename TypeTraits, typename GainTypes>
  bool MultiTryKWayFM<TypeTraits, GainTypes>::refineImpl(
              mt_kahypar_partitioned_hypergraph_t& hypergraph,
              const vec<HypernodeID>& refinement_nodes,
              Metrics& metrics,
              const double time_limit) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    if (!is_initialized) throw std::runtime_error("Call initialize on fm before calling refine");
    resizeDataStructuresForCurrentK();

    Gain overall_improvement = 0;
    size_t consecutive_rounds_with_too_little_improvement = 0;
    enable_light_fm = false;
    sharedData.release_nodes = context.refinement.fm.release_nodes;
    sharedData.perform_moves_global = context.refinement.fm.perform_moves_global;
    double current_time_limit = time_limit;
    tbb::task_group tg;
    vec<HypernodeWeight> initialPartWeights(size_t(context.partition.k));
    HighResClockTimepoint fm_start = std::chrono::high_resolution_clock::now();
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);

    for (size_t round = 0; round < context.refinement.fm.multitry_rounds; ++round) { // global multi try rounds
      for (PartitionID i = 0; i < context.partition.k; ++i) {
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
      if (context.type == ContextType::main
          && !refinement_nodes.empty()  /* n-level */
          && num_border_nodes < 20 * context.shared_memory.num_threads) {
        num_seeds = num_border_nodes / (4 * context.shared_memory.num_threads);
        num_seeds = std::min(num_seeds, context.refinement.fm.num_seed_nodes);
        num_seeds = std::max(num_seeds, UL(1));
      }

      timer.start_timer("find_moves", "Find Moves");
      sharedData.finishedTasks.store(0, std::memory_order_relaxed);

      auto task = [&](const size_t task_id) {
        auto& fm = ets_fm.local();
        while(sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
              && fm.findMoves(phg, task_id, num_seeds)) { /* keep running*/ }
        sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
      };
      size_t num_tasks = std::min(num_border_nodes, size_t(TBBInitializer::instance().total_number_of_threads()));
      for (size_t i = 0; i < num_tasks; ++i) {
        tg.run(std::bind(task, i));
      }
      tg.wait();
      timer.stop_timer("find_moves");

      timer.start_timer("rollback", "Rollback to Best Solution");
      HyperedgeWeight improvement = globalRollback.revertToBestPrefix(
        phg, sharedData, initialPartWeights);
      timer.stop_timer("rollback");

      const double roundImprovementFraction = improvementFraction(improvement,
        metrics.quality - overall_improvement);
      overall_improvement += improvement;
      if (roundImprovementFraction < context.refinement.fm.min_improvement) {
        consecutive_rounds_with_too_little_improvement++;
      } else {
        consecutive_rounds_with_too_little_improvement = 0;
      }

      HighResClockTimepoint fm_timestamp = std::chrono::high_resolution_clock::now();
      const double elapsed_time = std::chrono::duration<double>(fm_timestamp - fm_start).count();
      if (debug && context.type == ContextType::main) {
        FMStats stats;
        for (auto& fm : ets_fm) {
          fm.stats.merge(stats);
        }
        LOG << V(round) << V(improvement) << V(metrics::quality(phg, context))
            << V(metrics::imbalance(phg, context)) << V(num_border_nodes) << V(roundImprovementFraction)
            << V(elapsed_time) << V(current_time_limit) << stats.serialize();
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
        && context.type == ContextType::main
        && phg.initialNumNodes() == sharedData.moveTracker.moveOrder.size() /* top level */) {
      printMemoryConsumption();
    }

    if ( !context.isNLevelPartitioning() ) {
      is_initialized = false;
    }

    metrics.quality -= overall_improvement;
    metrics.imbalance = metrics::imbalance(phg, context);
    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(gain_cache));
    ASSERT(metrics.quality == metrics::quality(phg, context),
      V(metrics.quality) << V(metrics::quality(phg, context)));
    return overall_improvement > 0;
  }

  template<typename TypeTraits, typename GainTypes>
  void MultiTryKWayFM<TypeTraits, GainTypes>::roundInitialization(PartitionedHypergraph& phg,
                                                                  const vec<HypernodeID>& refinement_nodes) {
    // clear border nodes
    sharedData.refinementNodes.clear();

    if ( refinement_nodes.empty() ) {
      // log(n) level case
      // iterate over all nodes and insert border nodes into task queue
      tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
        [&](const tbb::blocked_range<HypernodeID>& r) {
          const int task_id = tbb::this_task_arena::current_thread_index();
          // In really rare cases, the tbb::this_task_arena::current_thread_index()
          // function a thread id greater than max_concurrency which causes an
          // segmentation fault if we do not perform the check here. This is caused by
          // our working queue for border nodes with which we initialize the localized
          // FM searches. For now, we do not know why this occurs but this prevents
          // the segmentation fault.
          if ( task_id >= 0 && task_id < TBBInitializer::instance().total_number_of_threads() ) {
            for (HypernodeID u = r.begin(); u < r.end(); ++u) {
              if (phg.nodeIsEnabled(u) && phg.isBorderNode(u) && !phg.isFixed(u)) {
                sharedData.refinementNodes.safe_push(u, task_id);
              }
            }
          }
        });
    } else {
      // n-level case
      tbb::parallel_for(UL(0), refinement_nodes.size(), [&](const size_t i) {
        const HypernodeID u = refinement_nodes[i];
        const int task_id = tbb::this_task_arena::current_thread_index();
        if ( task_id >= 0 && task_id < TBBInitializer::instance().total_number_of_threads() ) {
          if (phg.nodeIsEnabled(u) && phg.isBorderNode(u) && !phg.isFixed(u)) {
            sharedData.refinementNodes.safe_push(u, task_id);
          }
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


  template<typename TypeTraits, typename GainTypes>
  void MultiTryKWayFM<TypeTraits, GainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    if (!gain_cache.isInitialized()) {
      gain_cache.initializeGainCache(phg);
    }

    is_initialized = true;
  }

  template<typename TypeTraits, typename GainTypes>
  void MultiTryKWayFM<TypeTraits, GainTypes>::resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if ( current_k != context.partition.k ) {
      current_k = context.partition.k;
      // Note that in general changing the number of blocks in the
      // global rollback data structure should not resize any data structure
      // as we initialize them with the final number of blocks. This is just a fallback
      // if someone changes this in the future.
      globalRollback.changeNumberOfBlocks(current_k);
      for ( auto& localized_fm : ets_fm ) {
        localized_fm.changeNumberOfBlocks(current_k);
      }
      gain_cache.changeNumberOfBlocks(current_k);
    }
  }

  template<typename TypeTraits, typename GainTypes>
  void MultiTryKWayFM<TypeTraits, GainTypes>::printMemoryConsumption() {
    utils::MemoryTreeNode fm_memory("Multitry k-Way FM", utils::OutputType::MEGABYTE);

    for (const auto& fm : ets_fm) {
      fm.memoryConsumption(&fm_memory);
    }
    sharedData.memoryConsumption(&fm_memory);
    fm_memory.finalize();

    LOG << BOLD << "\n FM Memory Consumption" << END;
    LOG << fm_memory;
  }

  namespace {
  #define MULTITRY_KWAY_FM(X, Y) MultiTryKWayFM<X, Y>
  }

  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(MULTITRY_KWAY_FM)

} // namespace mt_kahypar
