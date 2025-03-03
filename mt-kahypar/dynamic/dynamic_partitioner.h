#pragma once

#include <mt-kahypar/dynamic/strategies/localFM_factor.h>
#include <mt-kahypar/dynamic/strategies/rebalance.h>
#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/connectivity.h"
#include "mt-kahypar/dynamic/strategies/localFM.h"
#include "mt-kahypar/dynamic/strategies/localFM_filter.h"
#include "mt-kahypar/dynamic/strategies/localFM_old.h"
#include "mt-kahypar/dynamic/strategies/localFM_incremental_gain.h"
#include "mt-kahypar/dynamic/strategies/localFM_filtered_gain.h"
#include "mt-kahypar/dynamic/strategies/localFM_v2.h"
#include "mt-kahypar/dynamic/strategies/localFM_small_blocks.h"
#include "mt-kahypar/dynamic/strategies/never_repartition.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance_vcycle.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance_debug.h"
#include "mt-kahypar/dynamic/strategies/v_cycle.h"
#include "mt-kahypar/dynamic/dynamic_io.h"

namespace mt_kahypar::dyn {

    void partition(Context& context) {

      context.partition.instance_type = InstanceType::hypergraph;
      context.partition.objective = Objective::km1;
      context.partition.gain_policy = GainPolicy::km1;

      // Read Hypergraph
      mt_kahypar_hypergraph_t hypergraph_t = io::readInputFile(
              context.partition.graph_filename, context.partition.preset_type,
              context.partition.instance_type, context.partition.file_format,
              context.preprocessing.stable_construction_of_incident_edges);
      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph_t);

      // Parse or generate changes
      std::vector<Change> changes;
      if (context.dynamic.changes_file.empty()) {
        changes = generateChanges(hypergraph_s, context);
      } else {
        changes = parseChanges(context.dynamic.changes_file);
        //TODO add reset flag
        //resetHypergraph(hypergraph_s, changes, context);
        for (size_t i = 0; i < context.dynamic.setup_moves_count; ++i) {
          Change& change = changes[i];
          DynamicStrategy::process_change(hypergraph_s, context, change);
        }
      }

      std::cout << "Number of changes: " << changes.size() << std::endl;

      // If the max_changes is not specified or is greater than the number of changes in the file, we process all the changes
      size_t max_changes = context.dynamic.max_changes == 0 ? changes.size() : std::min((size_t) context.dynamic.max_changes, changes.size());

      DynamicStrategy* strategy;

      if (context.dynamic.strategy == "connectivity") {
        strategy = new Connectivity();
      } else if (context.dynamic.strategy == "repartition") {
        strategy = new Repartition();
      } else if (context.dynamic.strategy == "localFM") {
        strategy = new LocalFM();
      } else if (context.dynamic.strategy == "localFM_filter") {
        strategy = new LocalFMFilter();
      } else if (context.dynamic.strategy == "localFM_old") {
        strategy = new LocalFMOld();
      } else if (context.dynamic.strategy == "never_repartition") {
        strategy = new NeverRepartition();
      } else if (context.dynamic.strategy == "localFM_factor") {
        strategy = new LocalFMFactor();
      } else if (context.dynamic.strategy == "localFM_incremental_gain") {
        strategy = new LocalFMIncGain();
      } else if (context.dynamic.strategy == "localFM_filtered_gain") {
        strategy = new LocalFMFilterGain();
      } else if (context.dynamic.strategy == "localFM_v2") {
        strategy = new LocalFMV2();
      } else if (context.dynamic.strategy == "localFM_small_blocks") {
        strategy = new LocalFMSBlocks();
      } else if (context.dynamic.strategy == "localFM_rebalance") {
        strategy = new LocalFMRebalance();
      } else if (context.dynamic.strategy == "localFM_rebalance_vcycle") {
        strategy = new LocalFMRebalanceVCycle();
      }else if (context.dynamic.strategy == "localFM_rebalance_debug") {
        strategy = new LocalFMRebalanceDebug();
      } else if (context.dynamic.strategy == "rebalance") {
        strategy = new Rebalance();
      } else if (context.dynamic.strategy == "v_cycle") {
        strategy = new VCycle();
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      mt_kahypar::LocalFMRound localFM_round = mt_kahypar::LocalFMRound();
      localFM_round.overall_improvement = 0;
      localFM_round.touched_nodes = 0;
      localFM_round.moved_nodes = 0;
      context.dynamic.localFM_round = &localFM_round;

      try {

        std::cout << "Processing " << max_changes << " changes" << std::endl;

        strategy->init(hypergraph_s, context);

        size_t log_step_size = max_changes * context.dynamic.logging_step_size_pct;

        auto duration_sum = std::chrono::high_resolution_clock::duration::zero();

        for (size_t i = context.dynamic.setup_moves_count; i < max_changes; ++i) {
          Change& change = changes[i];
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
          strategy->partition(hypergraph_s, context, change, max_changes);
          auto duration = std::chrono::high_resolution_clock::now() - start;
          duration_sum += duration;
          if (log_step_size != 0 && i % log_step_size != 0) {
            continue;
          }
          strategy->compute_km1_and_imbalance(hypergraph_s, context, change, strategy->history.back());
          log_km1_live(i, context, strategy->history.back(), duration_sum);
          if (!context.dynamic.server && *(&strategy->history.back().valid)) {
            print_progress_bar(i, max_changes, &strategy->history);
          }
        }

        strategy->printFinalStats(hypergraph_s, context);
        //log_km1(context, &strategy->history);

        utils::delete_hypergraph(hypergraph_t);

        } catch (std::exception& e) {
          std::cerr << "Error: " << e.what() << std::endl;
          generateErrorFile(context, strategy, e);
          exit(1);
        }
    }
}

