#pragma once

#include "mt-kahypar/dynamic/dynamic_io.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance_vcycle.h"
#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/streaming.h"
#include "mt-kahypar/partition/registries/registry.h"

namespace mt_kahypar::dyn {
  inline int partition(Context& context) {


      context.partition.instance_type = InstanceType::hypergraph;
      context.partition.objective = Objective::km1;
      context.partition.gain_policy = GainPolicy::km1;

      // Read Hypergraph
      mt_kahypar_hypergraph_t hypergraph_t = io::readInputFile(
              context.partition.graph_filename, context.partition.preset_type,
              context.partition.instance_type, context.partition.file_format,
              context.preprocessing.stable_construction_of_incident_edges);
      auto& hypergraph_m = utils::cast<ds::MutableHypergraph>(hypergraph_t);

      // Initialize Memory Pool and Algorithm/Policy Registries
      register_memory_pool(hypergraph_t, context);
      register_algorithms_and_policies();

      // Parse changes
      std::ifstream changes_file(context.dynamic.changes_file);
      if (!changes_file.is_open()) {
        throw std::runtime_error("Could not open file: " + context.dynamic.changes_file);
      }

      size_t num_changes = getNumChanges(changes_file);
      std::vector<Change> changes;
      if (!context.dynamic.stream_changes) {
          changes = parseChanges(std::move(changes_file), num_changes);
          // Process and delete setup changes
          DynamicStrategy::process_setup_changes(hypergraph_m, context, changes);
          changes.erase(changes.begin(), changes.begin() + context.dynamic.setup_moves_count);
      } else {
          // Process setup changes
          for (size_t i = 0; i < context.dynamic.setup_moves_count; ++i) {
            Change change = parseChange(changes_file);
            DynamicStrategy::process_change(hypergraph_m, context, change);
          }
      }

      mt_kahypar::dyn::DynamicStrategy* strategy;

      if (context.dynamic.strategy == "localFM_rebalance_vcycle") {
        strategy = new LocalFMRebalanceVCycleV4(hypergraph_m, context);
      } else if (context.dynamic.strategy == "streaming") {
        strategy = new Streaming(hypergraph_m, context);
      } else if (context.dynamic.strategy == "repartition") {
        strategy = new Repartition(hypergraph_m, context);
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      LocalFMRound localFM_round = LocalFMRound();
      context.dynamic.local_fm_round = &localFM_round;

      try {

        std::cout << "Processing " << num_changes << " changes" << std::endl;

        auto& hypergraph_p =  strategy->init();

        std::cout << "Initial km1: " << metrics::quality(hypergraph_p, context) << ", imbalance: " << metrics::imbalance(hypergraph_p, context) << std::endl;

        size_t log_step_size = num_changes * context.dynamic.logging_step_size_pct;

        auto duration_sum = std::chrono::high_resolution_clock::duration::zero();

        for (size_t i = 0; i < num_changes; ++i) {
          Change change;
          if (!context.dynamic.stream_changes) {
            change = changes[i];
          } else {
            change = parseChange(changes_file);
          }
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
            strategy->partition(change, num_changes);
          auto duration = std::chrono::high_resolution_clock::now() - start;
          duration_sum += duration;
          if (log_step_size != 0 && i % log_step_size != 0) {
            continue;
          }
          log_km1_live(i+1, num_changes, context, DynamicStrategy::getPartitionedHypergraphCopy(*strategy), duration_sum);
        }

        strategy->printAdditionalFinalStats();

        utils::delete_hypergraph(hypergraph_t);

        } catch (std::exception& e) {
          std::cerr << "Error: " << e.what() << std::endl;
          generateErrorFile(context, strategy, e);
          exit(1);
        }

      return 0;
    }
}

