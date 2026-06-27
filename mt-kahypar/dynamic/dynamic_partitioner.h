#pragma once

#include "mt-kahypar/dynamic/dynamic_io.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance_vcycle.h"
#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/streaming.h"
#include "mt-kahypar/dynamic/strategies/hashing.h"
#include "mt-kahypar/dynamic/competitors/hermes.h"
#include "mt-kahypar/dynamic/competitors/leopard.h"
#include "mt-kahypar/dynamic/strategies/greedy.h"
#include "mt-kahypar/partition/registries/registry.h"

namespace mt_kahypar::dyn {
  inline int partition(Context& context) {

      context.dynamic.other_timings = OtherTimings();
      auto total_time_start = std::chrono::high_resolution_clock::now();


      context.partition.instance_type = InstanceType::hypergraph;
      context.partition.objective = Objective::km1;
      context.partition.gain_policy = GainPolicy::km1;

      // Read Hypergraph
      auto graph_time_parsing_start = std::chrono::high_resolution_clock::now();
      mt_kahypar_hypergraph_t hypergraph_t = io::readInputFile(
              context.partition.graph_filename, context.partition.preset_type,
              context.partition.instance_type, context.partition.file_format,
              context.preprocessing.stable_construction_of_incident_edges);
      auto& hypergraph_m = utils::cast<ds::MutableHypergraph>(hypergraph_t);
      context.dynamic.other_timings.io_time_graph_parsing = std::chrono::high_resolution_clock::now() -
        graph_time_parsing_start;

      auto setup_time_start = std::chrono::high_resolution_clock::now();

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
      } else if (context.dynamic.strategy == "hashing") {
        strategy = new Hashing(hypergraph_m, context);
      } else if (context.dynamic.strategy == "repartition") {
        strategy = new Repartition(hypergraph_m, context);
      } else if (context.dynamic.strategy == "hermes") {
          strategy = new Hermes(hypergraph_m, context);
      } else if (context.dynamic.strategy == "leopard") {
        strategy = new Leopard(hypergraph_m, context);
      } else if (context.dynamic.strategy == "greedy") {
        strategy = new Greedy(hypergraph_m, context);
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      LocalFMRound localFM_round = LocalFMRound();
      context.dynamic.local_fm_round = &localFM_round;

      context.dynamic.other_timings.setup_time = std::chrono::high_resolution_clock::now() - setup_time_start;
      auto finalization_time_start = std::chrono::high_resolution_clock::now();

      try {

        auto init_time_start = std::chrono::high_resolution_clock::now();

        std::cout << "Processing " << num_changes << " changes" << std::endl;

        auto& hypergraph_p =  strategy->init();

        std::cout << "Initial km1: " << metrics::quality(hypergraph_p, context) << ", imbalance: " << metrics::imbalance(hypergraph_p, context) << std::endl;

        size_t log_step_size = num_changes * context.dynamic.logging_step_size_pct;

        auto duration_sum = std::chrono::high_resolution_clock::duration::zero();
        auto change_io_duration_sum = std::chrono::high_resolution_clock::duration::zero();
        auto log_duration_sum = std::chrono::high_resolution_clock::duration::zero();

        context.dynamic.other_timings.initialization_time = std::chrono::high_resolution_clock::now() - init_time_start - context.dynamic.other_timings.initial_partitioning_time;

        for (size_t i = 0; i < num_changes; ++i) {
          HighResClockTimepoint start_io = std::chrono::high_resolution_clock::now();
          Change change;
          if (!context.dynamic.stream_changes) {
            change = changes[i];
          } else {
            change = parseChange(changes_file);
          }
          change_io_duration_sum += std::chrono::high_resolution_clock::now() - start_io;
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
            strategy->partition(change, num_changes);
          auto duration = std::chrono::high_resolution_clock::now() - start;
          duration_sum += duration;
          if (log_step_size != 0 && i % log_step_size != 0 && i != num_changes - 1) {
            continue;
          }
          HighResClockTimepoint start_log = std::chrono::high_resolution_clock::now();
          log_km1_live(i+1, num_changes, context, DynamicStrategy::getPartitionedHypergraphCopy(*strategy), duration_sum);
          log_duration_sum += std::chrono::high_resolution_clock::now() - start_log;
        }

        context.dynamic.other_timings.strategy_time = duration_sum;
        context.dynamic.other_timings.io_time_change_parsing = change_io_duration_sum;
        context.dynamic.other_timings.logging_time = log_duration_sum;

        finalization_time_start = std::chrono::high_resolution_clock::now();

        strategy->printAdditionalFinalStats();

        if (context.partition.write_partition_file) {
          PartitionerFacade::writePartitionFile(
            utils::partitioned_hg_cast(strategy->getPartitionedHypergraphCopy(*strategy)), context.partition.graph_partition_filename);
          std::cout << "Final partition written to file: " << context.partition.graph_partition_filename << std::endl;
        }

        utils::delete_hypergraph(hypergraph_t);
        } catch (std::exception& e) {
          std::cerr << "Error: " << e.what() << std::endl;
          generateErrorFile(context, strategy, e);
          exit(1);
        }

    context.dynamic.other_timings.finalization_time = std::chrono::high_resolution_clock::now() - finalization_time_start;
    context.dynamic.other_timings.total_time = std::chrono::high_resolution_clock::now() - total_time_start;
    if (context.dynamic.save_other_timings) {
        generateOtherTimingsFile(context);
      }

      return 0;
    }
}

