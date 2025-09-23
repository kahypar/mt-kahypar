#pragma once

#include "mt-kahypar/dynamic/dynamic_io.h"
#include "mt-kahypar/dynamic/strategies/localFM_rebalance_vcycle.h"
#include "mt-kahypar/partition/registries/registry.h"

namespace mt_kahypar::dyn {

    int partition(Context& context) {


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

      // Parse or generate changes
      std::vector<Change> changes = parseChanges(context.dynamic.changes_file);

      std::cout << "Number of changes: " << changes.size() << std::endl;

      // If the max_changes is not specified or is greater than the number of changes in the file, we process all the changes
      size_t max_changes = context.dynamic.max_changes == 0 ? changes.size() : std::min((size_t) context.dynamic.max_changes, changes.size());

      mt_kahypar::dyn::DynamicStrategy* strategy;

      if (context.dynamic.strategy == "localFM_rebalance_vcycle") {
        strategy = new LocalFMRebalanceVCycleV4(hypergraph_m, context);
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      try {

        std::cout << "Processing " << max_changes << " changes" << std::endl;

        auto& hypergraph_p =  strategy->init();

        size_t log_step_size = max_changes * context.dynamic.logging_step_size_pct;

        auto duration_sum = std::chrono::high_resolution_clock::duration::zero();

        for (size_t i = 0; i < max_changes; ++i) {
          Change& change = changes[i];
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
          strategy->partition(change, max_changes);
          auto duration = std::chrono::high_resolution_clock::now() - start;
          duration_sum += duration;
          if (log_step_size != 0 && i % log_step_size != 0) {
            continue;
          }
          log_km1_live(i+1, max_changes, context, hypergraph_p, duration_sum);
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

