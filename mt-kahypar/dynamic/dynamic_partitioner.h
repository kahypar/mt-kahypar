#pragma once

#include <mt-kahypar/dynamic/strategies/localFM_factor.h>
#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/connectivity.h"
#include "mt-kahypar/dynamic/strategies/localFM.h"
#include "mt-kahypar/dynamic/strategies/never_repartition.h"
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
        resetHypergraph(hypergraph_s, changes, context);
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
      } else if (context.dynamic.strategy == "never_repartition") {
        strategy = new NeverRepartition();
      } else if (context.dynamic.strategy == "localFM_factor") {
        strategy = new LocalFMFactor();
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      mt_kahypar::LocalFMRound localFM_round = mt_kahypar::LocalFMRound();
      localFM_round.overall_improvement = 0;
      localFM_round.touched_nodes = 0;
      context.dynamic.localFM_round = &localFM_round;

      try {

        std::cout << "Processing " << max_changes << " changes" << std::endl;

        strategy->init(hypergraph_s, context);

        for (size_t i = 0; i < max_changes; ++i) {
          Change& change = changes[i];
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
          strategy->partition(hypergraph_s, context, change, max_changes);
          log_km1_live(i, context, strategy->history.back(), std::chrono::high_resolution_clock::now() - start);
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

