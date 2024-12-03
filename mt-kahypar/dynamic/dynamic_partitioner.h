#pragma once

#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/connectivity.h"
#include "mt-kahypar/dynamic/strategies/lokalFM.h"
#include "mt-kahypar/dynamic/strategies/never_repartition.h"
#include "mt-kahypar/dynamic/dynamic_io.h"

namespace mt_kahypar::dyn {

    void partition(Context& context) {

      auto [changes, hypergraph] = generateChanges(context);

      ds::StaticHypergraph& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);

      // If the number of changes is not specified or is greater than the number of changes in the file, we process all the changes
      size_t changes_size = context.dynamic.max_changes == 0 ? changes.size() : std::min((size_t) context.dynamic.max_changes, changes.size());

      DynamicStrategy* strategy;

      if (context.dynamic.strategy == "connectivity") {
        strategy = new Connectivity();
      } else if (context.dynamic.strategy == "repartition") {
        strategy = new Repartition();
      } else if (context.dynamic.strategy == "lokalFM") {
        strategy = new LokalFM();
      } else if (context.dynamic.strategy == "never_repartition") {
        strategy = new NeverRepartition();
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      initOutputFile(context);

      try {

        std::cout << "Processing " << context.dynamic.max_changes << " changes" << std::endl;

        for (size_t i = 0; i < changes_size; ++i) {
          Change& change = changes[i];
          strategy->partition(hypergraph_s, context, change);
          if (!context.dynamic.server && *(&strategy->history.back().valid)) {
            print_progress_bar(i, context.dynamic.max_changes, &strategy->history);
          }
          log_km1_live(i, context, strategy->history.back());
        }

        std::cout << std::endl;

        strategy->printFinalStats(hypergraph_s, context);
        //log_km1(context, &strategy->history);

        utils::delete_hypergraph(hypergraph);

        } catch (std::exception& e) {
          std::cerr << "Error: " << e.what() << std::endl;
          generateErrorFile(context, strategy, e);
          exit(1);
        }
    }
}

