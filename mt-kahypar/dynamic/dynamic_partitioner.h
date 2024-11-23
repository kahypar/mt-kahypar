#pragma once

#include "mt-kahypar/dynamic/strategies/repartition.h"
#include "mt-kahypar/dynamic/strategies/connectivity.h"
#include "mt-kahypar/dynamic/dynamic_io.h"

namespace mt_kahypar::dyn {

    void partition(Context& context) {

      auto [changes, hypergraph] = generateChanges(context);

      ds::StaticHypergraph& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);

      DynamicStrategy* strategy;

      if (context.dynamic.strategy == "connectivity") {
        strategy = new Connectivity();
      } else if (context.dynamic.strategy == "repartition") {
        strategy = new Repartition();
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      std::cout << "Processing " << changes.size() << " changes" << std::endl;

      for (size_t i = 0; i < changes.size(); ++i) {
        Change& change = changes[i];
        strategy->partition(hypergraph_s, context, change);
        if (*(&strategy->history.back().valid)) {
          print_progress_bar(i, changes.size(), &strategy->history);
        }
      }

      std::cout << std::endl;

      strategy->printFinalStats(hypergraph_s, context);
      log_km1(context, &strategy->history);

      utils::delete_hypergraph(hypergraph);
    }
}

