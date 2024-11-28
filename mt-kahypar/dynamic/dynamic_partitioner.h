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

      if (context.dynamic.max_changes == -1) context.dynamic.max_changes = changes.size();

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

      std::cout << "Processing " << context.dynamic.max_changes << " changes" << std::endl;

      for (size_t i = 0; i < context.dynamic.max_changes; ++i) {
        Change& change = changes[i];
        strategy->partition(hypergraph_s, context, change);
        if (*(&strategy->history.back().valid)) {
          print_progress_bar(i, context.dynamic.max_changes, &strategy->history);
        }
        log_km1_live(i, context, strategy->history.back());
      }

      std::cout << std::endl;

      strategy->printFinalStats(hypergraph_s, context);
      //log_km1(context, &strategy->history);

      utils::delete_hypergraph(hypergraph);
    }
}

