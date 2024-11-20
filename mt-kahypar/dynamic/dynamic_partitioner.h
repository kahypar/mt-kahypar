#pragma once

#include "mt-kahypar/dynamic/strategies/dynamic_strategy.h"
#include "mt-kahypar/dynamic/dynamic_io.h"

namespace mt_kahypar::dyn {

    void partition(Context& context) {

      auto [changes, hypergraph] = generateChanges(context);

      if (context.dynamic.strategy == "connectivity") {
        repartition_x_connectivity_partition_strategy(hypergraph, context, changes);
      } else if (context.dynamic.strategy == "repartition") {
        repartition_strategy(hypergraph, context, changes, 0.01);
      } else {
        throw std::runtime_error("Unknown dynamic strategy: " + context.dynamic.strategy);
      }

      utils::delete_hypergraph(hypergraph);
    }
}

