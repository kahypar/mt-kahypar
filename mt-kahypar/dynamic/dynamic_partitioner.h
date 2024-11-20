#pragma once

#include "mt-kahypar/dynamic/strategies/dynamic_strategy.h"
#include "mt-kahypar/dynamic/dynamic_io.h"

namespace mt_kahypar::dyn {

    void partition(Context& context) {

      auto [changes, hypergraph] = generateChanges(context);

      // first_fitting_partition_strategy(hypergraph, context, disabled_nodes);
      // repartition_strategy(hypergraph, context, disabled_nodes);
      // highest_connectivity_partition_strategy(hypergraph, context, disabled_nodes, start_id);
      //repartition_all(hypergraph, context, changes, 0.01);
      repartition_x_connectivity_partition_strategy(hypergraph, context, changes);

      utils::delete_hypergraph(hypergraph);
    }
}

