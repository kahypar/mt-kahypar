#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Repartition : public DynamicStrategy {
    private:
      size_t skipped_changes = 0;
    public:

      void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change) override {

        process_change(hypergraph, context, change);
        PartitionResult partition_result = *new PartitionResult();

        auto step_size = static_cast<size_t>(context.dynamic.step_size_pct * hypergraph.initialNumNodes());

        if (skipped_changes < step_size) {
          partition_result.valid = false;
          history.push_back(partition_result);
          skipped_changes++;
          return;
        }

        partition_result.valid = true;
        skipped_changes = 0;

        mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
        auto partitioned_hypergraph_s =
                std::move(utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph));

        //TODO: Is this necessary every time?
        utils::delete_partitioned_hypergraph(partitioned_hypergraph);

        partition_result.km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1);
        partition_result.imbalance = mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
        history.push_back(partition_result);
      }

      void printFinalStats(ds::StaticHypergraph& hypergraph, Context& context) override {
      };

    };
}