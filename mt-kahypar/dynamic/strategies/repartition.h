#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Repartition : public DynamicStrategy {
    private:
      size_t change_count = -1;
    public:

      Repartition(ds::MutableHypergraph& hypergraph_m, Context& context)
        : DynamicStrategy(hypergraph_m, context) {}


      MutablePartitionedHypergraph& init() override {
          partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          return partitioned_hypergraph_m;
      }

      void partition(Change& change, size_t changes_size) override {
        change_count++;

        process_change(hypergraph_m, context, change);

        size_t log_step_size = changes_size * context.dynamic.logging_step_size_pct;
        if (log_step_size != 0 && change_count % log_step_size != 0) {
          return;
        }

        context.partition.max_part_weights.clear();
        partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);

      }

      void printAdditionalFinalStats() override {
      }

    };
}