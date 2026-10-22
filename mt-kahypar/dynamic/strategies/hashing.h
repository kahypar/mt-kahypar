#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Hashing : public DynamicStrategy {

    private:
      std::vector<PartitionID> partition_assignment;
      HypernodeWeight global_partition_weight;
      bool first_change = true;

    public:

      Hashing(ds::MutableHypergraph& hypergraph_m, Context& context)
        : DynamicStrategy(hypergraph_m, context), global_partition_weight(0)
      {
      }

      MutablePartitionedHypergraph& init() override {
          partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          context.dynamic.incremental_km1 = metrics::quality(partitioned_hypergraph_m, Objective::km1);
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
           global_partition_weight = 0;
          for (HypernodeID hn = 0; hn < 100; ++hn)
          {
            partition_assignment.push_back(partitioned_hypergraph_m.partID(hn));
            global_partition_weight++;
          }
          return partitioned_hypergraph_m;
        }

        void partition(Change& change, size_t changes_size) override {
          (void) changes_size;
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
          // ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());

          // Removals are not supported in streaming mode

          ASSERT(change.removed_nodes.empty());
          ASSERT(change.removed_edges.empty());
          ASSERT(change.removed_pins.empty());

          // if (first_change) {
          //   first_change = false;
          //   partition_assignment.reserve(changes_size);
          // }

          for (const HypernodeID& hn : change.added_nodes) {
            PartitionID part = hn % context.partition.k;
            partition_assignment.push_back(part);
          }

      }

      void printAdditionalFinalStats() override {
        // write partitioning to context.partition.graph_partition_filename
        // auto filename = context.partition.graph_partition_filename;
        // if (filename.empty()) {
        //   throw InvalidInputException("No filename for output partition file specified");
        // }
        // std::cout << "Writing final partitioning to file: " << filename << std::endl;
        // std::ofstream out_stream(filename.c_str());
        // for (const PartitionID p : partition_assignment)
        // {
        //   out_stream << p << std::endl;
        // }
        // out_stream.close();
      }
    };
}
