#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

#include "mt-kahypar/io/command_line_options.h"

namespace mt_kahypar::dyn {

    class Greedy : public DynamicStrategy {

    private:
      std::vector<PartitionID> partition_assignment;
      std::vector<size_t> partition_connectivities;
      std::vector<HypernodeWeight> partition_weights;
      HypernodeWeight global_partition_weight;
      bool first_change = true;
    public:

      Greedy(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context) {}
      
        MutablePartitionedHypergraph& init() override {
            partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
          partition_connectivities = std::vector<size_t>(context.partition.k, 0);
          partition_weights = std::vector<HypernodeWeight>(context.partition.k, 0);
          global_partition_weight = 0;
          for (HypernodeID hn = 0; hn < 100; ++hn)
          {
            partition_assignment.push_back(partitioned_hypergraph_m.partID(hn));
            partition_weights[partitioned_hypergraph_m.partID(hn)]++;
            global_partition_weight++;
          }
          return partitioned_hypergraph_m;
        }

        void partition(Change& change, size_t changes_size) override {
          (void) changes_size;
          if (first_change) {
            first_change = false;
            partition_assignment.reserve(changes_size);
          }
          // ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
          // ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());

          partition_connectivities = std::vector<size_t>(context.partition.k, 0);
          // Removals are not supported in streaming mode

          ASSERT(change.removed_nodes.empty());
          ASSERT(change.removed_edges.empty());
          ASSERT(change.removed_pins.empty());

          HypernodeID added_node = change.added_nodes.front();


          for (const auto& [node, edge] : change.added_pins)
          {
            if (node != added_node) {
              partition_connectivities[partition_assignment[node]]++;
            }
          }

          PartitionID added_partition = 0;
          for (PartitionID p = 0; p < context.partition.k; ++p)
          {
            if ((partition_connectivities[p] > partition_connectivities[added_partition] && 1 + partition_weights[p] <= ((double)global_partition_weight / context.partition.k) * (1.0 + context.partition.epsilon)) || 1 + partition_weights[added_partition] > ((double)global_partition_weight / context.partition.k) * (1.0 + context.partition.epsilon)) {
              added_partition = p;
            }
          }
          partition_assignment.push_back(added_partition);
          partition_weights[added_partition]++;
          global_partition_weight++;

          // compute and print the current imbalance if the nodeid is a multiple of 1000
          // if (added_node % 1000 == 0)
          // {
          //   HypernodeWeight max_partition_weight = 0;
          //   for (PartitionID p = 0; p < context.partition.k; ++p) {
          //     max_partition_weight = std::max(max_partition_weight, partition_weights[p]);
          //   }
          //   double imbalance = (double)max_partition_weight / (double)global_partition_weight;
          //   std::cout << "Current imbalance after adding node " << added_node << ": " << imbalance << std::endl;
          //   std::cout << "Global partition weight: " << global_partition_weight << std::endl;
          //   std::cout << "Partition weights: ";
          //   for (PartitionID p = 0; p < context.partition.k; ++p)
          //   {
          //     std::cout << partition_weights[p] << " ";
          //   }
          //
          // }
      }
      void printAdditionalFinalStats() override {
        // write partitioning to context.partition.graph_partition_filename

        auto filename = context.partition.graph_partition_filename;
        if (filename.empty()) {
          throw InvalidInputException("No filename for output partition file specified");
        }
        std::ofstream out_stream(filename.c_str());
        for (const PartitionID p : partition_assignment)
        {
          out_stream << p << std::endl;
        }
        out_stream.close();
      }

    };
}
