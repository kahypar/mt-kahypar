#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Hashing : public DynamicStrategy {

    private:

        PartitionID hash_assign(const HypernodeID& hn) {
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            if (PartitionID p_prime = (p + hn) % context.partition.k; partitioned_hypergraph_m.partWeight(p_prime) + hypergraph_m.nodeWeight(hn) <=
                context.partition.max_part_weights[p_prime]) {
              partitioned_hypergraph_m.addNode(hn, p_prime);
              return p_prime;
            }
          }
          // if no partition could accommodate the node put in the first
          partitioned_hypergraph_m.addNode(hn, 0);
          return 0;
        }

    public:

      Hashing(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context) {}
      
        MutablePartitionedHypergraph& init() override {
            partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
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

          for (const HypernodeID& hn : change.added_nodes) {
            const HypernodeID new_hn = hypergraph_m.addHypernode({}, 1);
            (void) new_hn;
            ASSERT(hn == new_hn);
            updateMaxPartWeight(context, hypergraph_m);
            const PartitionID assigned_part = hash_assign(hn);
            (void) assigned_part;
            ASSERT(assigned_part != kInvalidPartition);
          }

          for (const HyperedgeID& he : change.added_edges) {
            hypergraph_m.addHyperedge({}, 1);
            partitioned_hypergraph_m.addEdge(he);
          }

          for (const auto& [node, edge] : change.added_pins)
          {
            hypergraph_m.addPin(edge, node);
            partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
          }
      }
        void printAdditionalFinalStats() override {
            assert(context.dynamic.incremental_km1 == mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) && ("Error: incremental_km1 does not match the quality metric. " + std::to_string(context.dynamic.incremental_km1) + " " + std::to_string(mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1))).c_str());
            std::cout << std::endl << "Final km1: " << context.dynamic.incremental_km1 << " Real km1: " << mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) << std::endl;
          }

    };
}
