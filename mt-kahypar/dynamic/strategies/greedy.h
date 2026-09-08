#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>
#include <mt-kahypar/partition/refinement/rebalancing/incremental_rebalancer.h>

namespace mt_kahypar::dyn {

    class Greedy : public DynamicStrategy {

    private:

      PartitionID move_node_to_best_partition(const HypernodeID& hn) {
          const PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);

          // compute for each block the number of nodes the new node connected to
          std::vector block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            ASSERT(partitioned_hypergraph_m.checkConnectivitySet(he, context.partition.k));
            for ( const PartitionID& p : partitioned_hypergraph_m.connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), [](const auto& a, const auto& b) {
            return std::get<0>(a) > std::get<0>(b);
          });

          //Add node to block with the highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (std::get<1>(block_connectivity) != current_part)
            {
              if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <=
                  context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
                partitioned_hypergraph_m.changeNodePart(hn, current_part,  std::get<1>(block_connectivity));
                return std::get<1>(block_connectivity);
              }
            }
          }
          // if no other partition could accommodate the node, stay in current partition
          return current_part;
        }

        PartitionID assign_node_first_free_partition(const HypernodeID& hn) {
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

      Greedy(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context)
        {
        }

        MutablePartitionedHypergraph& init() override {
            partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
            context.dynamic.incremental_km1 = metrics::quality(partitioned_hypergraph_m, Objective::km1);
          return partitioned_hypergraph_m;
        }

        void partition(Change& change, size_t changes_size) override {

          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

          // ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

          for (const auto& [hn, he] : change.removed_pins)
          {
            partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, partitioned_hypergraph_m.partID(hn));
            hypergraph_m.deletePin(he, hn);
          }

          for (const HypernodeID& hn : change.removed_nodes) {
            partitioned_hypergraph_m.removeNodePart(hn);
            hypergraph_m.deleteHypernode(hn);
            updateMaxPartWeight(context, hypergraph_m);
          }

          for (const HyperedgeID& he : change.removed_edges) {

            context.dynamic.incremental_km1 -= std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_m.pinCountInPart(he, p) > 0) {
                partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, p);
              }
            }
            hypergraph_m.deleteHyperedge(he);
          }

          for (const HypernodeID& hn : change.added_nodes) {
            const HypernodeID new_hn = hypergraph_m.addHypernode({}, 1);
            (void) new_hn;
            ASSERT(hn == new_hn);
            updateMaxPartWeight(context, hypergraph_m);
            const PartitionID assigned_part = assign_node_first_free_partition(hn);
            (void) assigned_part;
            ASSERT(assigned_part != kInvalidPartition);
          }

          for (const HyperedgeID& he : change.added_edges) {
            const HyperedgeID new_he = hypergraph_m.addHyperedge({}, 1);
            (void) new_he;
            ASSERT(he == new_he);
            partitioned_hypergraph_m.addEdge(he);
          }

          for (const auto& [node, edge] : change.added_pins)
          {
            hypergraph_m.addPin(edge, node);
            partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
          }

          for (const HypernodeID& hn : change.added_nodes)
          {
            move_node_to_best_partition(hn);
          }

          auto processing_duration_sum = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.processing_duration_sum += processing_duration_sum;

          // ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

        }

        void printAdditionalFinalStats() override {}
    };
}
