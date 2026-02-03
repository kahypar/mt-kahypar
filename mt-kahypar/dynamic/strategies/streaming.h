#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>
#include <mt-kahypar/partition/refinement/rebalancing/incremental_rebalancer.h>
#include <fstream>

namespace mt_kahypar::dyn {

    class Streaming : public DynamicStrategy {

    private:

        PartitionID move_node_to_best_partition(const HypernodeID& hn) {

          PartitionID current_part = partitioned_hypergraph_m.partID(hn);
          ASSERT(current_part != kInvalidPartition);

          //compute for each block the number of nodes the new node connected to
          std::vector block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_m.connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), std::greater<std::tuple<int,int>>());

          //Add node to block with the highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (std::get<1>(block_connectivity) != current_part)
            {
              if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <
                  context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
                partitioned_hypergraph_m.changeNodePart(hn, current_part,  std::get<1>(block_connectivity));
                return std::get<1>(block_connectivity);
              }
            } else
            {
              return current_part;
            }
          }
        }

        PartitionID assign_node_first_free_partition(const HypernodeID& hn) {
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            std::cout << "Checking partition " << p << " with current weight "
                      << partitioned_hypergraph_m.partWeight(p)
                      << " and max weight " << context.partition.max_part_weights[p] << std::endl;
            if (partitioned_hypergraph_m.partWeight(p) + hypergraph_m.nodeWeight(hn) <
                context.partition.max_part_weights[p]) {
              partitioned_hypergraph_m.addNode(hn, p);
              return p;
            }
          }
          // if no partition could accomodate the node put in the first
          partitioned_hypergraph_m.addNode(hn, 0);
          return 0;
        }

    public:

      Streaming(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context) {}
      
        MutablePartitionedHypergraph& init() override {
            partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
          return partitioned_hypergraph_m;
        }

        void partition(Change& change, size_t changes_size) override {
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());

          for (const auto& [hn, he] : change.removed_pins)
          {
            partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, partitioned_hypergraph_m.partID(hn));
            hypergraph_m.deletePin(he, hn);
          }
        std::cout << "Max part weights: ";
        for (PartitionID p = 0; p < context.partition.k; ++p) {
          std::cout << context.partition.max_part_weights[p] << " ";
        }
        std::cout << std::endl;

        std::cout << "Current part weights: ";
        for (PartitionID p = 0; p < context.partition.k; ++p) {
          std::cout << partitioned_hypergraph_m.partWeight(p) << " ";
        }
        std::cout << std::endl;

          for (const HypernodeID& hn : change.removed_nodes) {
            partitioned_hypergraph_m.removeNodePart(hn);
            hypergraph_m.deleteHypernode(hn);
            updateMaxPartWeight(context, hypergraph_m);
          }
        std::cout << "Max part weights: ";
        for (PartitionID p = 0; p < context.partition.k; ++p) {
          std::cout << context.partition.max_part_weights[p] << " ";
        }
        std::cout << std::endl;

        std::cout << "Current part weights: ";
        for (PartitionID p = 0; p < context.partition.k; ++p) {
          std::cout << partitioned_hypergraph_m.partWeight(p) << " ";
        }
        std::cout << std::endl;

          for (const HyperedgeID& he : change.removed_edges) {
            hypergraph_m.deleteHyperedge(he);
          }
        ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

        std::cout << "After removals:" << std::endl;
          for (const HypernodeID& hn : change.added_nodes) {
            std::cout << "Max part weights: ";
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              std::cout << context.partition.max_part_weights[p] << " ";
            }
            std::cout << std::endl;

            std::cout << "Current part weights: ";
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              std::cout << partitioned_hypergraph_m.partWeight(p) << " ";
            }
            std::cout << std::endl;
            const HypernodeID new_hn = hypergraph_m.addHypernode({}, 1);
            updateMaxPartWeight(context, hypergraph_m);
            const PartitionID assigned_part = assign_node_first_free_partition(hn);
            std::cout << "Max part weights: ";
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              std::cout << context.partition.max_part_weights[p] << " ";
            }
            std::cout << std::endl;

            std::cout << "Current part weights: ";
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              std::cout << partitioned_hypergraph_m.partWeight(p) << " ";
            }
            std::cout << std::endl;
          }

        // print max part weights and current part weights


        ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

          for (const HyperedgeID& he : change.added_edges) {
            hypergraph_m.addHyperedge({}, 1);
            partitioned_hypergraph_m.addEdge(he);
          }

          for (const auto& [node, edge] : change.added_pins)
          {
            hypergraph_m.addPin(edge, node);
            partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
          }


          // Reassign newly added nodes to partitions
          // for (const HypernodeID& hn : change.added_nodes) {
          //   const PartitionID assigned_part = add_node_to_partitioned_hypergraph(hn);
          //   ASSERT(assigned_part != kInvalidPartition);
          // }

          // add new pins to block counters
          // for (const auto& [node, edge] : change.added_pins)
          // {
          //   partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
          // }

          // print all part weights
          // for (PartitionID p = 0; p < context.partition.k; ++p) {
          //   std::cout << "Part " << p << " weight: " << partitioned_hypergraph_m.partWeight(p) << std::endl;
          // }

          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
          ASSERT(partitioned_hypergraph_m.checkAllConnectivitySets());
        }
    };
}
