#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>
#include <mt-kahypar/partition/refinement/rebalancing/incremental_rebalancer.h>

namespace mt_kahypar::dyn {

    class Greedy : public DynamicStrategy {

    private:

        PartitionID add_node_to_partitioned_hypergraph(const HypernodeID& hn) {

          // partitioned_hypergraph_m.addNode(hn, kInvalidPartition);

          //compute for each block the number of nodes the new node connected to
          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_m.connectivitySet(he) ) {
              ASSERT(partitioned_hypergraph_m.checkConnectivitySet(he, context.partition.k));
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), std::greater<std::tuple<int,int>>());

          //Add node to block with the highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <
                context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              // partitioned_hypergraph_m.setNodePart(hn, std::get<1>(block_connectivity));
              partitioned_hypergraph_m.addNode(hn, std::get<1>(block_connectivity));
              return std::get<1>(block_connectivity);
            }
          }
          // if no partition could accomodate the node put in the best
          partitioned_hypergraph_m.addNode(hn, std::get<1>(block_connectivities[0]));
          return std::get<1>(block_connectivities[0]);
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

          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1));
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

          for (const auto& [hn, he] : change.removed_pins)
          {
            size_t pin_count_in_part_prior_removal = partitioned_hypergraph_m.pinCountInPart(he, partitioned_hypergraph_m.partID(hn));

            // changed_weight += hypergraph_m.edgeWeight(he)/hypergraph_m.edgeSize(he);

            //decrement km1 if pin is single pin in partition for this edge prior to removal
            if (pin_count_in_part_prior_removal == 1 &&
                partitioned_hypergraph_m.connectivity(he) > 1)
            {
              context.dynamic.incremental_km1 -= partitioned_hypergraph_m.edgeWeight(he);
            }
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
                // ASSERT(false);
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
            const PartitionID assigned_part = add_node_to_partitioned_hypergraph(hn);
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

            //increment km1 if pin is single pin in partition for this edge after addition
            if (partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(node)) == 1 &&
                partitioned_hypergraph_m.connectivity(edge) > 1)
            {
              context.dynamic.incremental_km1 += partitioned_hypergraph_m.edgeWeight(edge);
            }
          }

          auto processing_duration_sum = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.processing_duration_sum += processing_duration_sum;

          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

        }

        void printAdditionalFinalStats() override {
          assert(context.dynamic.incremental_km1 == mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) && ("Error: incremental_km1 does not match the quality metric. " + std::to_string(context.dynamic.incremental_km1) + " " + std::to_string(mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1))).c_str());
          std::cout << std::endl << "Final km1: " << context.dynamic.incremental_km1 << " Real km1: " << mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) << std::endl;
        }
    };
}
