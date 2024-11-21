#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Connectivity : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        int repartition_count = 0;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph_s, context);
          partitioned_hypergraph_s = std::move(utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph));
          repartition_count++;
        }
    public:

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change) override {

          //on first call, initialize partitioned_hypergraph_s
          if (!partitioned_hypergraph_s) {
            repartition(hypergraph, context);
          }

          process_change(hypergraph, context, change);

          PartitionResult partition_result = *new PartitionResult();
          partition_result.valid = false;

          if (change.added_nodes.empty()) {
            history.push_back(partition_result);
            return;
          }

          if (change.added_nodes.size() > 1) {
            repartition(hypergraph, context);
            history.push_back(partition_result);
            return;
          }

          partition_result.valid = true;

          const HypernodeID& hn = change.added_nodes[0];

          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_s->connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          auto max_connectivity = std::max_element(block_connectivities.begin(), block_connectivities.end());


          partitioned_hypergraph_s->setNodePart(hn, std::get<1>(*max_connectivity));


          //check if imbalance is still within bounds else repartition
          if (partitioned_hypergraph_s->partWeight(std::get<1>(*max_connectivity)) /
              static_cast<double>(context.partition.perfect_balance_part_weights[std::get<1>(*max_connectivity)]) - 1.0 > context.partition.epsilon ) {
            repartition(hypergraph, context);
          }

          partition_result.km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
          partition_result.imbalance = mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context);

          history.push_back(partition_result);
        }

        void printFinalStats(ds::StaticHypergraph &hypergraph, Context &context) override {
          std::cout << "Final Stats for Connectivity Strategy" << std::endl;
          std::cout << "Repartition Count: " << repartition_count << std::endl;
        }
    };
}
