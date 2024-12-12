#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Connectivity : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        int repartition_count = 0;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          //ASSERT thet num_removed_nodes is correct
          ASSERT([&]() {
              size_t num_removed_nodes = 0;
              for (size_t hn = 0; hn < hypergraph_s.initialNumNodes(); hn++) {
                if (!hypergraph_s.nodeIsEnabled(hn)) {
                  num_removed_nodes++;
                }

              }
              std::cout << "num_removed_nodes: " << num_removed_nodes << " " << hypergraph_s.numRemovedHypernodes() << std::endl;
              return num_removed_nodes == hypergraph_s.numRemovedHypernodes();
          } (), "Number of removed nodes is not correct.");
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
          repartition_count++;

          if (!context.dynamic.use_final_weight) {
            // TODO check if this is necessary
            ASSERT(context.partition.use_individual_part_weights == false);
            context.partition.perfect_balance_part_weights.clear();
            context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
                    hypergraph_s.totalWeight()
                    / static_cast<double>(context.partition.k)));
            context.partition.max_part_weights.clear();
            context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
                                                                                                   * context.partition.perfect_balance_part_weights[0]);
          }
        }
    public:

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
            repartition(hypergraph, context);
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

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
          history.push_back(partition_result);
        }

        void compute_km1_and_imbalance(ds::StaticHypergraph& hypergraph, Context &context, Change change, PartitionResult& partition_result) override {
          partition_result.km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
          partition_result.imbalance = mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context);
        }

        void printFinalStats(ds::StaticHypergraph &hypergraph, Context &context) override {
          (void) hypergraph;
          (void) context;
          std::cout << std::endl;
          std::cout << "Final Stats for Connectivity Strategy" << std::endl;
          std::cout << "Repartition Count: " << repartition_count << std::endl;
        }
    };
}
