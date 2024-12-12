#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class NeverRepartition : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;

        size_t suboptimal_decision_count = 0;

        void innitial_partitioning(ds::StaticHypergraph& hypergraph_s, Context& context) {
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

        bool add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          //compute for each block the number of nodes it is connected to
          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_s->connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), std::greater<std::tuple<int,int>>());

          //Add node to block with highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (size_t i = 0; i < context.partition.k; i++) {
            std::tuple<int,int> block_connectivity = block_connectivities[i];
            if (partitioned_hypergraph_s->partWeight(std::get<1>(block_connectivity)) + hypergraph.nodeWeight(hn) < context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivity));
              if (i > 0) {
                suboptimal_decision_count++;
              }
              return true;
            }
          }

          //if no block can accomodate the node, add it to the block with the highest connectivity
          partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivities[0]));
          return false;
        }
    public:

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
          innitial_partitioning(hypergraph, context);
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          process_change(hypergraph, context, change);

          PartitionResult partition_result = *new PartitionResult();

          for (const HypernodeID& hn : change.added_nodes) {
            partition_result.valid = add_node_to_partitioned_hypergraph(hypergraph, context, hn);
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
        }
    };
}
