#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>

namespace mt_kahypar::dyn {

    class Repartition : public DynamicStrategy {
    private:
        ds::PartitionedHypergraph<ds::StaticHypergraph> partitioned_hypergraph_s;

      size_t skipped_changes = 0;
      size_t initial_num_enabled_nodes = 0;
      size_t step_size = 0;
    public:

      void init(ds::StaticHypergraph& hypergraph, Context& context) override {
        (void) hypergraph;
        (void) context;
      }

      void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

        process_change(hypergraph, context, change);
        PartitionResult partition_result = *new PartitionResult();

        if (initial_num_enabled_nodes == 0) {
          for (HypernodeID hn = 0; hn < hypergraph.initialNumNodes(); hn++) {
            if (hypergraph.nodeIsEnabled(hn)) {
              initial_num_enabled_nodes++;
            }
          }
        }

        if (step_size == 0) {
          size_t max_changes = context.dynamic.max_changes == 0 ? changes_size : std::min((size_t) context.dynamic.max_changes, changes_size);
          step_size = static_cast<size_t>(context.dynamic.step_size_pct * max_changes);
        }

        if (skipped_changes < step_size) {
          partition_result.valid = false;
          history.push_back(partition_result);
          skipped_changes++;
          return;
        }

        partition_result.valid = true;
        skipped_changes = 0;

        partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph, context);


        if (!context.dynamic.use_final_weight) {
          // TODO check if this is necessary
          ASSERT(context.partition.use_individual_part_weights == false);
          context.partition.perfect_balance_part_weights.clear();
          context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
                  hypergraph.totalWeight()
                  / static_cast<double>(context.partition.k)));
          context.partition.max_part_weights.clear();
          context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
                                                                                                 * context.partition.perfect_balance_part_weights[0]);
        }

        history.push_back(partition_result);
      }

      void compute_km1_and_imbalance(ds::StaticHypergraph& hypergraph, Context &context, Change change, PartitionResult& partition_result) override {
        partition_result.km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1);
        partition_result.imbalance = mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
      }

      void printFinalStats(ds::StaticHypergraph& hypergraph, Context& context) override {
        (void) hypergraph;
        (void) context;
      };

    };
}