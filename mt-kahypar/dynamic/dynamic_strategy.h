#pragma once

#include <mt-kahypar/partition/metrics.h>
#include <mt-kahypar/dynamic/dynamic_datastructures.h>
#include <mt-kahypar/partition/registries/register_memory_pool.h>
#include <mt-kahypar/partition/partitioner_facade.h>
#include <mt-kahypar/utils/cast.h>
#include <mt-kahypar/utils/delete.h>

namespace mt_kahypar::dyn {

    class DynamicStrategy {
    private:
        void activate_nodes(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          for (const HypernodeID &hn: change.added_nodes) {
            hypergraph.enableHypernodeWithEdges(hn);
            if (!context.dynamic.use_final_weight) {
              hypergraph.incrementTotalWeight(hn);

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
          }
        }
        void deactivate_nodes(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          for (const HypernodeID &hn: change.removed_nodes) {
            hypergraph.disableHypernodeWithEdges(hn);
            if (!context.dynamic.use_final_weight) {
              hypergraph.decrementTotalWeight(hn);
            }
          }
        }
        static void activate_edges(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          (void) context;
          for (const HyperedgeID &he: change.added_edges) {
            hypergraph.restoreEdge(he);
          }
        }
        static void deactivate_edges(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          (void) context;
          for (const HyperedgeID &he: change.removed_edges) {
            hypergraph.removeEdge(he);
          }
        }
        static void activate_pins(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          (void) context;
          for (PinChange &pin_change: change.added_pins) {
            hypergraph.restorePin(pin_change.node, pin_change.edge);
          }
        }
        static void deactivate_pins(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          (void) context;
          for (PinChange &pin_change: change.removed_pins) {
            hypergraph.removePin(pin_change.node, pin_change.edge);
          }
        }
    protected:

        /*
         * Partitions the hypergraph using the mt-kahypar partitioner.
         */
        static ds::PartitionedHypergraph<typename ds::StaticHypergraph>
        partition_hypergraph_km1(ds::StaticHypergraph &hypergraph, Context &context) {

          // copy the hypergraph to make sure we don't modify the original
          ds::StaticHypergraph hypergraph_copy = hypergraph.copy();

//          // Initialize Memory Pool
//          //TODO: do we have to repeat this for every partitioning?
//          register_memory_pool(hypergraph_copy, context);

          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
                  PartitionerFacade::partition(utils::hypergraph_cast(hypergraph_copy), context);

//          // TODO: Do we need to free memory chunks every time?
//          parallel::MemoryPool::instance().free_memory_chunks();
//          TBBInitializer::instance().terminate();


          auto partitioned_hypergraph_s = std::move(utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph));

          partitioned_hypergraph_s.setHypergraph(hypergraph);

          return partitioned_hypergraph_s;
        }


    public:
        struct PartitionResult {
            bool valid;
            HyperedgeWeight km1;
            double imbalance;
        };
        std::vector<PartitionResult> history{};

        /*
         * Process the whole Change object.
         */
        void process_change(ds::StaticHypergraph &hypergraph, Context &context, Change change) {
          activate_nodes(hypergraph, context, change);
          deactivate_nodes(hypergraph, context, change);
          activate_edges(hypergraph, context, change);
          deactivate_edges(hypergraph, context, change);
          activate_pins(hypergraph, context, change);
          deactivate_pins(hypergraph, context, change);

          context.setupPartWeights(hypergraph.totalWeight());
        }

        /*
         * partition() is called for each change in the dynamic hypergraph.
         * The strategy should partition the hypergraph according to the change.
         * The strategy should update the history vector with the result of the partitioning.
         */
        virtual void partition(ds::StaticHypergraph &hypergraph, Context &context, Change change, size_t changes_size) = 0;

        /*
         * A Strategy can print final statistics after the last partitioning step.
         */
        virtual void printFinalStats(ds::StaticHypergraph &hypergraph, Context &context) = 0;
    };

}