#pragma once

#include <mt-kahypar/partition/metrics.h>
#include <mt-kahypar/dynamic_save/dynamic_datastructures.h>
#include <mt-kahypar/partition/registries/register_memory_pool.h>
#include <mt-kahypar/partition/partitioner_facade.h>
#include <mt-kahypar/utils/cast.h>
#include <mt-kahypar/utils/delete.h>

#include "mt-kahypar/io/partitioning_output.h"

namespace mt_kahypar::dyn {

    class DynamicStrategy {
    private:
        static void activate_nodes(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          for (const HypernodeID &hn: change.added_nodes) {
          HypernodeID new_hn = hypergraph.addHypernode({}, 1);
          (void) new_hn;
          (void) hn;
          (void) context;
          ASSERT(new_hn == hn);
          // ASSERT(context.partition.use_individual_part_weights == false);
          // ASSERT(context.partition.perfect_balance_part_weights == std::vector<HypernodeWeight>(context.partition.k, ceil(
          //         hypergraph.totalWeight()
          //         / static_cast<double>(context.partition.k))));
          // ASSERT(context.partition.max_part_weights == std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
          //                                                                               * context.partition.perfect_balance_part_weights[0]));
          // updateMaxPartWeight(context, hypergraph);
        // if (!context.dynamic.use_final_weight) {
        //
        //   // TODO check if this is necessary
        //   ASSERT(context.partition.use_individual_part_weights == false);
        //   context.partition.perfect_balance_part_weights.clear();
        //   context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
        //           hypergraph.totalWeight()
        //           / static_cast<double>(context.partition.k)));
        //   context.partition.max_part_weights.clear();
        //   context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
        //                                                                                          * context.partition.perfect_balance_part_weights[0]);
        // }
          }
        }
        static void deactivate_nodes(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          (void) context;
          for (const HypernodeID &hn: change.removed_nodes) {
            hypergraph.deleteHypernode(hn);
            // updateMaxPartWeight(context, hypergraph);
            (void) hn;
            ASSERT(context.partition.use_individual_part_weights == false);
            //TODO check if assertion is relevant
            // ASSERT(context.partition.perfect_balance_part_weights == std::vector<HypernodeWeight>(context.partition.k, ceil(
            //         hypergraph.totalWeight()
            //         / static_cast<double>(context.partition.k))));
            // ASSERT(context.partition.max_part_weights == std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
            //                                                                               * context.partition.perfect_balance_part_weights[0]));
            // if (!context.dynamic.use_final_weight) {
              // hypergraph.decrementTotalWeight(hn);
              // // TODO check if this is necessary
              // ASSERT(context.partition.use_individual_part_weights == false);
              // context.partition.perfect_balance_part_weights.clear();
              // context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
              //         hypergraph.totalWeight()
              //         / static_cast<double>(context.partition.k)));
              // context.partition.max_part_weights.clear();
              // context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
              //                                                                                        * context.partition.perfect_balance_part_weights[0]);

            // }
          }
        }
        static void activate_edges(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          (void) context;
          for (const HyperedgeID &he: change.added_edges) {
              HyperedgeID new_he = hypergraph.addHyperedge({}, 1);
              (void) he;
              (void) new_he;
              ASSERT(he == new_he);
          }
        }
        static void deactivate_edges(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          (void) context;
          for (const HyperedgeID &he: change.removed_edges) {
            hypergraph.deleteHyperedge(he);
          }
        }
        static void activate_pins(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          (void) context;
          for (const PinChange &pin_change: change.added_pins) {
            hypergraph.addPin(pin_change.edge, pin_change.node);
          }
        }
        static void deactivate_pins(ds::MutableHypergraph &hypergraph, Context &context, const Change &change) {
          (void) context;
          for (const PinChange &pin_change: change.removed_pins) {
            hypergraph.deletePin(pin_change.edge, pin_change.node);
          }
        }
    protected:
        ds::MutableHypergraph& hypergraph_m;
        ds::PartitionedHypergraph<ds::MutableHypergraph> partitioned_hypergraph_m;
        Context& context;

        /*
         * Partitions the hypergraph using the mt-kahypar partitioner.
         */
        static ds::PartitionedHypergraph<ds::MutableHypergraph>
        partition_hypergraph_km1(ds::MutableHypergraph &hypergraph, Context &context) {

          // copy the hypergraph to make sure we don't modify the original
          // ds::MutableHypergraph hypergraph_copy = hypergraph.copy();

//          // Initialize Memory Pool
//          //TODO: do we have to repeat this for every partitioning?
//          register_memory_pool(hypergraph_copy, context);

          // io::printContext(context);
          // io::printHypergraphInfo(hypergraph, context, "Input Hypergraph", true);


          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
                  PartitionerFacade::partition(utils::hypergraph_cast(hypergraph), context);

//          // TODO: Do we need to free memory chunks every time?
//          parallel::MemoryPool::instance().free_memory_chunks();
//          TBBInitializer::instance().terminate();


          auto partitioned_hypergraph_s = std::move(utils::cast<ds::PartitionedHypergraph<typename ds::MutableHypergraph>>(partitioned_hypergraph));

          partitioned_hypergraph_s.setHypergraph(hypergraph);

          return partitioned_hypergraph_s;
        }

        // partition_hypergraph_km1_transform(ds::MutableHypergraph &hypergraph, Context &context) {
        //
        //   std::vector<HypernodeID> static_to_mut_hn;
        //   std::vector<HyperedgeID> static_to_mut_he;
        //   std::vector<HypernodeID> deleted_hn;
        //   std::vector<HyperedgeID> deleted_he;
        //
        //   ds::StaticHypergraph static_hg = hypergraph.toStaticHypergraph(&static_to_mut_hn, &static_to_mut_he, &deleted_hn, &deleted_he);
        //
        //   mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
        //           PartitionerFacade::partition(utils::hypergraph_cast(static_hg), context);
        //
        //   auto partitioned_hypergraph_s = utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);
        //
        //   ds::MutableHypergraph contracted_hg = hypergraph.fromStaticHypergraphSimple(partitioned_hypergraph_s);
        //
        // }


    public:
        virtual ~DynamicStrategy() = default;

        DynamicStrategy(ds::MutableHypergraph &hypergraph, Context &context) :
                hypergraph_m(hypergraph),
                context(context) {
        }

        /*
         * Process the whole setup-part of the Changelist.
         */
        static void process_setup_changes(ds::MutableHypergraph &hypergraph, Context &context, std::vector<Change> changes) {
          ASSERT(context.dynamic.setup_moves_count <= changes.size());
          for (size_t i = 0; i < context.dynamic.setup_moves_count; ++i)
          {
            const Change &change = changes[i];
            deactivate_pins(hypergraph, context, change);
            deactivate_nodes(hypergraph, context, change);
            deactivate_edges(hypergraph, context, change);

            activate_edges(hypergraph, context, change);
            activate_nodes(hypergraph, context, change);
            activate_pins(hypergraph, context, change);
          }
          // context.setupPartWeights(hypergraph.totalWeight());
        }

        static void process_change(ds::MutableHypergraph &hypergraph, Context &context, Change change) {
            deactivate_pins(hypergraph, context, change);
            deactivate_nodes(hypergraph, context, change);
            deactivate_edges(hypergraph, context, change);

            activate_edges(hypergraph, context, change);
            activate_nodes(hypergraph, context, change);
            activate_pins(hypergraph, context, change);
          }

        static void updateMaxPartWeight(Context &context, ds::MutableHypergraph &hypergraph)
        {
          //   // TODO check if this is necessary
          ASSERT(context.partition.use_individual_part_weights == false);
          context.partition.perfect_balance_part_weights.clear();
          context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
                  hypergraph.totalWeight()
                  / static_cast<double>(context.partition.k)));
          context.partition.max_part_weights.clear();
          context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
                                                                                                 * context.partition.perfect_balance_part_weights[0]);
        }

        virtual MutablePartitionedHypergraph& init() = 0;

        /*
         * partition() is called for each change in the dynamic hypergraph.
         * The strategy should partition the hypergraph according to the change.
         * The strategy should update the history vector with the result of the partitioning.
         */
        virtual void partition(Change &change, size_t change_id) = 0;

        /*
         * A Strategy can print final statistics after the last partitioning step.
         */
        virtual void printAdditionalFinalStats() = 0;

        static ds::PartitionedHypergraph<ds::MutableHypergraph>& getPartitionedHypergraphCopy(DynamicStrategy& strategy) {
          return strategy.partitioned_hypergraph_m;
        }
    };

}
