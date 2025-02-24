#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>

namespace mt_kahypar::dyn {

    class VCycle : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRebalancer> _rebalancer;
        std::unique_ptr<IRefiner> _fm;
        HyperedgeWeight prior_total_weight = 0;
        HyperedgeWeight changed_weight = 0;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          context.dynamic.repartition_count++;
          if (!context.dynamic.server) {
            std::cout << "Repartitioning" << std::endl;
          }
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
        }

        void init_rebalancer(ds::StaticHypergraph& hypergraph_s, Context& context) {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph_s.initialNumNodes(), context, _gain_cache);

          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = context.dynamic.multitry_localFM;

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_s.initialNumNodes(), hypergraph_s.initialNumEdges(), context, _gain_cache, *_rebalancer);

          GainCachePtr::resetGainCache(_gain_cache);
          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);
          _fm->initialize(partitioned_hypergraph);
        }

        //use rebalancer to rebalance partitioned_hypergraph_s
        void rebalance(Context& context) {
          if (!context.dynamic.server) {
            std::cout << "Rebalancing" << std::endl;
          }

          GainCachePtr::resetGainCache(_gain_cache);

          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);
          parallel::scalable_vector<parallel::scalable_vector<Move>> moves_by_part;
          Metrics best_Metrics = {mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1),
                                  mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};

          _fm->initialize(partitioned_hypergraph);

          _rebalancer->refineAndOutputMoves(partitioned_hypergraph, {}, moves_by_part, best_Metrics, std::numeric_limits<double>::max());
        }

        PartitionID add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          //compute for each block the number of nodes the new node connected to
          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
            ASSERT(partitioned_hypergraph_s->checkConnectivitySet(he, context.partition.k));
            for ( const PartitionID& p : partitioned_hypergraph_s->connectivitySet(he) ) {
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), std::greater<std::tuple<int,int>>());

          //Add node to block with highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_s->partWeight(std::get<1>(block_connectivity)) + hypergraph.nodeWeight(hn) <
                context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivity));
              return std::get<1>(block_connectivity);
            }
          }
          // if no partition could accomodate the node put in the best
          partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivities[0]));
        }

    public:

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
          repartition(hypergraph, context);
          init_rebalancer(hypergraph, context);
          prior_total_weight = hypergraph.totalWeight();
          changed_weight = 0;
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          for (const HypernodeID& hn : change.removed_nodes) {
            changed_weight += hypergraph.nodeWeight(hn);
            partitioned_hypergraph_s->removeNodePart(hn);
          }

          process_change(hypergraph, context, change);

          PartitionResult partition_result = *new PartitionResult();

          partition_result.valid = true;

          for (const HypernodeID& hn : change.added_nodes) {
            changed_weight += hypergraph.nodeWeight(hn);
            add_node_to_partitioned_hypergraph(hypergraph, context, hn);
          }

          //reset pin counts of added edges
          for (const HyperedgeID& he : change.added_edges) {
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_s->pinCountInPart(he, p) > 0) {
                partitioned_hypergraph_s->decrementPinCountOfBlockWrapper(he, p);
              }
            }
            for (const HypernodeID& hn : hypergraph.pins(he)) {
              partitioned_hypergraph_s->incrementPinCountOfBlockWrapper(he, partitioned_hypergraph_s->partID(hn));
            }
          }

          if (changed_weight > context.dynamic.step_size_pct * prior_total_weight) {
            mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = utils::partitioned_hg_cast(
                    *partitioned_hypergraph_s);
            HyperedgeWeight prior_km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);

            context.partition.num_vcycles = 1;

            PartitionerFacade::improve(partitioned_hypergraph, context);
            HyperedgeWeight post_km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
            if (!context.dynamic.server) {
              std::cout << std::endl << "Improvement: " << prior_km1 - post_km1 << std::endl;
            }

            context.partition.num_vcycles = 0;
            changed_weight = 0;
            prior_total_weight = hypergraph.totalWeight();
          }

          if (!metrics::isBalanced(*partitioned_hypergraph_s, context)) {
            rebalance(context);
          }

          ASSERT(metrics::isBalanced(*partitioned_hypergraph_s, context));

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
