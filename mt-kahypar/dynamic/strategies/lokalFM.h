#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>

namespace mt_kahypar::dyn {

    class LokalFM : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph;
        int repartition_count = 0;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _rebalancer;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          partitioned_hypergraph = partition_hypergraph_km1(hypergraph_s, context);
          partitioned_hypergraph_s = std::move(utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph));
          repartition_count++;
        }

        void local_fm(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          GainCachePtr::deleteGainCache(_gain_cache);

          _gain_cache = GainCachePtr::constructGainCache(context);
          _rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph.initialNumNodes(), context, _gain_cache);

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph.initialNumNodes(), hypergraph.initialNumEdges(), context, _gain_cache, *_rebalancer);

          _fm->initialize(partitioned_hypergraph);

          Metrics best_Metrics = {0, 0};

          _fm->refine(partitioned_hypergraph, {hn}, best_Metrics, context.refinement.fm.time_limit_factor);
        }

        PartitionID add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {
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
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_s->partWeight(std::get<1>(block_connectivity)) + hypergraph.nodeWeight(hn) < context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivity));
              return std::get<1>(block_connectivity);
            }
          }
        }

    public:

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change) override {

          //on first call, initialize partitioned_hypergraph_s
          if (!partitioned_hypergraph_s) {
            repartition(hypergraph, context);
          }

          process_change(hypergraph, context, change);

          PartitionResult partition_result = *new PartitionResult();
          partition_result.valid = true;

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

          add_node_to_partitioned_hypergraph(hypergraph, context, hn);

          local_fm(hypergraph, context, hn);

          //check if imbalance is still within bounds else repartition
          if (mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context) > context.partition.epsilon ) {
            repartition(hypergraph, context);
          }

          partition_result.km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
          partition_result.imbalance = mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context);

          history.push_back(partition_result);
        }

        void printFinalStats(ds::StaticHypergraph &hypergraph, Context &context) override {
          (void) hypergraph;
          (void) context;
        }
    };
}
