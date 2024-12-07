#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>

namespace mt_kahypar::dyn {

    class LocalFMFactor : public DynamicStrategy {

    private:
        size_t skipped_changes = 0;
        size_t step_size = 1;
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        int repartition_count = 0;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _rebalancer;
        parallel::scalable_vector<HypernodeID> nodes_to_partition;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          std::cout << "Repartitioning" << std::endl;
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
          _gain_cache = GainCachePtr::constructGainCache(context);
          _rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph_s.initialNumNodes(), context, _gain_cache);

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_s.initialNumNodes(), hypergraph_s.initialNumEdges(), context, _gain_cache, *_rebalancer);
          repartition_count++;
        }

        //use local_fm to refine partitioned_hypergraph_s
        void local_fm(ds::StaticHypergraph& hypergraph, Context& context) {

          //GainCachePtr::deleteGainCache(_gain_cache);
          //TODO maybe
          GainCachePtr::resetGainCache(_gain_cache);

//          _gain_cache = GainCachePtr::constructGainCache(context);

          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);

          _fm->initialize(partitioned_hypergraph);

          Metrics best_Metrics = {mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1),
                                  mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};

          _fm->refine(partitioned_hypergraph, nodes_to_partition, best_Metrics, std::numeric_limits<double>::max());
        }

        PartitionID add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          nodes_to_partition.push_back(hn);

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
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_s->partWeight(std::get<1>(block_connectivity)) + hypergraph.nodeWeight(hn) < context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              partitioned_hypergraph_s->setNodePart(hn, std::get<1>(block_connectivity));
              return std::get<1>(block_connectivity);
            }
          }
        }

    public:

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          //on first call, initialize partitioned_hypergraph_s
          if (!partitioned_hypergraph_s) {
            nodes_to_partition = parallel::scalable_vector<HypernodeID>(changes_size);
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

          //skip changes until step_size is reached
          if (skipped_changes >= step_size) {
            skipped_changes = 0;
            step_size *= 2;
            local_fm(hypergraph, context);
            nodes_to_partition = parallel::scalable_vector<HypernodeID>(changes_size);
          } else {
            skipped_changes++;
          }

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
