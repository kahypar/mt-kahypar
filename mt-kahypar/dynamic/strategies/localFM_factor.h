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
          context.dynamic.repartition_count++;
          if (!context.dynamic.server) {
            std::cout << "Repartitioning" << std::endl;
          }
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
          repartition_count++;
        }

        void init_local_fm(ds::StaticHypergraph& hypergraph_s, Context& context) {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph_s.initialNumNodes(), context, _gain_cache);

          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = context.dynamic.multitry_localFM;

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_s.initialNumNodes(), hypergraph_s.initialNumEdges(), context, _gain_cache, *_rebalancer);
        }

        //use local_fm to refine partitioned_hypergraph_s
        void local_fm(ds::StaticHypergraph& hypergraph, Context& context, parallel::scalable_vector<HypernodeID> local_fm_nodes) {

          GainCachePtr::resetGainCache(_gain_cache);

          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);

          if (!metrics::isBalanced(*partitioned_hypergraph_s, context)) {
            // use rebalancer to rebalance partitioned_hypergraph_s
            parallel::scalable_vector<parallel::scalable_vector<Move>> moves_by_part;
            Metrics best_Metrics = {mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1),
                                    mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};
            _rebalancer->refineAndOutputMoves(partitioned_hypergraph, {}, moves_by_part, best_Metrics, std::numeric_limits<double>::max());
          }

          //TODO: is second reset after rebalancing necessary?
          GainCachePtr::resetGainCache(_gain_cache);

          _fm->initialize(partitioned_hypergraph);

          Metrics best_Metrics = {mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1),
                                  mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};

          _fm->refine(partitioned_hypergraph, local_fm_nodes, best_Metrics, std::numeric_limits<double>::max());

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

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
          repartition(hypergraph, context);
          init_local_fm(hypergraph, context);
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          //on first call, initialize partitioned_hypergraph_s
          if (!partitioned_hypergraph_s) {
            nodes_to_partition = parallel::scalable_vector<HypernodeID>(changes_size);
            repartition(hypergraph, context);
          }

          for (const HypernodeID& hn : change.removed_nodes) {
            partitioned_hypergraph_s->removeNodePart(hn);
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
              for (const HypernodeID& hn2 : hypergraph.pins(he)) {
                if (std::find(nodes_to_partition.begin(), nodes_to_partition.end(), hn2) == nodes_to_partition.end()) {
                  nodes_to_partition.push_back(hn2);
                }
              }
            }
          }

          process_change(hypergraph, context, change);

          PartitionResult partition_result = *new PartitionResult();

          partition_result.valid = true;

          for (const HypernodeID &hn : change.added_nodes) {
            add_node_to_partitioned_hypergraph(hypergraph, context, hn);
            if (std::find(nodes_to_partition.begin(), nodes_to_partition.end(), hn) == nodes_to_partition.end()) {
              nodes_to_partition.push_back(hn);
            }
          }

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
