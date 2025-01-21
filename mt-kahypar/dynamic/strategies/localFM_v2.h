#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>

namespace mt_kahypar::dyn {

    class LocalFMV2 : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _rebalancer;
        vec<Gain>& _benefit_aggregator;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          context.dynamic.repartition_count++;
          if (!context.dynamic.server) {
            std::cout << "Repartitioning" << std::endl;
          }
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
        }

        void init_local_fm(ds::StaticHypergraph& hypergraph_s, Context& context) {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph_s.initialNumNodes(), context, _gain_cache);

          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = context.dynamic.multitry_localFM;

          _benefit_aggregator = vec<Gain>(context.partition.k, 0);

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_s.initialNumNodes(), hypergraph_s.initialNumEdges(), context, _gain_cache, *_rebalancer);

          GainCachePtr::resetGainCache(_gain_cache);
          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);
          _fm->initialize(partitioned_hypergraph);

          context.dynamic.localFM_round->incremental_km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
        }

        //use local_fm to refine partitioned_hypergraph_s
        void local_fm(Context& context, parallel::scalable_vector<HypernodeID> local_fm_nodes, std::vector<HypernodeID> gain_cache_nodes, parallel::scalable_vector<HypernodeID> rebalance_nodes) {
          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = utils::partitioned_hg_cast(
                  *partitioned_hypergraph_s);

          //update gain cache for all changed nodes
          for (const HypernodeID &hn: gain_cache_nodes) {
            ASSERT(_gain_cache.type == GainPolicy::km1);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCacheEntryForNode(
                    partitioned_hypergraph_s.value(), hn, _benefit_aggregator);
          }

          if (context.dynamic.custom_output_file == "_non-incremental-gain") {
            GainCachePtr::resetGainCache(_gain_cache);
            _fm->initialize(partitioned_hypergraph);
          }

          if (!metrics::isBalanced(*partitioned_hypergraph_s, context)) {
            // use rebalancer to rebalance partitioned_hypergraph_s
            parallel::scalable_vector<parallel::scalable_vector<Move>> moves_by_part;
            Metrics best_Metrics = {mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1),
                                    mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};
            _rebalancer->refineAndOutputMoves(partitioned_hypergraph, {}, moves_by_part, best_Metrics,
                                              std::numeric_limits<double>::max());

            // update gain cache for all nodes in moves_by_part
            for (const parallel::scalable_vector<Move> &moves: moves_by_part) {
              for (const Move &move: moves) {
                // TODO: Is it necessary to update all neighbors?
                for (const HyperedgeID &he: partitioned_hypergraph_s->incidentEdges(move.node)) {
                  for (const HypernodeID &hn2: partitioned_hypergraph_s->pins(he)) {
                    ASSERT(_gain_cache.type == GainPolicy::km1);
                    GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCacheEntryForNode(
                            partitioned_hypergraph_s.value(), hn2, _benefit_aggregator);
                  }
                }
              }
            }

            //deduplicate rebalance_nodes
            std::sort(rebalance_nodes.begin(), rebalance_nodes.end());
            rebalance_nodes.erase(std::unique(rebalance_nodes.begin(), rebalance_nodes.end()), rebalance_nodes.end());

            //append rebalance_nodes to local_fm_nodes
            local_fm_nodes.insert(local_fm_nodes.end(), rebalance_nodes.begin(), rebalance_nodes.end());

            context.dynamic.localFM_round->incremental_km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);

            if (context.dynamic.custom_output_file == "_non-incremental-gain") {
              GainCachePtr::resetGainCache(_gain_cache);
              _fm->initialize(partitioned_hypergraph);
            }
          }

          if (local_fm_nodes.size() == 0) {
            return;
          }

          Metrics best_Metrics = {context.dynamic.localFM_round->incremental_km1,
                                  mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};

          _fm->refine(partitioned_hypergraph, local_fm_nodes, best_Metrics, std::numeric_limits<double>::max());
        }

        PartitionID add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          //compute for each block the number of nodes the new node connected to
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

        LocalFMV2() : _benefit_aggregator(*new vec<Gain>()) {}

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
          repartition(hypergraph, context);
          init_local_fm(hypergraph, context);
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          parallel::scalable_vector<HypernodeID> local_fm_nodes;
          std::vector<HypernodeID> gain_cache_nodes;

          //account for removed edges in incremental km1
          for (const HyperedgeID& he : change.removed_edges) {
            context.dynamic.localFM_round->incremental_km1 -= std::max(partitioned_hypergraph_s->connectivity(he) - 1, 0) * partitioned_hypergraph_s->edgeWeight(he);
          }

          for (const HyperedgeID& he : change.removed_edges) {
            for (const HypernodeID& hn : hypergraph.pins(he)) {
              gain_cache_nodes.push_back(hn);
            }
          }

          for (const HypernodeID& hn : change.removed_nodes) {
            //TODO: mixed queries -> remove node from local_fm_nodes
            //TODO: does this work for multiple node removals?
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
              //skip removed edges
              bool skip = false;
              for (const HyperedgeID &he2 : change.removed_edges) {
                if (he == he2) {
                  skip = true;
                  break;
                }
              }
              if (skip) {
                continue;
              }
              size_t nodes_in_removed_partition_prior_removal = partitioned_hypergraph_s->pinCountInPart(he, partitioned_hypergraph_s->partID(hn));
              //gain increases for remaining node in partition because moving it to another partition will decrease the cut
              if (nodes_in_removed_partition_prior_removal == 2) {
                //append remaining node in partition to local_fm_nodes
                for (const HypernodeID& hn2 : hypergraph.pins(he)) {
                  if (hn2 != hn && partitioned_hypergraph_s->partID(hn2) == partitioned_hypergraph_s->partID(hn)) {
                    local_fm_nodes.push_back(hn2);
                    gain_cache_nodes.push_back(hn2);
//                  break;
                  }
                }
              //only negative gains => update gains but do not refine
              } else if (nodes_in_removed_partition_prior_removal == 1) {
                context.dynamic.localFM_round->incremental_km1 -= partitioned_hypergraph_s->connectivity(he) > 1 ? partitioned_hypergraph_s->connectivity(he) * hypergraph.edgeWeight(he): 0;
                context.dynamic.localFM_round->incremental_km1 += partitioned_hypergraph_s->connectivity(he) > 1 ? (partitioned_hypergraph_s->connectivity(he) - 1) * (hypergraph.edgeWeight(he)) : 0;
                for (const HypernodeID& hn2 : hypergraph.pins(he)) {
                  if (hn2 != hn) {
                    gain_cache_nodes.push_back(hn2);
                  }
                }
              }
            }
            partitioned_hypergraph_s->removeNodePart(hn);
          }


          parallel::scalable_vector<HypernodeID> rebalance_nodes;

          for (const HypernodeID& hn : change.removed_nodes) {
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
              rebalance_nodes.insert(rebalance_nodes.end(), hypergraph.pins(he).begin(), hypergraph.pins(he).end());
            }
          }

          process_change(hypergraph, context, change);

          for (const HypernodeID& hn : change.added_nodes) {
            add_node_to_partitioned_hypergraph(hypergraph, context, hn);
            local_fm_nodes.push_back(hn);
            gain_cache_nodes.push_back(hn);
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
              size_t nodes_in_added_partition = partitioned_hypergraph_s->pinCountInPart(he, partitioned_hypergraph_s->partID(hn));
              //gain increase for all connected nodes because they can be moved here without penalty
              if (nodes_in_added_partition == 1) {
                context.dynamic.localFM_round->incremental_km1 -= std::max(partitioned_hypergraph_s->connectivity(he) - 2, 0) * partitioned_hypergraph_s->edgeWeight(he);
                context.dynamic.localFM_round->incremental_km1 += std::max(partitioned_hypergraph_s->connectivity(he) - 1, 0) * partitioned_hypergraph_s->edgeWeight(he);
                local_fm_nodes.insert(local_fm_nodes.end(), hypergraph.pins(he).begin(), hypergraph.pins(he).end());
                gain_cache_nodes.insert(gain_cache_nodes.end(), hypergraph.pins(he).begin(), hypergraph.pins(he).end());
              //only negative gains => update gains but do not refine
              } else if (nodes_in_added_partition == 2) {
                //TODO: only change gains for node in partition (empiric says this breaks things)
                gain_cache_nodes.insert(gain_cache_nodes.end(), hypergraph.pins(he).begin(), hypergraph.pins(he).end());
              }
            }
          }

          //TODO fix pin counts of added edges
//          for (const HyperedgeID& he : change.added_edges) {
//            std::vector<size_t> pin_counts(context.partition.k, 0);
//            size_t connectivity = 0;
//            for (const HypernodeID& hn : hypergraph.pins(he)) {
//              pin_counts[partitioned_hypergraph_s->partID(hn)]++;
//            }
//            for (PartitionID p = 0; p < context.partition.k; ++p) {
//              if (pin_counts[p] > 0) {
//                connectivity++;
//              }
//              partitioned_hypergraph_s->setPinCountInPart(he, p, pin_counts[p]);
//            }
//          }

          //remove duplicates
          std::sort(gain_cache_nodes.begin(), gain_cache_nodes.end());
          gain_cache_nodes.erase(std::unique(gain_cache_nodes.begin(), gain_cache_nodes.end()), gain_cache_nodes.end());
          std::sort(local_fm_nodes.begin(), local_fm_nodes.end());
          local_fm_nodes.erase(std::unique(local_fm_nodes.begin(), local_fm_nodes.end()), local_fm_nodes.end());

          local_fm(context, local_fm_nodes, gain_cache_nodes, rebalance_nodes);

          ASSERT(metrics::isBalanced(*partitioned_hypergraph_s, context));

          PartitionResult partition_result = *new PartitionResult();
          partition_result.valid = true;
          history.push_back(partition_result);
        }

        void compute_km1_and_imbalance(ds::StaticHypergraph& hypergraph, Context &context, Change change, PartitionResult& partition_result) override {
          partition_result.km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
          if (!context.dynamic.server) {
            std::cout << "km1: " << partition_result.km1 << " Incremental km1: "
                      << context.dynamic.localFM_round->incremental_km1 << std::endl;
          }
          partition_result.imbalance = mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context);
        }

        void printFinalStats(ds::StaticHypergraph &hypergraph, Context &context) override {
          assert(context.dynamic.localFM_round->incremental_km1 == mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1) && ("Error: incremental_km1 does not match the quality metric. " + std::to_string(context.dynamic.localFM_round->incremental_km1) + " " + std::to_string(mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1))).c_str());
          std::cout << std::endl << "Final km1: " << context.dynamic.localFM_round->incremental_km1 << " Real km1: " << mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1) << std::endl;
          (void) hypergraph;
          (void) context;
        }
    };
}
