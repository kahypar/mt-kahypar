#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>
#include <mt-kahypar/partition/refinement/rebalancing/incremental_rebalancer.h>

namespace mt_kahypar::dyn {

    class LocalFMRebalanceVCycle : public DynamicStrategy {

    private:
        std::optional<ds::PartitionedHypergraph<ds::StaticHypergraph>> partitioned_hypergraph_s;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _global_rebalancer;
        vec<Gain>& _benefit_aggregator;
        IncrementalRebalancer _rebalancer;
        HyperedgeWeight prior_total_weight = 0;
        HyperedgeWeight changed_weight = 0;

        void repartition(ds::StaticHypergraph& hypergraph_s, Context& context) {
          context.dynamic.repartition_count++;
          if (!context.dynamic.server) {
            std::cout << "Repartitioning" << std::endl;
          }
          partitioned_hypergraph_s = partition_hypergraph_km1(hypergraph_s, context);
        }

        void init_local_fm(ds::StaticHypergraph& hypergraph_s, Context& context) {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _global_rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancer, hypergraph_s.initialNumNodes(), context, _gain_cache);

          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = context.dynamic.multitry_localFM;

          _benefit_aggregator = vec<Gain>(context.partition.k, 0);

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_s.initialNumNodes(), hypergraph_s.initialNumEdges(), context, _gain_cache, *_global_rebalancer);

          GainCachePtr::resetGainCache(_gain_cache);
          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(*partitioned_hypergraph_s);
          _fm->initialize(partitioned_hypergraph);

          context.dynamic.localFM_round->incremental_km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);

          ASSERT(partitioned_hypergraph_s->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
        }

        //use local_fm to refine partitioned_hypergraph_s
        void local_fm(Context& context, parallel::scalable_vector<HypernodeID> local_fm_nodes, std::vector<HypernodeID> gain_cache_nodes, Change change, vec<PartitionID> empty_blocks) {

          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = utils::partitioned_hg_cast(
                  *partitioned_hypergraph_s);

          //update gain cache for all changed nodes
          for (const HypernodeID &hn: gain_cache_nodes) {
            ASSERT(_gain_cache.type == GainPolicy::km1);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCacheEntryForNode(
                    partitioned_hypergraph_s.value(), hn, _benefit_aggregator);
            _rebalancer.insertOrUpdateNode(hn, partitioned_hypergraph_s->partID(hn));
          }

          //pull into emptier blocks
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            auto [gain, moved_nodes] = _rebalancer.pullAndUpdateGainCache(p);
            context.dynamic.localFM_round->incremental_km1 -= gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }

          ASSERT(partitioned_hypergraph_s->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));

          if (!metrics::isBalanced(*partitioned_hypergraph_s, context)) {
            // use rebalancer to rebalance partitioned_hypergraph_s
            auto [gain, moved_nodes] = _rebalancer.rebalanceAndUpdateGainCache();
            context.dynamic.localFM_round->incremental_km1 -= gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }

          if (local_fm_nodes.size() == 0) {
            return;
          }

          Metrics best_Metrics = {context.dynamic.localFM_round->incremental_km1,
                                  mt_kahypar::metrics::imbalance(*partitioned_hypergraph_s, context)};

          _fm->refine(partitioned_hypergraph, local_fm_nodes, best_Metrics, std::numeric_limits<double>::max());

          for (Move move : context.dynamic.localFM_round->moves) {
            if (move.to != partitioned_hypergraph_s->partID(move.node)) {
              continue;
            }
            _rebalancer.applyMove(move);
          }

          _rebalancer.updateGainForMoves(context.dynamic.localFM_round->moves);

          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains())
        }

        PartitionID add_node_to_partitioned_hypergraph(ds::StaticHypergraph& hypergraph, Context& context, const HypernodeID& hn) {

          //compute for each block the number of nodes the new node connected to
          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_s->connectivitySet(he) ) {
              ASSERT(partitioned_hypergraph_s->checkConnectivitySet(he, context.partition.k));
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
          return std::get<1>(block_connectivities[0]);
        }

    public:

        LocalFMRebalanceVCycle() : _benefit_aggregator(*new vec<Gain>()) {}

        void init(ds::StaticHypergraph& hypergraph, Context& context) override {
          repartition(hypergraph, context);
          init_local_fm(hypergraph, context);
          _rebalancer.init(partitioned_hypergraph_s.value(), context, _gain_cache, _benefit_aggregator);
          prior_total_weight = hypergraph.totalWeight();
          changed_weight = 0;
        }

        void partition(ds::StaticHypergraph& hypergraph, Context& context, Change change, size_t changes_size) override {

          parallel::scalable_vector<HypernodeID> local_fm_nodes;
          std::vector<HypernodeID> gain_cache_nodes;

          ASSERT(partitioned_hypergraph_s->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
          ASSERT(context.dynamic.localFM_round->incremental_km1 == metrics::quality(*partitioned_hypergraph_s, Objective::km1));

          //account for removed edges in incremental km1
//          for (const HyperedgeID& he : change.removed_edges) {
//            context.dynamic.localFM_round->incremental_km1 -= std::max(partitioned_hypergraph_s->connectivity(he) - 1, 0) * partitioned_hypergraph_s->edgeWeight(he);
//          }
//          context.dynamic.localFM_round->incremental_km1 = metrics::quality(*partitioned_hypergraph_s, Objective::km1);

          vec<PartitionID> empty_blocks;

          for (const HypernodeID& hn : change.removed_nodes) {
            changed_weight += hypergraph.nodeWeight(hn);
            empty_blocks.push_back(partitioned_hypergraph_s->partID(hn));
            //TODO: mixed queries -> remove node from local_fm_nodes
            //TODO: does this work for multiple node removals?
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {

              for (const HypernodeID& hn2 : hypergraph.pins(he)) {
                if (partitioned_hypergraph_s->pinCountInPart(he, partitioned_hypergraph_s->partID(hn2)) <= context.dynamic.small_blocks_threshold) {
                  local_fm_nodes.push_back(hn2);
                }
              }

              size_t nodes_in_removed_partition_prior_removal = partitioned_hypergraph_s->pinCountInPart(he, partitioned_hypergraph_s->partID(hn));
              //gain increases for remaining node in partition because moving it to another partition will decrease the cut
              if (nodes_in_removed_partition_prior_removal == 2) {
                //append remaining node in partition to gain_cache_nodes
                for (const HypernodeID& hn2 : hypergraph.pins(he)) {
                  if (hn2 != hn && partitioned_hypergraph_s->partID(hn2) == partitioned_hypergraph_s->partID(hn)) {
                    gain_cache_nodes.push_back(hn2);
                    break;
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

          process_change(hypergraph, context, change);

          for (const HypernodeID& hn : change.added_nodes) {
            changed_weight += hypergraph.nodeWeight(hn);
            PartitionID assigned_part = add_node_to_partitioned_hypergraph(hypergraph, context, hn);
            _rebalancer.insertOrUpdateNode(hn, assigned_part);
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
                // Empiric: is not so sure about this
                gain_cache_nodes.insert(gain_cache_nodes.end(), hypergraph.pins(he).begin(), hypergraph.pins(he).end());
              }
            }
          }

          ASSERT(context.dynamic.localFM_round->incremental_km1 == metrics::quality(*partitioned_hypergraph_s, Objective::km1));

          //reset pin counts of added edges
          for (const HyperedgeID& he : change.added_edges) {
            context.dynamic.localFM_round->incremental_km1 -= std::max(partitioned_hypergraph_s->connectivity(he) - 1, 0) * partitioned_hypergraph_s->edgeWeight(he);
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_s->pinCountInPart(he, p) > 0) {
                partitioned_hypergraph_s->decrementPinCountOfBlockWrapper(he, p);
              }
            }
            for (const HypernodeID& hn : hypergraph.pins(he)) {
              partitioned_hypergraph_s->incrementPinCountOfBlockWrapper(he, partitioned_hypergraph_s->partID(hn));
            }
            context.dynamic.localFM_round->incremental_km1 += std::max(partitioned_hypergraph_s->connectivity(he) - 1, 0) * partitioned_hypergraph_s->edgeWeight(he);
          }

          //remove duplicates
          std::sort(gain_cache_nodes.begin(), gain_cache_nodes.end());
          gain_cache_nodes.erase(std::unique(gain_cache_nodes.begin(), gain_cache_nodes.end()), gain_cache_nodes.end());
          std::sort(local_fm_nodes.begin(), local_fm_nodes.end());
          local_fm_nodes.erase(std::unique(local_fm_nodes.begin(), local_fm_nodes.end()), local_fm_nodes.end());

          local_fm(context, local_fm_nodes, gain_cache_nodes, change, empty_blocks);

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

          ASSERT(metrics::isBalanced(*partitioned_hypergraph_s, context));

          PartitionResult partition_result = *new PartitionResult();
          partition_result.valid = true;
          history.push_back(partition_result);

        }

        void compute_km1_and_imbalance(ds::StaticHypergraph& hypergraph, Context &context, Change change, PartitionResult& partition_result) override {
          partition_result.km1 = mt_kahypar::metrics::quality(*partitioned_hypergraph_s, Objective::km1);
          if (!context.dynamic.server) {
            std::cout << "  km1: " << partition_result.km1 << " Incremental km1: "
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
