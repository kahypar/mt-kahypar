#pragma once

#include <mt-kahypar/dynamic/dynamic_strategy.h>
#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>
#include <mt-kahypar/partition/factories.h>
#include <mt-kahypar/partition/refinement/rebalancing/incremental_rebalancer.h>
#include <fstream>

namespace mt_kahypar::dyn {

    class LocalFMRebalanceVCycleV4 : public DynamicStrategy {

    private:
        ds::PartitionedHypergraph<ds::MutableHypergraph> partitioned_hypergraph_m;
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _global_rebalancer;
        vec<Gain>& _benefit_aggregator;
        IncrementalRebalancer _rebalancer;
        HyperedgeWeight prior_total_weight = 0;
        HyperedgeWeight changed_weight = 0;
        size_t change_count = 0;

        void init_local_fm() {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _global_rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancing.algorithm, hypergraph_m.initialNumNodes(), context, _gain_cache);

          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = 1;
          // ASSERT(context.refinement.fm.multitry_rounds == 1, context.refinement.fm.multitry_rounds);

          _benefit_aggregator = vec<Gain>(context.partition.k, 0);

          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_m.initialNumNodes(), hypergraph_m.initialNumEdges(), context, _gain_cache, *_global_rebalancer);

          GainCachePtr::resetGainCache(_gain_cache);
          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(partitioned_hypergraph_m);
          _fm->initialize(partitioned_hypergraph);

          context.dynamic.incremental_km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1);

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
        }

        //use local_fm to refine partitioned_hypergraph_m
        void local_fm(parallel::scalable_vector<HypernodeID> local_fm_nodes, std::vector<HypernodeID> gain_cache_nodes, Change change, vec<PartitionID> empty_blocks) {

          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = utils::partitioned_hg_cast(
                  partitioned_hypergraph_m);

          //update gain cache for all changed nodes
          for (const HypernodeID &hn: gain_cache_nodes) {
            ASSERT(_gain_cache.type == GainPolicy::km1);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCacheEntryForNode(
                    partitioned_hypergraph_m, hn, _benefit_aggregator);
            _rebalancer.insertOrUpdateNode(hn, partitioned_hypergraph_m.partID(hn));
          }

          //pull into emptier blocks
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            auto [gain, moved_nodes] = _rebalancer.pullAndUpdateGainCache(p);
            context.dynamic.incremental_km1 -= gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));

          if (!metrics::isBalanced(partitioned_hypergraph_m, context)) {
            // use rebalancer to rebalance partitioned_hypergraph_m
            auto [gain, moved_nodes] = _rebalancer.rebalanceAndUpdateGainCache();
            context.dynamic.incremental_km1 -= gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }

          if (local_fm_nodes.size() == 0) {
            return;
          }

          Metrics best_Metrics = {context.dynamic.incremental_km1,
                                  mt_kahypar::metrics::imbalance(partitioned_hypergraph_m, context)};
          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph_t = utils::partitioned_hg_cast(partitioned_hypergraph_m);
          _fm->refine(partitioned_hypergraph_t, local_fm_nodes, best_Metrics, std::numeric_limits<double>::max());

          for (Move move : context.dynamic.moves) {
            if (move.to != partitioned_hypergraph_m.partID(move.node)) {
              continue;
            }
            _rebalancer.updateHeapsForMove(move);
          }

          _rebalancer.updateGainForMoves(context.dynamic.moves);

          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains());
        }

        PartitionID add_node_to_partitioned_hypergraph(const HypernodeID& hn) {

          //compute for each block the number of nodes the new node connected to
          std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            block_connectivities[p] = std::make_tuple(0, p);
          }
          for ( const HyperedgeID& he : hypergraph_m.incidentEdges(hn) ) {
            for ( const PartitionID& p : partitioned_hypergraph_m.connectivitySet(he) ) {
              ASSERT(partitioned_hypergraph_m.checkConnectivitySet(he, context.partition.k));
              block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
            }
          }

          // sort block_connectivities in descending order
          std::sort(block_connectivities.begin(), block_connectivities.end(), std::greater<std::tuple<int,int>>());

          //Add node to block with highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <
                context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              partitioned_hypergraph_m.setNodePart(hn, std::get<1>(block_connectivity));
              return std::get<1>(block_connectivity);
            }
          }
          // if no partition could accomodate the node put in the best
          partitioned_hypergraph_m.setNodePart(hn, std::get<1>(block_connectivities[0]));
          return std::get<1>(block_connectivities[0]);
        }

    public:

      LocalFMRebalanceVCycleV4(ds::MutableHypergraph& hypergraph_m, Context& context)
          : DynamicStrategy(hypergraph_m, context), _benefit_aggregator(*new vec<Gain>()) {}
      
        MutablePartitionedHypergraph& init() override {
            partitioned_hypergraph_m = partition_hypergraph_km1(hypergraph_m, context);
          init_local_fm();
          _rebalancer.init(partitioned_hypergraph_m, context, _gain_cache, _benefit_aggregator);
          prior_total_weight = hypergraph_m.totalWeight();
          changed_weight = 0;
          return partitioned_hypergraph_m;
        }

        void partition(Change& change, size_t changes_size) override {

          change_count++;

          parallel::scalable_vector<HypernodeID> local_fm_nodes;
          std::vector<HypernodeID> gain_cache_nodes;

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1));

          vec<PartitionID> empty_blocks;

          for (const HypernodeID& hn : change.removed_nodes) {
            changed_weight += hypergraph_m.nodeWeight(hn);
            empty_blocks.push_back(partitioned_hypergraph_m.partID(hn));
            //TODO: mixed queries -> remove node from local_fm_nodes
            //TODO: does this work for multiple node removals?
            for (const HyperedgeID& he : hypergraph_m.incidentEdges(hn)) {

              size_t nodes_in_removed_partition_prior_removal = 0;

              for (const HypernodeID& hn2 : hypergraph_m.pins(he)) {
                if (partitioned_hypergraph_m.partID(hn2) != kInvalidPartition &&
                    partitioned_hypergraph_m.pinCountInPart(he, partitioned_hypergraph_m.partID(hn2)) <=
                    context.dynamic.small_blocks_threshold) {
                  local_fm_nodes.push_back(hn2);
                }
                if (partitioned_hypergraph_m.partID(hn2) == partitioned_hypergraph_m.partID(hn)) {
                  nodes_in_removed_partition_prior_removal++;
                }
              }

              //gain increases for remaining node in partition because moving it to another partition will decrease the cut
              if (nodes_in_removed_partition_prior_removal == 2) {
                //append remaining node in partition to gain_cache_nodes
                for (const HypernodeID& hn2 : hypergraph_m.pins(he)) {
                  if (hn2 != hn && partitioned_hypergraph_m.partID(hn2) == partitioned_hypergraph_m.partID(hn)) {
                    gain_cache_nodes.push_back(hn2);
                    break;
                  }
                }
              //only negative gains => update gains but do not refine
              } else if (nodes_in_removed_partition_prior_removal == 1) {
                context.dynamic.incremental_km1 -= partitioned_hypergraph_m.connectivity(he) > 1 ? partitioned_hypergraph_m.connectivity(he) * hypergraph_m.edgeWeight(he): 0;
                context.dynamic.incremental_km1 += partitioned_hypergraph_m.connectivity(he) > 1 ? (partitioned_hypergraph_m.connectivity(he) - 1) * (hypergraph_m.edgeWeight(he)) : 0;
                for (const HypernodeID& hn2 : hypergraph_m.pins(he)) {
                  if (hn2 != hn) {
                    gain_cache_nodes.push_back(hn2);
                  }
                }
              }
            }
            partitioned_hypergraph_m.removeNodePart(hn);
          }

          for (const HyperedgeID& he : change.removed_edges) {
            context.dynamic.incremental_km1 -= std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
            for (const HypernodeID& hn : hypergraph_m.pins(he)) {
              gain_cache_nodes.push_back(hn);
            }
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_m.pinCountInPart(he, p) > 0) {
                partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, p);
              }
            }
          }

          for (const PinChange& pin_change : change.removed_pins) {
            //decrement km1 if pin is single pin in partition for this edge prior to removal
            if (partitioned_hypergraph_m.pinCountInPart(pin_change.edge, partitioned_hypergraph_m.partID(pin_change.node)) == 1)
            {
              context.dynamic.incremental_km1 -= partitioned_hypergraph_m.edgeWeight(pin_change.edge);
            }
            partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(pin_change.edge, partitioned_hypergraph_m.partID(pin_change.node));
            gain_cache_nodes.push_back(pin_change.node);
            for (const HypernodeID& hn : hypergraph_m.pins(pin_change.edge)) {
              if (hn != pin_change.node) {
                gain_cache_nodes.push_back(hn);
              }
            }
          }

          process_change(hypergraph_m, context, change);

          for (const HypernodeID& hn : change.added_nodes) {
            changed_weight += hypergraph_m.nodeWeight(hn);
            PartitionID assigned_part = add_node_to_partitioned_hypergraph(hn);
            _rebalancer.insertOrUpdateNode(hn, assigned_part);
            local_fm_nodes.push_back(hn);
            gain_cache_nodes.push_back(hn);
            for (const HyperedgeID& he : hypergraph_m.incidentEdges(hn)) {
              size_t nodes_in_added_partition = partitioned_hypergraph_m.pinCountInPart(he, partitioned_hypergraph_m.partID(hn));
              //gain increase for all connected nodes because they can be moved here without penalty
              if (nodes_in_added_partition == 1) {
                context.dynamic.incremental_km1 -= std::max(partitioned_hypergraph_m.connectivity(he) - 2, 0) * partitioned_hypergraph_m.edgeWeight(he);
                context.dynamic.incremental_km1 += std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
                local_fm_nodes.insert(local_fm_nodes.end(), hypergraph_m.pins(he).begin(), hypergraph_m.pins(he).end());
                gain_cache_nodes.insert(gain_cache_nodes.end(), hypergraph_m.pins(he).begin(), hypergraph_m.pins(he).end());
              //only negative gains => update gains but do not refine
              } else if (nodes_in_added_partition == 2) {
                //TODO: only change gains for node in partition (empiric says this breaks things)
                // Empiric: is not so sure about this
                gain_cache_nodes.insert(gain_cache_nodes.end(), hypergraph_m.pins(he).begin(), hypergraph_m.pins(he).end());
              }
            }
          }

          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1));

          //reset pin counts of added edges
          for (const HyperedgeID& he : change.added_edges) {
            context.dynamic.incremental_km1 -= std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_m.pinCountInPart(he, p) > 0) {
                partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, p);
              }
            }
            for (const HypernodeID& hn : hypergraph_m.pins(he)) {
              partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(he, partitioned_hypergraph_m.partID(hn));
            }
            context.dynamic.incremental_km1 += std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
          }

          //remove duplicates
          std::sort(gain_cache_nodes.begin(), gain_cache_nodes.end());
          gain_cache_nodes.erase(std::unique(gain_cache_nodes.begin(), gain_cache_nodes.end()), gain_cache_nodes.end());
          std::sort(local_fm_nodes.begin(), local_fm_nodes.end());
          local_fm_nodes.erase(std::unique(local_fm_nodes.begin(), local_fm_nodes.end()), local_fm_nodes.end());

          local_fm(local_fm_nodes, gain_cache_nodes, change, empty_blocks);

          if (changed_weight > context.dynamic.step_size_pct * prior_total_weight && change_count <= changes_size * (static_cast<float>(context.dynamic.stop_vcycle_at_pct) / 100)) {
            mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = utils::partitioned_hg_cast(
                    partitioned_hypergraph_m);

            HyperedgeWeight prior_km1 = 0;
            if (!context.dynamic.server) {
              prior_km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1);
            }

            //save partition of nodes
            std::vector<PartitionID> partition_of_nodes(hypergraph_m.initialNumNodes());
            for (HypernodeID hn = 0; hn < hypergraph_m.initialNumNodes(); ++hn) {
              partition_of_nodes[hn] = partitioned_hypergraph_m.partID(hn);
            }

            context.partition.num_vcycles = context.dynamic.vcycle_num;

            context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
            mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph_t = utils::partitioned_hg_cast(partitioned_hypergraph_m);

            if (context.dynamic.vcycle_algorithm == "kway_fm") {
              PartitionerFacade::improve(partitioned_hypergraph_t, context);
            } else if (context.dynamic.vcycle_algorithm == "unconstrained") {
              // TODO PartitionerFacade::improve(partitioned_hypergraph_m, *static_cast<Context*>(context.dynamic.old_context));
              PartitionerFacade::improve(partitioned_hypergraph_t, context);
            } else {
              throw std::runtime_error("Unknown vcycle algorithm: " + context.dynamic.vcycle_algorithm);
            }

            HyperedgeWeight post_km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1);
            context.dynamic.incremental_km1 = post_km1;

            GainCachePtr::resetGainCache(_gain_cache);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCache(partitioned_hypergraph_m);

            if (!context.dynamic.server) {
              std::cout << std::endl << "Improvement: " << prior_km1 - post_km1 << std::endl;
            }

            context.partition.num_vcycles = 0;
            changed_weight = 0;
            prior_total_weight = hypergraph_m.totalWeight();

            for (HypernodeID hn = 0; hn < hypergraph_m.initialNumNodes(); ++hn) {
              if (partitioned_hypergraph_m.partID(hn) != partition_of_nodes[hn]) {
                //TODO why is this worse than reset?
                _rebalancer.insertOrUpdateNode(hn, partitioned_hypergraph_m.partID(hn));
              }
            }
            _rebalancer.reset();
          }

          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
        }

        void printAdditionalFinalStats() override {
          assert(context.dynamic.incremental_km1 == mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) && ("Error: incremental_km1 does not match the quality metric. " + std::to_string(context.dynamic.incremental_km1) + " " + std::to_string(mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1))).c_str());
          std::cout << std::endl << "Final km1: " << context.dynamic.incremental_km1 << " Real km1: " << mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) << std::endl;
        }
    };
}
