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
        gain_cache_t _gain_cache;
        std::unique_ptr<IRefiner> _fm;
        std::unique_ptr<IRebalancer> _global_rebalancer;
        vec<Gain>& _benefit_aggregator;
        IncrementalRebalancer _rebalancer;
        HyperedgeWeight prior_total_weight = 0;
        HyperedgeWeight changed_weight = 0;
        size_t change_count = 0;
        size_t _fm_num_nodes = 0;
        size_t _fm_num_edges = 0;

        void init_local_fm() {
          _gain_cache = GainCachePtr::constructGainCache(context);
          _global_rebalancer = RebalancerFactory::getInstance().createObject(
                  context.refinement.rebalancing.algorithm, hypergraph_m.initialNumNodes(), context, _gain_cache);
          context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
          context.refinement.fm.multitry_rounds = 1;
          _benefit_aggregator = vec<Gain>(context.partition.k, 0);
          _fm = FMFactory::getInstance().createObject(
                  context.refinement.fm.algorithm,
                  hypergraph_m.initialNumNodes(), hypergraph_m.initialNumEdges(), context, _gain_cache, *_global_rebalancer);
          _fm_num_nodes = hypergraph_m.initialNumNodes();
          _fm_num_edges = hypergraph_m.initialNumEdges();

          GainCachePtr::resetGainCache(_gain_cache, hypergraph_m.initialNumNodes(), context.partition.k);
          mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph = utils::partitioned_hg_cast(partitioned_hypergraph_m);
          _fm->initialize(partitioned_hypergraph);

          context.dynamic.incremental_km1 = mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1);

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
        }

        //use local_fm to refine partitioned_hypergraph_m
        void local_fm(parallel::scalable_vector<HypernodeID> local_fm_nodes, std::vector<HypernodeID> gain_cache_nodes, Change change, const vec<PartitionID>& empty_blocks) {
          (void) change;
          (void) empty_blocks;

          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

          //update gain cache for all changed nodes
          for (const HypernodeID &hn: gain_cache_nodes) {
            ASSERT(_gain_cache.type == GainPolicy::km1);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCacheEntryForNode(
                    partitioned_hypergraph_m, hn, _benefit_aggregator);
            _rebalancer.insertOrUpdateNode(hn);
          }

          auto gain_cache_update_duration = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.gain_cache_update_duration_sum += gain_cache_update_duration;

          ASSERT(_rebalancer.checkBlockQueues());
          // ASSERT(_rebalancer.checkPushQueueGains());
          ASSERT(_rebalancer.checkPullQueueGains());

          start = std::chrono::high_resolution_clock::now();

          //pull into emptier blocks
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            auto [gain, moved_nodes] = _rebalancer.pullAndUpdateGainCache(p);
            context.dynamic.incremental_km1 -= gain;
            context.dynamic.km1_gain_rebalance_pull += gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }

          auto rebalance_pull_duration = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.rebalance_duration_sum_pull += rebalance_pull_duration;

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains());

          start = std::chrono::high_resolution_clock::now();

          if (!metrics::isBalanced(partitioned_hypergraph_m, context)) {
            // use rebalancer to rebalance partitioned_hypergraph_m
            auto [gain, moved_nodes] = _rebalancer.rebalanceAndUpdateGainCache();
            context.dynamic.incremental_km1 -= gain;
            context.dynamic.km1_gain_rebalance_push += gain;
            local_fm_nodes.insert(local_fm_nodes.end(), moved_nodes.begin(), moved_nodes.end());
          }
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));
          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains());

          if (local_fm_nodes.size() == 0) {
            return;
          }

          auto rebalance_push_duration = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.rebalance_duration_sum_push += rebalance_push_duration;

          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));

          start = std::chrono::high_resolution_clock::now();

          Metrics best_Metrics = {context.dynamic.incremental_km1,
                                  mt_kahypar::metrics::imbalance(partitioned_hypergraph_m, context)};
          mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph_t = utils::partitioned_hg_cast(partitioned_hypergraph_m);

          if (_fm_num_nodes < hypergraph_m.initialNumNodes()) {
            // in case nodes were added to the hypergraph, we need to initialize the gain cache entries for these nodes
            _fm = FMFactory::getInstance().createObject(
                    context.refinement.fm.algorithm,
                    hypergraph_m.initialNumNodes() * 2, _fm_num_edges, context, _gain_cache, *_global_rebalancer);
            _fm_num_nodes = hypergraph_m.initialNumNodes() * 2;
          } else if (_fm_num_edges < hypergraph_m.initialNumEdges()) {
            // in case edges were added to the hypergraph, we need to initialize the gain cache entries for these edges
            _fm = FMFactory::getInstance().createObject(
                    context.refinement.fm.algorithm,
                    _fm_num_nodes, hypergraph_m.initialNumEdges() * 2, context, _gain_cache, *_global_rebalancer);
            _fm_num_edges = hypergraph_m.initialNumEdges() * 2;
          }

          _fm->refine(partitioned_hypergraph_t, local_fm_nodes, best_Metrics, std::numeric_limits<double>::max());

          auto local_fm_duration = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.localFM_duration_sum += local_fm_duration;

          start = std::chrono::high_resolution_clock::now();

          for (Move move : context.dynamic.local_fm_round->moves) {
            if (move.to != partitioned_hypergraph_m.partID(move.node)) {
              continue;
            }
            _rebalancer.insertOrUpdateNode(move.node);
          }

          for (Move move : context.dynamic.local_fm_round->moves) {
            if (move.to != partitioned_hypergraph_m.partID(move.node)) {
              continue;
            }
            _rebalancer.updateHeapsForMove(move);
            context.dynamic.incremental_km1 -= move.gain;
            context.dynamic.km1_gain_localFM += move.gain;
          }

          _rebalancer.updateGainForMoves(context.dynamic.local_fm_round->moves);

          auto rebalancer_duration = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.rebalance_duration_sum_push += rebalancer_duration;

          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains());
          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));
        }

        PartitionID add_node_to_partitioned_hypergraph(const HypernodeID& hn) {

          // partitioned_hypergraph_m.addNode(hn, kInvalidPartition);

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

          //Add node to block with the highest connectivity if it doesn't violate max_part_weights (imbalance)
          for (const auto& block_connectivity : block_connectivities) {
            if (partitioned_hypergraph_m.partWeight(std::get<1>(block_connectivity)) + hypergraph_m.nodeWeight(hn) <
                context.partition.max_part_weights[std::get<1>(block_connectivity)]) {
              // partitioned_hypergraph_m.setNodePart(hn, std::get<1>(block_connectivity));
              partitioned_hypergraph_m.addNode(hn, std::get<1>(block_connectivity));
              // add node to gain cache and rebalancer
              ASSERT(_gain_cache.type == GainPolicy::km1);
              GainCachePtr::cast<Km1GainCache>(_gain_cache).addNode(hn);
              return std::get<1>(block_connectivity);
            }
          }
          // if no partition could accomodate the node put in the best
          partitioned_hypergraph_m.addNode(hn, std::get<1>(block_connectivities[0]));
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

          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

          change_count++;

          parallel::scalable_vector<HypernodeID> local_fm_nodes;
          std::vector<HypernodeID> gain_cache_nodes;

          ASSERT(partitioned_hypergraph_m.checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(_gain_cache)));
          ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1));
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

          vec<PartitionID> empty_blocks;

          ASSERT(_rebalancer.checkBlockQueues());
          ASSERT(_rebalancer.checkPullQueueGains());

          for (const auto& [hn, he] : change.removed_pins)
          {
            size_t pin_count_in_part_prior_removal = partitioned_hypergraph_m.pinCountInPart(he, partitioned_hypergraph_m.partID(hn));

            //decrement km1 if pin is single pin in partition for this edge prior to removal
            if (pin_count_in_part_prior_removal == 1 &&
                partitioned_hypergraph_m.connectivity(he) > 1)
            {
              context.dynamic.incremental_km1 -= partitioned_hypergraph_m.edgeWeight(he);
            }
            partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, partitioned_hypergraph_m.partID(hn));
            gain_cache_nodes.push_back(hn);

             if (pin_count_in_part_prior_removal == 2) {
              // decrease penalty for remaining node in partition because moving it to another partition will decrease the cut
              for (const HypernodeID& hn2 : hypergraph_m.pins(he))
              {
                if (hn2 != hn)
                {
                  if (partitioned_hypergraph_m.partID(hn2) == partitioned_hypergraph_m.partID(hn))
                  {
                    GainCachePtr::cast<Km1GainCache>(_gain_cache).changePenalty(hn2, -hypergraph_m.edgeWeight(he));
                    _rebalancer.insertOrUpdateNode(hn2);
                  }
                  // if hn2 is in small block, add to local_fm_nodes
                  if (partitioned_hypergraph_m.partID(hn2) != kInvalidPartition &&
                      partitioned_hypergraph_m.pinCountInPart(he, partitioned_hypergraph_m.partID(hn2)) <=
                      context.dynamic.small_blocks_threshold)
                  {
                    local_fm_nodes.push_back(hn2);
                  }
                }
              }
          } else if (pin_count_in_part_prior_removal == 1)
            {
              // reduce benefit for remaining nodes in other partitions because moving them to this partition will increase the cut
              for (const HypernodeID& hn2 : hypergraph_m.pins(he)) {
                if (hn2 != hn) {
                  GainCachePtr::cast<Km1GainCache>(_gain_cache).changeBenefit(hn2, -hypergraph_m.edgeWeight(he), partitioned_hypergraph_m.partID(hn));
                  // not necessary since it would only reduce gain
                  if (!context.dynamic.lazy_pull_updates)
                  {
                    _rebalancer.insertOrUpdateNode(hn2, partitioned_hypergraph_m.partID(hn2), partitioned_hypergraph_m.partID(hn), -hypergraph_m.edgeWeight(he));
                  }
                }
              }
            }

            hypergraph_m.deletePin(he, hn);
          }

          for (const HypernodeID& hn : change.removed_nodes) {
            changed_weight += hypergraph_m.nodeWeight(hn);
            empty_blocks.push_back(partitioned_hypergraph_m.partID(hn));
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
            hypergraph_m.deleteHypernode(hn);
            updateMaxPartWeight(context, hypergraph_m);
          }

          for (const HyperedgeID& he : change.removed_edges) {

            context.dynamic.incremental_km1 -= std::max(partitioned_hypergraph_m.connectivity(he) - 1, 0) * partitioned_hypergraph_m.edgeWeight(he);
            for (const HypernodeID& hn : hypergraph_m.pins(he)) {
              // ASSERT(false);
              gain_cache_nodes.push_back(hn);
            }
            for (PartitionID p = 0; p < context.partition.k; ++p) {
              while(partitioned_hypergraph_m.pinCountInPart(he, p) > 0) {
                // ASSERT(false);
                partitioned_hypergraph_m.decrementPinCountOfBlockWrapper(he, p);
              }
            }
            hypergraph_m.deleteHyperedge(he);
          }

          for (const HypernodeID& hn : change.added_nodes) {
            const HypernodeID new_hn = hypergraph_m.addHypernode({}, 1);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).addNode(hn);
            ASSERT(hn == new_hn);
            changed_weight += hypergraph_m.nodeWeight(hn);
            updateMaxPartWeight(context, hypergraph_m);
            const PartitionID assigned_part = add_node_to_partitioned_hypergraph(hn);
            ASSERT(assigned_part != kInvalidPartition);
            _rebalancer.insertOrUpdateNode(hn);
            local_fm_nodes.push_back(hn);
            gain_cache_nodes.push_back(hn);
          }

          //reset pin counts of added edges
          for (const HyperedgeID& he : change.added_edges) {
            hypergraph_m.addHyperedge({}, 1);
            partitioned_hypergraph_m.addEdge(he);
          }

          for (const auto& [node, edge] : change.added_pins)
          {
            hypergraph_m.addPin(edge, node);
            local_fm_nodes.push_back(node);
            partitioned_hypergraph_m.incrementPinCountOfBlockWrapper(edge, partitioned_hypergraph_m.partID(node));
            gain_cache_nodes.push_back(node);
            PartitionID part_id = partitioned_hypergraph_m.partID(node);
            HyperedgeWeight edge_weight = partitioned_hypergraph_m.edgeWeight(edge);
            if (partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(node)) == 1) {
              for (const HypernodeID& hn2 : hypergraph_m.pins(edge)) {
                if (hn2 != node) {
                  GainCachePtr::cast<Km1GainCache>(_gain_cache).changeBenefit(hn2, edge_weight, part_id);
                  _rebalancer.insertOrUpdateNode(hn2, partitioned_hypergraph_m.partID(hn2), part_id, edge_weight);
                  // if hn2 is in small block, add to local_fm_nodes
                  if (partitioned_hypergraph_m.partID(hn2) != kInvalidPartition &&
                      partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(hn2)) <=
                      context.dynamic.small_blocks_threshold)
                  {
                    local_fm_nodes.push_back(hn2);
                  }
                }
              }
            } else if (partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(node)) == 2) {
              for (const HypernodeID& hn2 : hypergraph_m.pins(edge)) {
                if (hn2 != node) {
                  // gain_cache_nodes.push_back(hn2);
                  if (partitioned_hypergraph_m.partID(hn2) == partitioned_hypergraph_m.partID(node))
                  {
                    GainCachePtr::cast<Km1GainCache>(_gain_cache).changePenalty(hn2, edge_weight);

                  _rebalancer.insertOrUpdateNode(hn2);
                  }
                  // if hn2 is in small block, add to local_fm_nodes
                  if (partitioned_hypergraph_m.partID(hn2) != kInvalidPartition &&
                      partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(hn2)) <=
                      context.dynamic.small_blocks_threshold)
                  {
                    local_fm_nodes.push_back(hn2);
                  }
                }
              }
            }

            //increment km1 if pin is single pin in partition for this edge after addition
            if (partitioned_hypergraph_m.pinCountInPart(edge, partitioned_hypergraph_m.partID(node)) == 1 &&
                partitioned_hypergraph_m.connectivity(edge) > 1)
            {
              context.dynamic.incremental_km1 += partitioned_hypergraph_m.edgeWeight(edge);
            }
          }

          auto processing_duration_sum = std::chrono::high_resolution_clock::now() - start;
          context.dynamic.processing_duration_sum += processing_duration_sum;

          //remove duplicates
          std::sort(gain_cache_nodes.begin(), gain_cache_nodes.end());
          gain_cache_nodes.erase(std::unique(gain_cache_nodes.begin(), gain_cache_nodes.end()), gain_cache_nodes.end());
          std::sort(local_fm_nodes.begin(), local_fm_nodes.end());
          local_fm_nodes.erase(std::unique(local_fm_nodes.begin(), local_fm_nodes.end()), local_fm_nodes.end());

          // remove all disabled nodes from local_fm_nodes and gain_cache_nodes
          local_fm_nodes.erase(std::remove_if(local_fm_nodes.begin(), local_fm_nodes.end(),
                                              [&](const HypernodeID& hn) { return !hypergraph_m.nodeIsEnabled(hn); }),
                               local_fm_nodes.end());
          gain_cache_nodes.erase(std::remove_if(gain_cache_nodes.begin(), gain_cache_nodes.end(),
                                                [&](const HypernodeID& hn) { return !hypergraph_m.nodeIsEnabled(hn); }),
                                 gain_cache_nodes.end());

          auto sorting_duration_sum = std::chrono::high_resolution_clock::now() - start - processing_duration_sum;
          context.dynamic.sorting_duration_sum += sorting_duration_sum;

        // ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));
        local_fm(local_fm_nodes, gain_cache_nodes, change, empty_blocks);
        ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));

          if (changed_weight > context.dynamic.vcycle_step_size_pct * prior_total_weight && change_count <= changes_size * (static_cast<float>(context.dynamic.stop_vcycle_at_pct) / 100)) {
            // std::cout << "Starting v-cycle " << change_count << "/" << changes_size << " after processing " << changed_weight << " weight changes (" << (100.0 * changed_weight / prior_total_weight) << "% of total weight)" << std::endl;

            HighResClockTimepoint vcycle_start = std::chrono::high_resolution_clock::now();

            HyperedgeWeight prior_km1 = context.dynamic.incremental_km1;

            //save partition of nodes
            std::vector<PartitionID> partition_of_nodes(hypergraph_m.initialNumNodes());
            for (HypernodeID hn = 0; hn < hypergraph_m.initialNumNodes(); ++hn) {
              partition_of_nodes[hn] = partitioned_hypergraph_m.partID(hn);
            }

            context.partition.num_vcycles = context.dynamic.vcycle_num;

            context.refinement.fm.algorithm = FMAlgorithm::kway_fm;
            mt_kahypar_partitioned_hypergraph_t  partitioned_hypergraph_t = utils::partitioned_hg_cast(partitioned_hypergraph_m);

            if (!context.partition.use_individual_part_weights)
            {
              context.partition.max_part_weights.clear();
            }

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

            GainCachePtr::resetGainCache(_gain_cache, hypergraph_m.initialNumNodes(), context.partition.k);
            GainCachePtr::cast<Km1GainCache>(_gain_cache).initializeGainCache(partitioned_hypergraph_m);

            context.dynamic.km1_gain_vcycle += prior_km1 - post_km1;

            if (!context.dynamic.server) {
              std::cout << std::endl << "Improvement: " << prior_km1 - post_km1 << std::endl;
            }

            context.partition.num_vcycles = 0;
            changed_weight = 0;
            prior_total_weight = hypergraph_m.totalWeight();

            // std::cout << "Verifying v-cycle partition for " << hypergraph_m.initialNumNodes() << " nodes." << std::endl;
            for (HypernodeID hn = 0; hn < hypergraph_m.initialNumNodes(); ++hn) {
              // std::cout << "Node " << hn << " in part " << partitioned_hypergraph_m.partID(hn) << " should be in part " << partition_of_nodes[hn] << std::endl;
              if (partitioned_hypergraph_m.partID(hn) != partition_of_nodes[hn]) {
                //TODO why is this worse than reset?
                // std::cout << "Resetting node " << hn << " to part " << partition_of_nodes[hn] << std::endl;
                _rebalancer.insertOrUpdateNode(hn);
                // std::cout << "Before: " << partitioned_hypergraph_m.partID(hn) << " After: ";
              }
            }
            // std::cout << "V-cycle partition verified." << std::endl;
            _rebalancer.reset();
            updateMaxPartWeight(context, hypergraph_m);
            // std::cout << "Finished v-cycle " << change_count << "/" << changes_size << std::endl;
            auto vcycle_duration = std::chrono::high_resolution_clock::now() - vcycle_start;
            context.dynamic.vcycle_duration_sum += vcycle_duration;
          }

        ASSERT(context.dynamic.incremental_km1 == metrics::quality(partitioned_hypergraph_m, Objective::km1), context.dynamic.incremental_km1 << " vs. " << metrics::quality(partitioned_hypergraph_m, Objective::km1));
          ASSERT(metrics::isBalanced(partitioned_hypergraph_m, context));

        }

        void printAdditionalFinalStats() override {
          assert(context.dynamic.incremental_km1 == mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) && ("Error: incremental_km1 does not match the quality metric. " + std::to_string(context.dynamic.incremental_km1) + " " + std::to_string(mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1))).c_str());
          std::cout << std::endl << "Final km1: " << context.dynamic.incremental_km1 << " Real km1: " << mt_kahypar::metrics::quality(partitioned_hypergraph_m, Objective::km1) << std::endl;
        }
    };
}
