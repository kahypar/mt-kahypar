#pragma once

#include <mt-kahypar/partition/metrics.h>
#include <mt-kahypar/dynamic/dynamic_datastructures.h>

namespace mt_kahypar::dyn {

    mt_kahypar_partitioned_hypergraph_t partition_hypergraph_km1(mt_kahypar_hypergraph_t hypergraph, Context& context) {

      // Initialize Memory Pool
      register_memory_pool(hypergraph, context);

      mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
              PartitionerFacade::partition(hypergraph, context, nullptr);

      parallel::MemoryPool::instance().free_memory_chunks();
      TBBInitializer::instance().terminate();

      return partitioned_hypergraph;
    }

    void repartition_all(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<Change>& changes, double step_size) {
      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
      size_t re_enabled_nodes = 0;

      std::ofstream output_file(context.dynamic.result_folder + "rep_all_" + std::to_string(context.partition.k) + "k_" + std::to_string(step_size) + "s" + (context.dynamic.use_final_weight ? "_final_weight" : ""));

      //enable nodes in steps
      for ( size_t j = 0; j < changes.size(); j += static_cast<size_t>(step_size * hypergraph_s.initialNumNodes()) ) {
        //enable step_size amount of nodes
        while ( re_enabled_nodes <= j && re_enabled_nodes < changes.size() ) {
          HypernodeID hn = changes[re_enabled_nodes].added_nodes[0];
          hypergraph_s.enableHypernodeWithEdges(hn);
          if ( !context.dynamic.use_final_weight ) {
            hypergraph_s.incrementTotalWeight(hn);
          }
          re_enabled_nodes++;
        }

        mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
        auto& partitioned_hypergraph_s =
                utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

        // print "re_enabled_nodes, km1"
        std::cout << re_enabled_nodes << ", " << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1) << std::endl;
        //print "re_enabled_nodes, imbalance"
        std::cout << re_enabled_nodes << ", " << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context) << std::endl;

        //write "re_enabled_nodes, km1" to file
        output_file << re_enabled_nodes << ", " << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1) << std::endl;

        utils::delete_partitioned_hypergraph(partitioned_hypergraph);
      }
    }

    void repartition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<Change>& changes) {
      for ( const Change& change : changes ) {
        const HypernodeID& hn = change.added_nodes[0];
        auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
        hypergraph_s.enableHypernodeWithEdges(hn);
        mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
        auto& partitioned_hypergraph_s =
                utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);
        LOG << " " << std::left << std::setw(20) << "imbalance" << "=" << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
        utils::delete_partitioned_hypergraph(partitioned_hypergraph);
      }
    }

    void first_fitting_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<Change>& changes) {
      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
      mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
      auto& partitioned_hypergraph_s =
              utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

      std::vector<int> block_sizes(context.partition.k, 0);
      for (PartitionID block = 0; block < context.partition.k; ++block ) {
        block_sizes[block] = partitioned_hypergraph_s.partWeight(block);
        std::cout << "Partition " << block << " has size " << block_sizes[block] << std::endl;
      }

      for ( const Change& change : changes ) {
        const HypernodeID& hn = change.added_nodes[0];
        PartitionID block = 0;
        for (PartitionID i = 1; i < context.partition.k; ++i) {
          if ( block_sizes[i] < block_sizes[block] ) {
            block = i;
          }
        }
        partitioned_hypergraph_s.setNodePart(hn, block);
        hypergraph_s.enableHypernodeWithEdges(hn);

        //compute new km1
        LOG << " " << std::left << std::setw(20) << "km1" << "=" << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1);
        //compute new imbalance factor
        LOG << " " << std::left << std::setw(20) << "imbalance" << "=" << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
      }
    }

    void highest_connectivity_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<Change>& changes, size_t start_id = 0) {
      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
      mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);

      auto& partitioned_hypergraph_s =
              utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

      std::ofstream output_file(context.dynamic.result_folder + "/connectivity_" + std::to_string(start_id) + "_" + std::to_string(context.partition.k) + "k" + (context.dynamic.use_final_weight ? "_final_weight" : ""));
      if ( !output_file.is_open() ) {
        throw std::runtime_error("Could not open output file");
      }

      for ( size_t i = start_id; i < changes.size(); ++i ) {
        const HypernodeID& hn = changes[i].added_nodes[0];
        hypergraph_s.enableHypernodeWithEdges(hn);

        std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
        for ( PartitionID p = 0; p < context.partition.k; ++p ) {
          block_connectivities[p] = std::make_tuple(0, p);
        }
        for ( const HyperedgeID& he : hypergraph_s.incidentEdges(hn) ) {
          for ( const PartitionID& p : partitioned_hypergraph_s.connectivitySet(he) ) {
            block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
          }
        }
        std::sort(block_connectivities.begin(), block_connectivities.end());
        std::reverse(block_connectivities.begin(), block_connectivities.end());

        for (auto & block_connectivity : block_connectivities) {
          //try adding to the block
          partitioned_hypergraph_s.setNodePart(hn, std::get<1>(block_connectivity));
          //check if imbalance is still within bounds
          if ( mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context) <= context.partition.epsilon ) {
            break;
          }
          partitioned_hypergraph_s.removeNodePart(hn);
        }

        // quit if no block was found
        if ( partitioned_hypergraph_s.partID(hn) == kInvalidPartition ) {
          std::cout << "No block found for node " << hn << std::endl;
          break;
        }

        //TODO: fix underlying metric
        //ASSERT(mt_kahypar::metrics::isBalanced(partitioned_hypergraph_s, context))

        // write "i, km1" to file
        output_file << i << ", " << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1) << std::endl;

        //compute new km1
        LOG << " " << std::left << std::setw(20) << "km1" << "=" << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1);
        //compute new imbalance factor
        LOG << " " << std::left << std::setw(20) << "imbalance" << "=" << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
      }
    }

    void repartition_x_connectivity_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<Change>& changes, size_t start_id = 0) {
      auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
      mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);

      auto& partitioned_hypergraph_s =
              utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

      HypernodeWeight total_hypergraph_weight = 0;

      if ( !context.dynamic.use_final_weight ) {
        ASSERT(context.partition.use_individual_part_weights == false);
        context.partition.max_part_weights = std::vector<int>(context.partition.k, 0);
        context.partition.perfect_balance_part_weights = std::vector<int>(context.partition.k, 0);

        for ( HypernodeID hn = 0; hn < hypergraph_s.initialNumNodes(); ++hn ) {
          if ( hypergraph_s.nodeIsEnabled(hn) ) {
            total_hypergraph_weight += hypergraph_s.nodeWeight(hn);
          }
        }
      }

      std::ofstream output_file(context.dynamic.result_folder + "connectivity_x_rep_" + std::to_string(start_id) + "_" + std::to_string(context.partition.k) + "k" + (context.dynamic.use_final_weight ? "_final_weight" : ""));
      if ( !output_file.is_open() ) {
        throw std::runtime_error("Could not open output file");
      }

      int repartition_count = 0;

      for ( size_t i = start_id; i < changes.size(); ++i ) {
        const HypernodeID& hn = changes[i].added_nodes[0];
        hypergraph_s.enableHypernodeWithEdges(hn);
        if (!context.dynamic.use_final_weight) {
          total_hypergraph_weight += hypergraph_s.nodeWeight(hn);
          hypergraph_s.incrementTotalWeight(hn);

          context.partition.perfect_balance_part_weights.clear();
          context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
                  total_hypergraph_weight
                  / static_cast<double>(context.partition.k)));
          context.partition.max_part_weights.clear();
          context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, ((1 + context.partition.epsilon)
                                                                                                 * context.partition.perfect_balance_part_weights[0]));

        }

        std::vector<std::tuple<int,int>> block_connectivities(context.partition.k, std::make_tuple(0,0));
        for ( PartitionID p = 0; p < context.partition.k; ++p ) {
          block_connectivities[p] = std::make_tuple(0, p);
        }
        for ( const HyperedgeID& he : hypergraph_s.incidentEdges(hn) ) {
          for ( const PartitionID& p : partitioned_hypergraph_s.connectivitySet(he) ) {
            block_connectivities[p] = std::make_tuple(std::get<0>(block_connectivities[p]) + 1, p);
          }
        }

        auto max_connectivity = std::max_element(block_connectivities.begin(), block_connectivities.end());


        partitioned_hypergraph_s.setNodePart(hn, std::get<1>(*max_connectivity));

        //print total hypergraph weight
        std::cout << "Total hypergraph weight " << total_hypergraph_weight <<  " " << hypergraph_s.totalWeight() << std::endl;

        //print all values relevant for the imbalance calculation
        std::cout << "Partition " << std::get<1>(*max_connectivity) << " has size " << partitioned_hypergraph_s.partWeight(std::get<1>(*max_connectivity)) << std::endl;
        std::cout << "Perfect balance " << context.partition.perfect_balance_part_weights[std::get<1>(*max_connectivity)] << std::endl;
        std::cout << "Imbalance " << partitioned_hypergraph_s.partWeight(std::get<1>(*max_connectivity)) /
                                     static_cast<double>(context.partition.perfect_balance_part_weights[std::get<1>(*max_connectivity)]) - 1.0 << std::endl;


        //check if imbalance is still within bounds else repartition
        if (partitioned_hypergraph_s.partWeight(std::get<1>(*max_connectivity)) /
            static_cast<double>(context.partition.perfect_balance_part_weights[std::get<1>(*max_connectivity)]) - 1.0 > context.partition.epsilon ) {
          partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
          partitioned_hypergraph_s = std::move(utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph));
          LOG << " " << std::left << std::setw(20) << "Repartitioned";
          // print partition sizes after repartition
          for ( PartitionID p = 0; p < context.partition.k; ++p ) {
            std::cout << "Partition " << p << " has size " << partitioned_hypergraph_s.partWeight(p) << std::endl;
          }
          repartition_count++;
        }

        // write "i, km1" to file
        output_file << i << ", " << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1) << std::endl;

        //compute new km1
        LOG << " " << std::left << std::setw(20) << "km1" << "=" << mt_kahypar::metrics::quality(partitioned_hypergraph_s, Objective::km1);
        //compute new imbalance factor
        LOG << " " << std::left << std::setw(20) << "imbalance" << "=" << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
      }

      std::cout << "Repartitioned " << repartition_count << " times" << std::endl;
    }
}

