/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <iostream>
#include <mt-kahypar/partition/metrics.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/partition/partitioner_facade.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/exception.h"

using namespace mt_kahypar;

mt_kahypar_partitioned_hypergraph_t partition_hypergraph_km1(mt_kahypar_hypergraph_t hypergraph, Context& context) {

  // Initialize Memory Pool
  register_memory_pool(hypergraph, context);

  mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
          PartitionerFacade::partition(hypergraph, context, nullptr);

  parallel::MemoryPool::instance().free_memory_chunks();
  TBBInitializer::instance().terminate();

  return partitioned_hypergraph;
}

void repartition_all(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<HypernodeID>& disabled_nodes, double step_size) {
  auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
  size_t re_enabled_nodes = 0;

  std::ofstream output_file(context.dynamic.result_folder + "rep_all_" + std::to_string(context.partition.k) + "k_" + std::to_string(step_size) + "s" + (context.dynamic.use_final_weight ? "_final_weight" : ""));

  //enable nodes in steps
  for ( size_t j = 0; j < disabled_nodes.size(); j += static_cast<size_t>(step_size * hypergraph_s.initialNumNodes()) ) {
    //enable step_size amount of nodes
    while ( re_enabled_nodes <= j && re_enabled_nodes < disabled_nodes.size() ) {
      HypernodeID hn = disabled_nodes[re_enabled_nodes];
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

void repartition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<HypernodeID>& disabled_nodes) {
  for ( const HypernodeID& hn : disabled_nodes ) {
    auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
    hypergraph_s.enableHypernodeWithEdges(hn);
    mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
    auto& partitioned_hypergraph_s =
    utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);
    LOG << " " << std::left << std::setw(20) << "imbalance" << "=" << mt_kahypar::metrics::imbalance(partitioned_hypergraph_s, context);
    utils::delete_partitioned_hypergraph(partitioned_hypergraph);
  }
}

void first_fitting_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<HypernodeID>& disabled_nodes) {
  auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
  mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);
  auto& partitioned_hypergraph_s =
  utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

  std::vector<int> block_sizes(context.partition.k, 0);
  for (PartitionID block = 0; block < context.partition.k; ++block ) {
    block_sizes[block] = partitioned_hypergraph_s.partWeight(block);
    std::cout << "Partition " << block << " has size " << block_sizes[block] << std::endl;
  }

  for ( const HypernodeID& hn : disabled_nodes ) {
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

void highest_connectivity_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<HypernodeID>& disabled_nodes, size_t start_id = 0) {
  auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
  mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph = partition_hypergraph_km1(hypergraph, context);

  auto& partitioned_hypergraph_s =
  utils::cast<ds::PartitionedHypergraph<typename ds::StaticHypergraph>>(partitioned_hypergraph);

  std::ofstream output_file(context.dynamic.result_folder + "/connectivity_" + std::to_string(start_id) + "_" + std::to_string(context.partition.k) + "k" + (context.dynamic.use_final_weight ? "_final_weight" : ""));
  if ( !output_file.is_open() ) {
    throw std::runtime_error("Could not open output file");
  }

  for ( size_t i = start_id; i < disabled_nodes.size(); ++i ) {
    HypernodeID hn = disabled_nodes[i];
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

    for (auto & block_connectivitie : block_connectivities) {
      //try adding to the block
      partitioned_hypergraph_s.setNodePart(hn, std::get<1>(block_connectivitie));
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

void repartition_x_connectivity_partition_strategy(mt_kahypar_hypergraph_t hypergraph, Context& context, const std::vector<HypernodeID>& disabled_nodes, size_t start_id = 0) {
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

  for ( size_t i = start_id; i < disabled_nodes.size(); ++i ) {
    HypernodeID hn = disabled_nodes[i];
    hypergraph_s.enableHypernodeWithEdges(hn);
    if (!context.dynamic.use_final_weight) {
      total_hypergraph_weight += hypergraph_s.nodeWeight(hn);
      hypergraph_s.incrementTotalWeight(hn);

      context.partition.perfect_balance_part_weights.clear();
      context.partition.perfect_balance_part_weights = std::vector<HypernodeWeight>(context.partition.k, ceil(
              total_hypergraph_weight
              / static_cast<double>(context.partition.k)));
      context.partition.max_part_weights.clear();
      context.partition.max_part_weights = std::vector<HypernodeWeight>(context.partition.k, (1 + context.partition.epsilon)
                                                                                             * context.partition.perfect_balance_part_weights[0]);

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

int main(int argc, char* argv[]) {

  std::cout << "Starting dynamic partitioning" << std::endl;

  Context context(false);
  processCommandLineInput(context, argc, argv, nullptr);

  if ( context.partition.preset_file.empty() ) {
    if ( context.partition.preset_type != PresetType::UNDEFINED ) {
      // Only a preset type specified => load according preset
      auto preset_option_list = loadPreset(context.partition.preset_type);
      processCommandLineInput(context, argc, argv, &preset_option_list);
    } else {
      throw InvalidInputException("No preset specified");
    }
  }

  // Determine instance (graph or hypergraph) and partition type
  if ( context.partition.instance_type == InstanceType::UNDEFINED ) {
    context.partition.instance_type = to_instance_type(context.partition.file_format);
  }
  context.partition.partition_type = to_partition_c_type(
          context.partition.preset_type, context.partition.instance_type);


  context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
  if (context.partition.verbose_output) {
    io::printBanner();
  }

  utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    utils::Randomize::instance().enableLocalizedParallelShuffle(
            context.shared_memory.shuffle_block_size);
  }

  size_t num_available_cpus = HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < context.shared_memory.num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
                                       << "Setting number of threads from" << context.shared_memory.num_threads
                                       << "to" << num_available_cpus);
    context.shared_memory.num_threads = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  TBBInitializer::instance(context.shared_memory.num_threads);

  // We set the membind policy to interleaved allocations in order to
  // distribute allocations evenly across NUMA nodes
  hwloc_cpuset_t cpuset = TBBInitializer::instance().used_cpuset();
  parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
  hwloc_bitmap_free(cpuset);

  // Read Hypergraph
  utils::Timer& timer =
          utils::Utilities::instance().getTimer(context.utility_id);
  timer.start_timer("io_hypergraph", "I/O Hypergraph");
  mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
          context.partition.graph_filename, context.partition.preset_type,
          context.partition.instance_type, context.partition.file_format,
          context.preprocessing.stable_construction_of_incident_edges);
  timer.stop_timer("io_hypergraph");

  auto& hypergraph_s = utils::cast<ds::StaticHypergraph>(hypergraph);
  std::vector<HypernodeID> disabling_order(hypergraph_s.initialNumNodes());
  std::iota(disabling_order.begin(), disabling_order.end(), 0);

  // Shuffle the order of the nodes using seed
  utils::Randomize::instance().shuffleVector(disabling_order);

  std::vector<HypernodeID> disabled_nodes;
  size_t start_id = context.dynamic.initial_partitioning_size;
  std::cout << "initial partitioning size: " << start_id << std::endl;

  // Disable all the nodes using seed
  for ( const HypernodeID& hn : disabling_order ) {
    hypergraph_s.disableHypernodeWithEdges(hn);
    if ( !context.dynamic.use_final_weight ) {
      hypergraph_s.decrementTotalWeight(hn);
    }
    disabled_nodes.push_back(hn);
  }

  // re-enable all nodes until start_id

  for ( size_t i = 0; i < start_id; ++i ) {
    HypernodeID hn = disabled_nodes.at(i);
    ASSERT(hn < hypergraph_s.initialNumNodes());
    hypergraph_s.enableHypernodeWithEdges(hn);
    if ( !context.dynamic.use_final_weight ) {
      hypergraph_s.incrementTotalWeight(hn);
    }
  }

  //first_fitting_partition_strategy(hypergraph, context, disabled_nodes);
  // repartition_strategy(hypergraph, context, disabled_nodes);
  // highest_connectivity_partition_strategy(hypergraph, context, disabled_nodes, start_id);
   repartition_all(hypergraph, context, disabling_order, 0.01);
//  repartition_x_connectivity_partition_strategy(hypergraph, context, disabled_nodes, start_id);

  utils::delete_hypergraph(hypergraph);

  return 0;
}