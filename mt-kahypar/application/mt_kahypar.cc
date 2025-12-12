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
#include <chrono>

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/parallel/thread_management.h"
#include "mt-kahypar/partition/partitioner_facade.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/registries/registry.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/exception.h"

using namespace mt_kahypar;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

int main(int argc, char* argv[]) {

  Context context(false);
  processCommandLineInput(context, argc, argv);

  if ( context.partition.preset_type == PresetType::UNDEFINED ) {
    ERR("No preset specified (--preset-type)");
  }

  // Determine instance (graph or hypergraph) and partition type
  if ( context.partition.instance_type == InstanceType::UNDEFINED ) {
    context.partition.instance_type = to_instance_type(context.partition.file_format);
  }
  context.partition.partition_type = to_partition_c_type(
    context.partition.preset_type, context.partition.instance_type);


  context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
  if (context.partition.enable_logging) {
    io::printBanner();
  }

  utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    utils::Randomize::instance().enableLocalizedParallelShuffle(
      context.shared_memory.shuffle_block_size);
  }

  if constexpr (parallel::provides_hardware_information) {
    size_t num_available_cpus = parallel::num_hardware_cpus();
    if ( num_available_cpus < context.shared_memory.num_threads ) {
      WARNING("There are currently only" << num_available_cpus << "cpus available."
        << "Setting number of threads from" << context.shared_memory.num_threads
        << "to" << num_available_cpus);
      context.shared_memory.num_threads = num_available_cpus;
    }
  }

  // Initialize TBB task arenas on numa nodes
  parallel::initialize_tbb(context.shared_memory.num_threads);

  if constexpr (parallel::provides_hardware_information) {
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    parallel::activate_interleaved_membind_policy();
  }

  // Read Hypergraph
  utils::Timer& timer =
    utils::Utilities::instance().getTimer(context.utility_id);
  timer.start_timer("io_hypergraph", "I/O Hypergraph");
  mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
      context.partition.graph_filename, context.partition.preset_type,
      context.partition.instance_type, context.partition.file_format,
      context.preprocessing.stable_construction_of_incident_edges);
  timer.stop_timer("io_hypergraph");

  // Read Target Graph
  std::unique_ptr<TargetGraph> target_graph;
  if ( context.partition.objective == Objective::steiner_tree ) {
    if ( context.mapping.target_graph_file != "" ) {
      target_graph = std::make_unique<TargetGraph>(
        io::readInputFile<ds::StaticGraph>(
          context.mapping.target_graph_file, FileFormat::Metis, true));
    } else {
      throw InvalidInputException("No target graph file specified (use -g <file> or --target-graph-file=<file>)!");
    }
  }

  if ( context.partition.fixed_vertex_filename != "" ) {
    timer.start_timer("read_fixed_vertices", "Read Fixed Vertex File");
    io::addFixedVerticesFromFile(hypergraph,
      context.partition.fixed_vertex_filename, context.partition.k);
    timer.stop_timer("read_fixed_vertices");
  }

  // Initialize Memory Pool and Algorithm/Policy Registries
  register_memory_pool(hypergraph, context);
  register_algorithms_and_policies();

  // Partition Hypergraph
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
    PartitionerFacade::partition(hypergraph, context, target_graph.get());
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  PartitionerFacade::printPartitioningResults(
    partitioned_hypergraph, context, elapsed_seconds);

  if ( context.partition.sp_process_output ) {
    std::cout << PartitionerFacade::serializeResultLine(
      partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if ( context.partition.csv_output ) {
    std::cout << PartitionerFacade::serializeCSV(
      partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if (context.partition.write_partition_file) {
    PartitionerFacade::writePartitionFile(
      partitioned_hypergraph, context.partition.graph_partition_filename);
  }

  parallel::MemoryPool::instance().free_memory_chunks();
  parallel::terminate_tbb();

  utils::delete_hypergraph(hypergraph);
  utils::delete_partitioned_hypergraph(partitioned_hypergraph);

  return 0;
}
