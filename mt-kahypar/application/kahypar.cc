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


#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/io/csv_output.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"

int main(int argc, char* argv[]) {

  mt_kahypar::Context context(false);
  mt_kahypar::processCommandLineInput(context, argc, argv);
  context.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();
  if (context.partition.verbose_output) {
    mt_kahypar::io::printBanner();
  }

  mt_kahypar::utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    mt_kahypar::utils::Randomize::instance().enableLocalizedParallelShuffle(
      context.shared_memory.shuffle_block_size);
  }

  size_t num_available_cpus = mt_kahypar::HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < context.shared_memory.num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
      << "Setting number of threads from" << context.shared_memory.num_threads
      << "to" << num_available_cpus);
    context.shared_memory.num_threads = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  mt_kahypar::TBBInitializer::instance(context.shared_memory.num_threads);

  // We set the membind policy to interleaved allocations in order to
  // distribute allocations evenly across NUMA nodes
  hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
  mt_kahypar::parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
  hwloc_bitmap_free(cpuset);

  // Read Hypergraph
  mt_kahypar::utils::Timer& timer =
    mt_kahypar::utils::Utilities::instance().getTimer(context.utility_id);
  timer.start_timer("io_hypergraph", "I/O Hypergraph");
  mt_kahypar::Hypergraph hypergraph = mt_kahypar::io::readInputFile(
      context.partition.graph_filename,
      context.partition.file_format,
      context.preprocessing.stable_construction_of_incident_edges);
  timer.stop_timer("io_hypergraph");

  // Initialize Memory Pool
  mt_kahypar::register_memory_pool(hypergraph, context);

  // Partition Hypergraph
  mt_kahypar::HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  mt_kahypar::PartitionedHypergraph partitioned_hypergraph = mt_kahypar::partition(hypergraph, context);
  mt_kahypar::HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  mt_kahypar::io::printPartitioningResults(partitioned_hypergraph, context, elapsed_seconds);


  if ( context.partition.sp_process_output ) {
    std::cout << mt_kahypar::io::serializer::serialize(partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if ( context.partition.csv_output ) {
    std::cout << mt_kahypar::io::csv::serialize(partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if (context.partition.write_partition_file) {
    mt_kahypar::io::writePartitionFile(partitioned_hypergraph, context.partition.graph_partition_filename);
  }

  mt_kahypar::parallel::MemoryPool::instance().free_memory_chunks();
  mt_kahypar::TBBInitializer::instance().terminate();

  return 0;
}
