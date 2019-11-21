/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <iostream>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/application/command_line_options.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/partition/partitioner.h"

#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/randomize.h"

int main(int argc, char* argv[]) {

  mt_kahypar::Context context;
  mt_kahypar::processCommandLineInput(context, argc, argv);
  mt_kahypar::io::printBanner(context);

  mt_kahypar::utils::Randomize::instance().setSeed(context.partition.seed);

  // Initialize TBB task arenas on numa nodes
  mt_kahypar::TBBNumaArena::instance(context.shared_memory.num_threads);

  // Read Hypergraph
  mt_kahypar::Hypergraph hypergraph = mt_kahypar::io::readHypergraphFile(
    context.partition.graph_filename, context.partition.k,
    context.shared_memory.initial_distribution);

  // Partition Hypergraph
  mt_kahypar::HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  mt_kahypar::partition::Partitioner().partition(hypergraph, context);
  mt_kahypar::HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  mt_kahypar::io::printPartitioningResults(hypergraph, context, elapsed_seconds);
  mt_kahypar::io::serializer::serialize(hypergraph, context, elapsed_seconds);
  if ( context.partition.write_partition_file ) {
    mt_kahypar::io::writePartitionFile(hypergraph, context.partition.graph_partition_filename);
  }

  mt_kahypar::TBBNumaArena::instance().terminate();
  return 0;
}
