/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Sebastian Schlag <tobias.heuer@kit.edu>
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

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/hypergraph_io.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

void readBipartPartitionFile(const std::string& bipart_partition_file,
                             PartitionedHypergraph& hypergraph,
                             const PartitionID k) {
  ASSERT(!bipart_partition_file.empty(), "No filename for partition file specified");
  std::ifstream file(bipart_partition_file);
  if (file) {
    for ( PartitionID block = 0; block < k; ++block ) {
      std::string line;
      std::getline(file, line);
      std::istringstream line_stream(line);
      HypernodeID hn = 0;
      PartitionID bipart_block = 0;
      line_stream >> bipart_block;
      ASSERT(block == bipart_block - 1);
      while ( line_stream >> hn ) {
        hypergraph.setOnlyNodePart(hn - 1, block);
      }
    }
    hypergraph.initializePartition(0);
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

int main(int argc, char* argv[]) {
  Context context;

  po::options_description options("Options");
  options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&context.partition.graph_filename)->value_name("<string>")->required(),
    "Hypergraph Filename")
    ("bipart-partition-file,b",
    po::value<std::string>(&context.partition.graph_partition_filename)->value_name("<string>")->required(),
    "BiPart Partition Filename")
    ("blocks,k",
    po::value<PartitionID>(&context.partition.k)->value_name("<int>")->required(),
    "Number of Blocks");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  // Read Hypergraph
  Hypergraph hg =
    mt_kahypar::io::readHypergraphFile(context.partition.graph_filename, 0, true);
  PartitionedHypergraph phg(context.partition.k, 0, hg);

  // Setup Context
  context.partition.epsilon = 0.03;
  context.setupPartWeights(hg.totalWeight());

  // Read Bipart Partition File
  readBipartPartitionFile(context.partition.graph_partition_filename, phg,
    context.partition.k);

  std::cout << "RESULT"
            << " graph=" << context.partition.graph_filename
            << " k=" << context.partition.k
            << " imbalance=" << metrics::imbalance(phg, context)
            << " cut=" << metrics::hyperedgeCut(phg)
            << " km1=" << metrics::km1(phg)
            << " soed=" << metrics::soed(phg) << std::endl;
  return 0;
}
