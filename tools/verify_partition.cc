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

bool readPartitionFile(const std::string& partition_file, PartitionedHypergraph& hypergraph) {
  bool success = true;
  std::vector<PartitionID> partition;
  mt_kahypar::io::readPartitionFile(partition_file, partition);
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    if ( partition[hn] == kInvalidPartition ) {
      LOG << RED << "[ERROR]" << END << "Hypernode" << hn << "is not assigned to a block";
      success = false;
    } else if ( partition[hn] >= hypergraph.k() ) {
      LOG << RED << "[ERROR]" << END << "Hypernode" << hn << "is assigned to block"
          << ( partition[hn] + 1 ) << ", but there are only" << hypergraph.k() << "blocks";
      success = false;
    }
    hypergraph.setOnlyNodePart(hn, partition[hn]);
  }
  hypergraph.initializePartition();
  return success;
}

int main(int argc, char* argv[]) {
  Context context;

  po::options_description options("Options");
  options.add_options()
          ("hypergraph,h",
           po::value<std::string>(&context.partition.graph_filename)->value_name("<string>")->required(),
           "Hypergraph Filename")
          ("partition-file,b",
           po::value<std::string>(&context.partition.graph_partition_filename)->value_name("<string>")->required(),
           "Partition Filename")
          ("blocks,k",
           po::value<PartitionID>(&context.partition.k)->value_name("<int>")->required(),
           "Number of Blocks")
           ("epsilon,e",
           po::value<double>(&context.partition.epsilon)->value_name("<double>")->required(),
           "Imbalance");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  // Read Hypergraph
  Hypergraph hg = mt_kahypar::io::readHypergraphFile(context.partition.graph_filename, true);
  PartitionedHypergraph phg(context.partition.k, hg, parallel_tag_t());

  // Setup Context
  context.setupPartWeights(hg.totalWeight());

  // Read Partition File
  bool success = readPartitionFile(context.partition.graph_partition_filename, phg);

  for ( PartitionID i = 0; i < context.partition.k; ++i ) {
    if ( phg.partWeight(i) == 0 ) {
      LOG << RED << "[ERROR]" << END << "Block" << (i + 1) << "is empty";
      success = false;
    } else if ( phg.partWeight(i) > context.partition.max_part_weights[i] ) {
      LOG << RED << "[ERROR]" << END << "Block" << (i + 1) << "has weight"
          << phg.partWeight(i) << ", but maximum allowed block weight is"
          << context.partition.max_part_weights[i];
      success = false;
    }
  }

  return success ? 0 : -1;
}
