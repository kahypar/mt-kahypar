/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <CLI/CLI.hpp>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/connectivity_info.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"


using namespace mt_kahypar;

using Hypergraph = ds::StaticHypergraph;
using PartitionedHypergraph = ds::PartitionedHypergraph<Hypergraph, ds::ConnectivityInfo>;

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
    hypergraph.initializePartition();
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

int main(int argc, char* argv[]) {
  Context context;

  CLI::App app;
  app.set_help_flag("--help");
  app.add_option(
    "-h,--hypergraph,hypergraph",
    context.partition.graph_filename,
    "Hypergraph filename"
  )->required()->check(CLI::ExistingFile);
  app.add_option(
    "-b,--bipart-partition-file,bipart-partition-file",
    context.partition.graph_partition_filename,
    "Partition Filename"
  )->required()->check(CLI::ExistingFile);
  app.add_option(
    "--k,--blocks",
    context.partition.k,
    "Number of Blocks"
  )->required();
  CLI11_PARSE(app, argc, argv);

  // Read Hypergraph
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar::io::readInputFile(
      context.partition.graph_filename, PresetType::default_preset,
      InstanceType::hypergraph, FileFormat::hMetis, true);
  Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);
  PartitionedHypergraph phg(context.partition.k, hg, parallel_tag_t());

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
            << " cut=" << metrics::quality(phg, Objective::cut)
            << " km1=" << metrics::quality(phg, Objective::km1) << std::endl;

  utils::delete_hypergraph(hypergraph);

  return 0;
}
