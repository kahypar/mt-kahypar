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
#include <string>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

using HypernodeID = mt_kahypar::HypernodeID;
using HyperedgeID = mt_kahypar::HyperedgeID;

static void writeParkwayHypergraphForProc(const Hypergraph& hypergraph,
                                          const std::string& hgr_filename,
                                          const size_t num_procs,
                                          const int rank) {
  const size_t num_vertices = hypergraph.initialNumNodes();
  const size_t num_edges = hypergraph.initialNumEdges();
  const size_t vertices_per_proc = num_vertices / num_procs;
  const size_t hyperedges_per_proc = num_edges / num_procs;

  const HypernodeID hn_start = rank * vertices_per_proc;
  const HypernodeID hn_end = static_cast<size_t>(rank) != num_procs - 1 ?
    (rank + 1) * vertices_per_proc : num_vertices;
  const HyperedgeID he_start = rank * hyperedges_per_proc;
  const HyperedgeID he_end = static_cast<size_t>(rank) != num_procs - 1 ?
    (rank + 1) * hyperedges_per_proc : num_edges;

  std::vector<int> hypernode_weight;
  for ( HypernodeID hn = hn_start; hn < hn_end; ++hn ) {
    hypernode_weight.push_back(hypergraph.nodeWeight(hn));
  }

  std::vector<int> hyperedge_data;
  for ( HyperedgeID he = he_start; he < he_end; ++he ) {
    hyperedge_data.push_back(static_cast<int>(hypergraph.edgeSize(he)) + 2);
    hyperedge_data.push_back(static_cast<int>(hypergraph.edgeWeight(he)));
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      hyperedge_data.push_back(static_cast<int>(pin));
    }
  }

  int num_hypernodes = static_cast<int>(num_vertices);
  int num_local_hypernodes = static_cast<int>(hn_end - hn_start);
  int hyperedge_data_length = static_cast<int>(hyperedge_data.size());

  std::string out_file = hgr_filename + "-" + std::to_string(rank);
  std::ofstream out_stream(out_file, std::ofstream::out | std::ofstream::binary);

  out_stream.write((char *) &num_hypernodes, sizeof(int));
  out_stream.write((char *) &num_local_hypernodes, sizeof(int));
  out_stream.write((char *) &hyperedge_data_length, sizeof(int));
  out_stream.write((char *) hypernode_weight.data(), sizeof(int) * num_local_hypernodes);
  out_stream.write((char *) hyperedge_data.data(), sizeof(int) * hyperedge_data_length);
  out_stream.close();
}

int main(int argc, char* argv[]) {
  std::string hgr_filename;
  std::string out_filename;
  int num_procs;

  po::options_description options("Options");
  options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph filename")
    ("num-procs,p",
    po::value<int>(&num_procs)->value_name("<int>")->required(),
    "Number of Processor Parkway will be called with")
    ("out-file,o",
    po::value<std::string>(&out_filename)->value_name("<string>")->required(),
    "Hypergraph Output Filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  Hypergraph hypergraph =
    mt_kahypar::io::readHypergraphFile(hgr_filename, 0, true);

  for ( int p = 0; p < num_procs; ++p ) {
    writeParkwayHypergraphForProc(hypergraph, out_filename, num_procs, p);
  }

  return 0;
}
