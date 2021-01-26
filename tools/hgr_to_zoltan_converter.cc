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

#include "mt-kahypar/io/hypergraph_io.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

static void writeZoltanHypergraph(const Hypergraph& hypergraph,
                                  const std::string& hgr_filename) {
  std::ofstream out_stream(hgr_filename.c_str());
  out_stream << 0;                     // 0-based indexing
  out_stream << " " << hypergraph.initialNumNodes() << " " << hypergraph.initialNumEdges() << " " << hypergraph.initialNumPins();
  out_stream << " " << 3 << std::endl;  // weighting scheme: both edge and node weights

  for (const HyperedgeID& he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      out_stream << pin << " ";
    }
    out_stream << "\n";
  }

  for (const HypernodeID& hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << "\n";
  }
  out_stream << std::endl;
  out_stream.close();
}

int main(int argc, char* argv[]) {
  std::string hgr_filename;
  std::string out_filename;

  po::options_description options("Options");
  options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph filename")
    ("out-file,o",
    po::value<std::string>(&out_filename)->value_name("<string>")->required(),
    "Hypergraph Output Filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  Hypergraph hypergraph =
    mt_kahypar::io::readHypergraphFile(hgr_filename, 0, true);

  writeZoltanHypergraph(hypergraph, out_filename);
  return 0;
}
