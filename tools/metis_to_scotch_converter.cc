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

static void writeScotchGraphFile(const Hypergraph& graph,
                                 const std::string& hgr_filename) {
  std::ofstream out(hgr_filename.c_str());
  out << "0" << std::endl;
  out << graph.initialNumNodes() << " " << ( 2 * graph.initialNumEdges() ) << std::endl;
  out << "0 000" << std::endl; // we only support conversion of unweighted instances here

  for ( const HypernodeID& u : graph.nodes() ) {
    out << graph.nodeDegree(u);
    for ( const HyperedgeID& e : graph.incidentEdges(u) ) {
      for ( const HypernodeID& v : graph.pins(e) ) {
        if ( u != v ) {
          out << " " << v;
        }
      }
    }
    out << std::endl;
  }

  out.close();
}

int main(int argc, char* argv[]) {
  std::string graph_filename;
  std::string out_filename;

  po::options_description options("Options");
  options.add_options()
    ("graph,g",
    po::value<std::string>(&graph_filename)->value_name("<string>")->required(),
    "Metis filename")
    ("out-file,o",
    po::value<std::string>(&out_filename)->value_name("<string>")->required(),
    "Graph Output Filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  Hypergraph graph = mt_kahypar::io::readMetisFile(graph_filename, true);
  writeScotchGraphFile(graph, out_filename);

  return 0;
}
