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

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  std::string graph_filename;
  std::string hgr_filename;


  po::options_description options("Options");
  options.add_options()
    ("graph,g",
    po::value<std::string>(&graph_filename)->value_name("<string>")->required(),
    "Graph filename")
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  std::ofstream out_stream(hgr_filename.c_str());

  std::ifstream in_stream(graph_filename);
  std::string line;
  std::getline(in_stream, line);

  // Read header
  int num_nodes;
  int num_edges;
  {
    std::stringstream sstream(line);
    sstream >> num_nodes >> num_edges;
  }

  std::vector<std::vector<int>> adj_list(num_nodes + 1);
  int u = 1;
  while ( std::getline(in_stream, line) ) {
    std::istringstream sstream(line);
    int v;
    while ( sstream >> v ) {
      adj_list[u].push_back(v);
    }
    ++u;
  }

  // Write header
  out_stream << num_edges << " " << num_nodes << " 0"  /* Unweighted */ << std::endl;

  // Write hyperedges
  for ( int u = 1; u <= num_nodes; ++u ) {
    for ( const int v : adj_list[u]  ) {
      if ( u < v ) {
        out_stream << u << " " << v << std::endl;
      }
    }
  }

  in_stream.close();
  out_stream.close();

  return 0;
}
