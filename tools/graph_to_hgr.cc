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
