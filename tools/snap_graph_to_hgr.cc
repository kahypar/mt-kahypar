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
  int num_nodes;
  int num_edges;

  po::options_description options("Options");
  options.add_options()
    ("graph,g",
    po::value<std::string>(&graph_filename)->value_name("<string>")->required(),
    "Graph filename")
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph filename")
    ("num-nodes",
    po::value<int>(&num_nodes)->value_name("<int>")->required(),
    "Number of Nodes")
    ("num-edges",
    po::value<int>(&num_edges)->value_name("<int>")->required(),
    "Number of Edges");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  std::ofstream out_stream(hgr_filename.c_str());

  std::ifstream in_stream(graph_filename);
  std::string line;
  std::getline(in_stream, line);

  // skip any comment
  while (line[0] == '#') {
    std::getline(in_stream, line);
  }

  // Read graph edges
  std::vector<std::pair<int, int> > edges;
  int max_node_id = 0;
  for (int i = 0; i < num_edges; ++i) {
    std::istringstream sstream(line);
    int u, v;
    sstream >> u >> v;
    max_node_id = std::max(max_node_id, std::max(u, v));
    edges.emplace_back(std::make_pair(u, v));
    std::getline(in_stream, line);
  }

  int node_id = 0;
  std::vector<int> node_mapping(max_node_id + 1, -1);
  for (const auto& edge : edges) {
    int u = edge.first;
    int v = edge.second;
    if (node_mapping[u] == -1) {
      node_mapping[u] = ++node_id;
    }
    if (node_mapping[v] == -1) {
      node_mapping[v] = ++node_id;
    }
  }

  // Write header
  out_stream << num_edges << " " << num_nodes << " 0"  /* Unweighted */ << std::endl;

  // Write hyperedges
  for (const auto& edge : edges) {
    int u = edge.first;
    int v = edge.second;
    if (u != v) {
      out_stream << node_mapping[u] << " " << node_mapping[v] << std::endl;
    } else {
      out_stream << node_mapping[u] << std::endl;
    }
  }

  in_stream.close();
  out_stream.close();

  return 0;
}
