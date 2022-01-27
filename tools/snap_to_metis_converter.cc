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
#include <vector>
#include <string>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

using HypernodeID = mt_kahypar::HypernodeID;
using HyperedgeID = mt_kahypar::HyperedgeID;


int main(int argc, char* argv[]) {
  std::string graph_filename;
  std::string out_filename;

  po::options_description options("Options");
  options.add_options()
    ("snap-file,s",
    po::value<std::string>(&graph_filename)->value_name("<string>")->required(),
    "Snap filename")
    ("out-file,o",
    po::value<std::string>(&out_filename)->value_name("<string>")->required(),
    "Graph Output Filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  std::ifstream in_stream(graph_filename);
  std::string line;

  // Read graph edges
  std::vector<std::pair<int, int> > edges;
  std::unordered_map<int, int> node2graph;
  std::vector<int> distinct_nodes;
  int num_nodes = 0;
  while ( std::getline(in_stream, line) ) {
    if ( line[0] == '#' || line[0] == '\n') continue;
    std::istringstream sstream(line);
    int u, v;
    sstream >> u >> v;
    if ( u != v ) {
      if ( node2graph.count(u) == 0 ) {
        node2graph[u] = num_nodes++;
        distinct_nodes.push_back(u);
      }
      if ( node2graph.count(v) == 0 ) {
        node2graph[v] = num_nodes++;
        distinct_nodes.push_back(v);
      }
      edges.emplace_back(std::make_pair(std::min(u,v), std::max(u,v)));
    }
  }
  in_stream.close();

  std::sort(edges.begin(), edges.end());
  std::sort(distinct_nodes.begin(), distinct_nodes.end());
  for ( size_t i = 0; i < distinct_nodes.size(); ++i ) {
    node2graph[distinct_nodes[i]] = i;
  }

  int num_edges = 0;
  std::vector<std::vector<int>> adj_list(num_nodes, std::vector<int>());
  for ( size_t i = 0; i < edges.size(); ++i ) {
    if ( i > 0 && edges[i] == edges[i - 1] ) continue;
    int u = node2graph[edges[i].first];
    int v = node2graph[edges[i].second];
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
    num_edges += 1;
  }


  std::ofstream out(out_filename.c_str());
  out << num_nodes << " " << num_edges << " 0" << std::endl;
  for ( size_t i = 0; i < adj_list.size(); ++i ) {
    std::sort(adj_list[i].begin(), adj_list[i].end());
    for ( const int& v : adj_list[i] ) {
      out << (v + 1) << " ";
    }
    out << std::endl;
  }
  out.close();

  return 0;
}
