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
#include <vector>
#include <string>

#include <CLI/CLI.hpp>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using namespace mt_kahypar;

using HypernodeID = mt_kahypar::HypernodeID;
using HyperedgeID = mt_kahypar::HyperedgeID;

int main(int argc, char* argv[]) {
  std::string graph_filename;
  std::string out_filename;

  CLI::App app;
  app.add_option(
    "-s,--snap-file,snap-file",
    graph_filename,
    "Snap filename"
  )->required()->check(CLI::ExistingFile);
  app.add_option(
    "-o,--out-file,out-file",
    out_filename,
    "Graph Output Filename"
  )->required();
  CLI11_PARSE(app, argc, argv);

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
