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
#include <functional>

#include <CLI/CLI.hpp>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/utils/randomize.h"

using namespace mt_kahypar;

int main(int argc, char* argv[]) {
  int N, M, MAX_WEIGHT;
  std::string out_filename;

  CLI::App app;
  app.add_option(
    "-o,--out-file,out-file",
    out_filename,
    "Target Graph Output Filename"
  )->required();
  app.add_option(
    "-n,--rows",
    N,
    "Number of rows"
  )->required();
  app.add_option(
    "-m,--columns",
    M,
    "Number of columns"
  )->required();
  app.add_option(
    "--max-weight",
    MAX_WEIGHT,
    "Maximum weight of an edge in the target graph"
  )->required();
  CLI11_PARSE(app, argc, argv);

  const PartitionID k = N * M;
  out_filename = out_filename + ".k" + std::to_string(k);

  auto up = [&](const PartitionID& u) {
    PartitionID v = u - M;
    return v >= 0 ? v : kInvalidPartition;
  };

  auto down = [&](const PartitionID& u) {
    PartitionID v = u + M;
    return v < k ? v : kInvalidPartition;
  };

  auto left = [&](const PartitionID& u) {
    PartitionID v = u - 1;
    return ( u / M ) == (v / M) ? v : kInvalidPartition;
  };

  auto right = [&](const PartitionID& u) {
    PartitionID v = u + 1;
    return ( u / M ) == (v / M) ? v : kInvalidPartition;
  };

  std::ofstream out(out_filename.c_str());
  int num_nodes = k, num_edges = 0;
  for ( PartitionID u = 0; u < k; ++u ) {
    num_edges += (up(u) != kInvalidPartition);
    num_edges += (right(u) != kInvalidPartition);
    num_edges += (down(u) != kInvalidPartition);
    num_edges += (left(u) != kInvalidPartition);
  }
  out << num_nodes << " " << (num_edges/2) << " 1" << std::endl;

  utils::Randomize& rand = utils::Randomize::instance();
  rand.setSeed(std::hash<std::string>{}(out_filename));
  for ( PartitionID u = 0; u < k; ++u ) {
    std::vector<PartitionID> neighbors = { up(u), right(u), down(u), left(u) };
    std::sort(neighbors.begin(), neighbors.end());
    for ( const PartitionID v : neighbors ) {
      if ( v != kInvalidPartition ) {
        out << (v + 1) << " " << rand.getRandomInt(1, MAX_WEIGHT, 0) << " ";
      }
    }
    out << std::endl;
  }
  out.close();

  return 0;
}
