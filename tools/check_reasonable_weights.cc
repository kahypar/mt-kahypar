/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <string>

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

template<typename Hypergraph>
double maxRelativeWeight(const Hypergraph& hg) {
  double max = 0;
  for (HypernodeID hn : hg.nodes()) {
    for (Dimension d = 0; d < hg.dimension(); ++d) {
      max = std::max(max, hg.nodeWeight(hn).at(d) / static_cast<double>(hg.totalWeight().at(d)));
    }
  }
  return max;
}

int main(int argc, char* argv[]) {
  std::string hgr_filename;
  FileFormat file_format = FileFormat::hMetis;
  PartitionID num_blocks;
  double fraction;

  po::options_description options("Options");
  options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph (or graph) filename")
    ("input-file-format",
      po::value<std::string>()->value_name("<string>")->notifier([&](const std::string& s) {
        if (s == "hmetis") {
          file_format = FileFormat::hMetis;
        } else if (s == "metis") {
          file_format = FileFormat::Metis;
        }
      }),
      "Input file format: \n"
      " - hmetis : hMETIS hypergraph file format \n"
      " - metis : METIS graph file format")
    ("num-blocks,k",
    po::value<int>(&num_blocks)->value_name("<int>")->required(),
    "Number of blocks for partitioning")
    ("fraction,f",
    po::value<double>(&fraction)->value_name("<double>")->required(),
    "Tolerable fraction of block weight");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  double max_relative_weight;
  if (file_format == FileFormat::Metis) {
    auto graph = io::readInputFile<ds::StaticGraph>(hgr_filename, file_format, true);
    max_relative_weight = maxRelativeWeight(graph);
  } else {
    auto hypergraph = io::readInputFile<ds::StaticHypergraph>(hgr_filename, FileFormat::hMetis, true);
    max_relative_weight = maxRelativeWeight(hypergraph);
  }

  if (max_relative_weight > fraction * (1.0 / static_cast<double>(num_blocks))) {
    return 1;
  } else {
    return 0;
  }
}
