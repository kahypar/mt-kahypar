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

#include <string>

#include <CLI/CLI.hpp>

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

using namespace mt_kahypar;

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

  CLI::App app;
  app.set_help_flag("--help");
  app.add_option(
    "-h,--hypergraph,hypergraph",
    hgr_filename,
    "Hypergraph (or graph) filename"
  )->required()->check(CLI::ExistingFile);
  app.add_option(
    "-k,--blocks",
    num_blocks,
    "Number of blocks for partitioning"
  )->required();
  app.add_option(
    "-f,--fraction",
    fraction,
    "Tolerable fraction of block weight"
  )->required();
  app.add_option_function<std::string>(
    "--file-format,--input-file-format",
    [&](const std::string& s) {
      file_format = fileFormatFromString(s);
    },
    "Input file format:\n"
    " - hmetis: hMETIS hypergraph file format\n"
    " - metis: METIS graph file format"
  )->check(CLI::IsMember({"hmetis", "metis"}))->default_str("hmetis");
  CLI11_PARSE(app, argc, argv);

  double max_relative_weight;
  if (file_format == FileFormat::Metis) {
    auto graph = io::readInputFile<ds::StaticGraph>(hgr_filename, file_format, true, false, false);
    max_relative_weight = maxRelativeWeight(graph);
  } else {
    auto hypergraph = io::readInputFile<ds::StaticHypergraph>(hgr_filename, FileFormat::hMetis, true, false, false);
    max_relative_weight = maxRelativeWeight(hypergraph);
  }

  if (max_relative_weight > fraction * (1.0 / static_cast<double>(num_blocks))) {
    return 1;
  } else {
    return 0;
  }
}
