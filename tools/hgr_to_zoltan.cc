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
#include <string>

#include <CLI/CLI.hpp>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"

using namespace mt_kahypar;

using Hypergraph = ds::StaticHypergraph;

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

  CLI::App app;
  app.set_help_flag("--help");
  app.add_option(
    "-h,--hypergraph,hypergraph",
    hgr_filename,
    "Hypergraph filename"
  )->required()->check(CLI::ExistingFile);
  app.add_option(
    "-o,--out-file,out-file",
    out_filename,
    "Hypergraph Output Filename"
  )->required();
  CLI11_PARSE(app, argc, argv);

  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar::io::readInputFile(
      hgr_filename, PresetType::default_preset,
      InstanceType::hypergraph, FileFormat::hMetis, true);
  Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);

  writeZoltanHypergraph(hg, out_filename);

  utils::delete_hypergraph(hypergraph);

  return 0;
}
