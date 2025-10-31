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

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

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
bool hasNodeWeights(const Hypergraph& hg) {
  for (HypernodeID hn = 0; hn < hg.initialNumNodes(); ++hn) {
    if (hg.nodeWeight(hn) != weight::broadcast(1, hg.dimension())) {
      return true;
    }
  }
  return false;
}

template<typename Hypergraph>
bool hasHyperedgeWeights(const Hypergraph& hg) {
  for (HyperedgeID he = 0; he < hg.initialNumEdges(); ++he) {
    if (hg.edgeWeight(he) != 1) {
      return true;
    }
  }
  return false;
}

template<typename Hypergraph>
HypernodeWeightArray generateWeights(const Hypergraph& hg, bool natural_weights, bool unit_weights, bool degree_weights) {
  if (natural_weights && unit_weights && !hasNodeWeights(hg)) {
    WARNING("Natural weights are equal to unit weights for unweighted graph! Continuing without natural weights.");
    natural_weights = false;
  }

  Dimension dimension = 0;
  dimension += natural_weights * hg.dimension();
  dimension += unit_weights;
  dimension += degree_weights;

  if (dimension == 0) {
    ERROR("No weights specified!");
  }

  HypernodeWeightArray result(hg.initialNumNodes(), dimension, 0, true);
  for (HypernodeID hn = 0; hn < hg.initialNumNodes(); ++hn) {
    Dimension d = 0;
    if (natural_weights) {
      for (Dimension j = 0; j < hg.dimension(); ++j) {
        result[hn].set(j, hg.nodeWeight(hn).at(j));
      }
      d += hg.dimension();
    }
    if (unit_weights) {
      result[hn].set(d, 1);
      ++d;
    }
    if (degree_weights) {
      result[hn].set(d, hg.nodeDegree(hn));
      ++d;
    }
  }
  return result;
}

void writeHMetisOutput(std::ofstream& out_stream,
                       const ds::StaticHypergraph& hg,
                       const HypernodeWeightArray& hn_weights) {
  bool has_edge_weights = hasHyperedgeWeights(hg);

  // Write header
  out_stream << hg.initialNumEdges() << " " << hg.initialNumNodes() << " ";
  if (hn_weights.dimension() > 1) {
    out_stream << "1";
  }
  out_stream << "1";
  out_stream << (has_edge_weights ? "1" : "0");
  if (hn_weights.dimension() > 1) {
    out_stream << " " << hn_weights.dimension();
  }
  out_stream << std::endl;

  // Write data
  for (HyperedgeID he = 0; he < hg.initialNumEdges(); ++he) {
    if (has_edge_weights) {
      out_stream << hg.edgeWeight(he) << " ";
    }
    for (HypernodeID pin: hg.pins(he)) {
      out_stream << (pin + 1) << " ";
    }
    out_stream << std::endl;
  }
  for (HypernodeID hn = 0; hn < hg.initialNumNodes(); ++hn) {
    for (Dimension d = 0; d < hn_weights.dimension(); ++d) {
      out_stream << hn_weights[hn].at(d) << " ";
    }
    out_stream << std::endl;
  }
  out_stream.close();
}

void writeMetisOutput(std::ofstream& out_stream,
                      const ds::StaticGraph& graph,
                      const HypernodeWeightArray& hn_weights) {
  bool has_edge_weights = hasHyperedgeWeights(graph);

  // Write header
  out_stream << graph.initialNumNodes() << " " << (graph.initialNumEdges() / 2) << " ";
  if (hn_weights.dimension() > 1) {
    out_stream << "1";
  }
  out_stream << "1";
  out_stream << (has_edge_weights ? "1" : "0");
  if (hn_weights.dimension() > 1) {
    out_stream << " " << hn_weights.dimension();
  }
  out_stream << std::endl;

  // Write data
  for (HypernodeID hn = 0; hn < graph.initialNumNodes(); ++hn) {
    for (Dimension d = 0; d < hn_weights.dimension(); ++d) {
      out_stream << hn_weights[hn].at(d) << " ";
    }
    for (HyperedgeID he: graph.incidentEdges(hn)) {
      out_stream << (graph.edgeTarget(he) + 1) << " ";
      if (has_edge_weights) {
        out_stream << graph.edgeWeight(he) << " ";
      }
    }
    out_stream << std::endl;
  }
  out_stream.close();
}

int main(int argc, char* argv[]) {
  std::string hgr_filename;
  std::string out_filename;
  bool is_metis = false;

  bool natural_weights = false;
  bool unit_weights = false;
  bool degree_weights = false;

  po::options_description options("Options");
  options.add_options()
    ("hypergraph,h",
    po::value<std::string>(&hgr_filename)->value_name("<string>")->required(),
    "Hypergraph (or graph) filename")
    ("out,o",
    po::value<std::string>(&out_filename)->value_name("<string>")->required(),
    "Output filename")
    ("is-metis,m",
    po::value<bool>(&is_metis)->value_name("<bool>"),
    "Set to true if using graphs in METIS format")
    ("natural-weights,n",
    po::value<bool>(&natural_weights)->value_name("<bool>"), "Use natural weights.")
    ("unit-weights,u",
    po::value<bool>(&unit_weights)->value_name("<bool>"), "Use unit weights.")
    ("degree-weights,d",
    po::value<bool>(&degree_weights)->value_name("<bool>"), "Use degree weights.");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);
  po::notify(cmd_vm);

  std::ofstream out_stream(out_filename.c_str());

  // Read Hypergraph
  if (is_metis) {
    auto graph = io::readInputFile<ds::StaticGraph>(hgr_filename, FileFormat::Metis, true);
    HypernodeWeightArray hn_weights = generateWeights(graph, natural_weights, unit_weights, degree_weights);
    writeMetisOutput(out_stream, graph, hn_weights);
  } else {
    auto hypergraph = io::readInputFile<ds::StaticHypergraph>(hgr_filename, FileFormat::hMetis, true);
    HypernodeWeightArray hn_weights = generateWeights(hypergraph, natural_weights, unit_weights, degree_weights);
    writeHMetisOutput(out_stream, hypergraph, hn_weights);
  }

  return 0;
}
