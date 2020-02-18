/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace tmp_io {
namespace {

using Hyperedge = parallel::scalable_vector<HypernodeID>;
using HyperedgeVector = parallel::scalable_vector<Hyperedge>;

static inline void readHGRHeader(std::ifstream& file, HyperedgeID& num_hyperedges,
                                 HypernodeID& num_hypernodes, mt_kahypar::Type& hypergraph_type) {
  std::string line;
  std::getline(file, line);

  // skip any comments
  while (line[0] == '%') {
    std::getline(file, line);
  }

  std::istringstream sstream(line);
  int i = 0;
  sstream >> num_hyperedges >> num_hypernodes >> i;
  hypergraph_type = static_cast<mt_kahypar::Type>(i);
}

static inline void readStringAsHyperedge(const std::string& hyperedge_line,
                                         const bool has_hyperedge_weights,
                                         Hyperedge& hyperedge,
                                         HyperedgeWeight& hyperedge_weight) {
  std::istringstream line_stream(hyperedge_line);

  // Read weight of hyperedge
  hyperedge_weight = 1;
  if (has_hyperedge_weights) {
    line_stream >> hyperedge_weight;
  }

  // Read pins of hyperedges
  HypernodeID pin;
  while (line_stream >> pin) {
    hyperedge.push_back(--pin);
  }
}

static inline void readHyperedges(std::ifstream& file,
                                  const HyperedgeID num_hyperedges,
                                  const mt_kahypar::Type type,
                                  HyperedgeVector& hyperedges,
                                  parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight) {
  // Read input file line by line
  utils::Timer::instance().start_timer("parse_hyperedges", "Parse Hyperedges");
  parallel::scalable_vector<std::string> lines;
  tbb::parallel_invoke([&] {
    lines.reserve(num_hyperedges);
    std::string he_line;
    while (lines.size() < num_hyperedges) {
      std::getline(file, he_line);
      // skip any comments
      while (he_line[0] == '%') {
        std::getline(file, he_line);
      }
      lines.push_back(he_line);
    }
  }, [&] {
    hyperedges.resize(num_hyperedges);
    hyperedges_weight.resize(num_hyperedges);
  });

  ASSERT(lines.size() == num_hyperedges);
  const bool has_hyperedge_weights = type == mt_kahypar::Type::EdgeWeights ||
                                     type == mt_kahypar::Type::EdgeAndNodeWeights ?
                                     true : false;
  tbb::parallel_for(0UL, num_hyperedges, [&](const HyperedgeID i) {
    readStringAsHyperedge(lines[i], has_hyperedge_weights,
      hyperedges[i], hyperedges_weight[i]);
  });
  utils::Timer::instance().stop_timer("parse_hyperedges");
}

static inline void readHypernodeWeights(std::ifstream& file,
                                        const HypernodeID num_hypernodes,
                                        const mt_kahypar::Type type,
                                        parallel::scalable_vector<HypernodeWeight>& hypernodes_weight) {
  utils::Timer::instance().start_timer("parse_hypernode_weights", "Parse Hypernode Weights");
  bool has_hypernode_weights = type == mt_kahypar::Type::NodeWeights ||
                               type == mt_kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  hypernodes_weight.assign(num_hypernodes, 1);
  if (has_hypernode_weights) {
    std::string line;
    for (HypernodeID hn = 0; hn < num_hypernodes; ++hn) {
      std::getline(file, line);
      // skip any comments
      while (line[0] == '%') {
        std::getline(file, line);
      }
      std::istringstream line_stream(line);
      line_stream >> hypernodes_weight[hn];
    }
  }

  utils::Timer::instance().stop_timer("parse_hypernode_weights");
}

}  // namespace

template <typename HyperGraph,
          typename HyperGraphFactory>
static inline HyperGraph readHypergraphFile(const std::string& filename,
                                            const TaskGroupID task_group_id) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  utils::Timer::instance().start_timer("construct_hypergraph_from_file", "Construct Hypergraph from File");

  // Read Hypergraph file
  std::ifstream file(filename);
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  HyperedgeVector hyperedges;
  parallel::scalable_vector<HyperedgeWeight> hyperedges_weight;
  parallel::scalable_vector<HypernodeWeight> hypernodes_weight;
  mt_kahypar::Type type = mt_kahypar::Type::Unweighted;
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, type);
    readHyperedges(file, num_hyperedges, type, hyperedges, hyperedges_weight);
    readHypernodeWeights(file, num_hypernodes, type, hypernodes_weight);
    file.close();
  } else {
    ERROR("Error: File not found: " + filename);
  }

  // Construct Hypergraph
  utils::Timer::instance().start_timer("construct_hypergraph", "Construct Hypergraph");
  HyperGraph hypergraph = HyperGraphFactory::construct(task_group_id,
    num_hypernodes, num_hyperedges, hyperedges,
    hyperedges_weight.data(), hypernodes_weight.data());
  utils::Timer::instance().stop_timer("construct_hypergraph");
  utils::Timer::instance().stop_timer("construct_hypergraph_from_file");
  return hypergraph;
}

static inline void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  ASSERT(partition.empty(), "Partition vector is not empty");
  std::ifstream file(filename);
  if (file) {
    int part;
    while (file >> part) {
      partition.push_back(part);
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

template<typename HyperGraph>
static inline void writePartitionFile(const HyperGraph& hypergraph, const std::string& filename) {
  if (filename.empty()) {
    LOG << "No filename for partition file specified";
  } else {
    std::ofstream out_stream(filename.c_str());
    std::vector<PartitionID> partition(hypergraph.initialNumNodes(), -1);
    for (const HypernodeID& hn : hypergraph.nodes()) {
      ASSERT(hypergraph.originalNodeID(hn) < partition.size());
      partition[hypergraph.originalNodeID(hn)] = hypergraph.partID(hn);
    }
    for (const PartitionID& part : partition) {
      out_stream << part << std::endl;
    }
    out_stream.close();
  }
}
}  // namespace tmp_io
}  // namespace mt_kahypar
