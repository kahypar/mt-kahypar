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

#include "tbb/pipeline.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"

namespace kahypar {
namespace io {


static inline void readHGRHeader(std::ifstream& file, HyperedgeID& num_hyperedges,
                                 HypernodeID& num_hypernodes, kahypar::Type& hypergraph_type) {
  std::string line;
  std::getline(file, line);

  // skip any comments
  while (line[0] == '%') {
    std::getline(file, line);
  }

  std::istringstream sstream(line);
  int i = 0;
  sstream >> num_hyperedges >> num_hypernodes >> i;
  hypergraph_type = static_cast<kahypar::Type>(i);
}

template< typename HwTopology >
static inline void readHyperedges(std::ifstream& file, 
                                  const HypernodeID num_hypernodes,
                                  const HyperedgeID num_hyperedges,
                                  const kahypar::Type type) {

  // Allocate numa hypergraph on their corresponding numa nodes
  int used_numa_nodes = TBBNumaArena::instance().num_used_numa_nodes();
  std::vector<StreamingHypergraph> numa_hypergraphs;
  tbb::task_group group;
  for ( int node = 0; node < used_numa_nodes; ++node ) {
    TBBNumaArena::instance().numa_task_arena(node).execute([&] {
      group.run([&] {
        numa_hypergraphs.emplace_back(node);
      });
    });
    TBBNumaArena::instance().wait(node, group);
  }

  // We read file sequential before to 
  // fully utilize parallel pipeline afterwards
  std::vector<std::string> lines;
  std::string he_line;
  bool has_hyperedge_weights = type == kahypar::Type::EdgeWeights ||
                               type == kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  while ( lines.size() < num_hyperedges ) {
    std::getline(file, he_line);
    lines.push_back(he_line);
  }

  // Parallel Pipeline
  // 1.) Choose a line containing a hyperedge
  // 2.) Constructs hyperedge vector
  // 3.) Streams it into the corresponding numa hypergraph
  // TODO(heuer): It looks like that this pipeline does not evenly distribute load
  // across the different numa hypergraphs
  using WeightedHyperedge = std::pair<std::vector<HypernodeID>, HyperedgeWeight>;
  size_t num_threads = TBBNumaArena::instance().total_number_of_threads();
  HwTopology& topology = HwTopology::instance();
  std::atomic<size_t> idx(0);
  tbb::parallel_pipeline( 128 * num_threads,
    tbb::make_filter<void, std::string>(
      tbb::filter::parallel,
      [&](tbb::flow_control& fc) -> std::string {
        size_t current_idx = idx++;
        if ( current_idx < num_hyperedges ) {
          return lines[current_idx];
        } else {
          fc.stop();
          return "";
        }
      }
    ) &
    tbb::make_filter<std::string, WeightedHyperedge>(
      tbb::filter::parallel,
      [&](std::string& line) {
        std::istringstream line_stream(line);
        std::vector<HypernodeID> hyperedge;
        HyperedgeWeight weight = 1;
        if (has_hyperedge_weights) {
          line_stream >> weight;
        }
        HypernodeID pin;
        while ( line_stream >> pin ) {
          hyperedge.push_back(--pin);
        }
        return std::make_pair(hyperedge, weight);
      }
    ) &
    tbb::make_filter<WeightedHyperedge, void>(
      tbb::filter::parallel,
      [&](WeightedHyperedge& hyperedge) {
        int node = topology.numa_node_of_cpu(sched_getcpu());
        numa_hypergraphs[node].stream(hyperedge.first, hyperedge.second);
      }
    ) );

  // Initialize numa hypergraph
  // Involves to memcpy streamed hyperedges of each cpu into
  // global data structure
  for ( int node = 0; node < used_numa_nodes; ++node ) {
    TBBNumaArena::instance().numa_task_arena(node).execute([&] {
      numa_hypergraphs[node].initialize(num_hypernodes);
    });
  }
}

template< typename HwTopology >
static inline void readHypergraphFile(const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  std::ifstream file(filename);
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  kahypar::Type type = kahypar::Type::Unweighted;
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, type);
    readHyperedges<HwTopology>(file, num_hypernodes, num_hyperedges, type);
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}


}
} // namespace kahypar
