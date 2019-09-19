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
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace io {


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

static inline Hypergraph readHyperedges(std::ifstream& file, 
                                                    const HypernodeID num_hypernodes,
                                                    const HyperedgeID num_hyperedges,
                                                    const mt_kahypar::Type type) {

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

  // Read input file line by line
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::vector<std::string> lines;
  std::string he_line;
  bool has_hyperedge_weights = type == mt_kahypar::Type::EdgeWeights ||
                               type == mt_kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  while ( lines.size() < num_hyperedges ) {
    std::getline(file, he_line);
    lines.push_back(he_line);
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("sequential_read", "Sequential Line Reading",
    "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());


  // Parallel for, for reading hyperedges
  start = std::chrono::high_resolution_clock::now();
  HardwareTopology& topology = HardwareTopology::instance();
  tbb::parallel_for(tbb::blocked_range<size_t>(0UL, lines.size()),
    [&](const tbb::blocked_range<size_t> range) {
    for ( size_t pos = range.begin(); pos < range.end(); ++pos ) {
      std::istringstream line_stream(lines[pos]);

      // Read weight of hyperedge
      std::vector<HypernodeID> hyperedge;
      HyperedgeWeight weight = 1;
      if (has_hyperedge_weights) {
        line_stream >> weight;
      }

      // Read pins of hyperedges
      HypernodeID pin;
      while ( line_stream >> pin ) {
        hyperedge.push_back(--pin);
      }

      // Stream hyperedge into numa hypergraph on which this cpu
      // is part of
      int node = topology.numa_node_of_cpu(sched_getcpu());
      numa_hypergraphs[node].streamHyperedge(hyperedge, weight);
    }
  });
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("stream_hyperedges", "Stream hyperedges",
    "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());

  // Initialize numa hypergraph
  // Involves to memcpy streamed hyperedges of each cpu into
  // global data structure
  start = std::chrono::high_resolution_clock::now();
  for ( int node = 0; node < used_numa_nodes; ++node ) {
    group.run([&, node] {
      TBBNumaArena::instance().numa_task_arena(node).execute([&] {
        numa_hypergraphs[node].initializeHyperedges(num_hypernodes);
      });
    });
  }
  group.wait();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_hyperedges", "Initialize Hyperedges",
    "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());

  start = std::chrono::high_resolution_clock::now();
  Hypergraph hypergraph(num_hypernodes, std::move(numa_hypergraphs));
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_hypernodes", "Initialize Hypernodes",
    "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 3, std::chrono::duration<double>(end - start).count());
  return hypergraph;
}

static inline Hypergraph readHypergraphFile(const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  std::ifstream file(filename);
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  mt_kahypar::Type type = mt_kahypar::Type::Unweighted;
  Hypergraph hypergraph;
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, type);
    hypergraph = readHyperedges(file, num_hypernodes, num_hyperedges, type);
    // TODO(heuer): Read hypernodes weight
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
    exit(EXIT_FAILURE);
  }
  return hypergraph;
}


} // namespace io
} // namespace mt_kahypar
