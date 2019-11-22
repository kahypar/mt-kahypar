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
namespace {
using Hyperedge = std::pair<parallel::scalable_vector<HypernodeID>, HyperedgeWeight>;

static inline Hyperedge readStringAsHyperedge(const std::string& hyperedge_line,
                                              bool has_hyperedge_weights) {
  std::istringstream line_stream(hyperedge_line);

  // Read weight of hyperedge
  parallel::scalable_vector<HypernodeID> hyperedge;
  HyperedgeWeight weight = 1;
  if (has_hyperedge_weights) {
    line_stream >> weight;
  }

  // Read pins of hyperedges
  HypernodeID pin;
  while (line_stream >> pin) {
    hyperedge.push_back(--pin);
  }

  return std::make_pair(hyperedge, weight);
}

template <typename StreamingHyperGraph = StreamingHypergraph,
          typename TBB = TBBNumaArena,
          typename HwTopology = HardwareTopology>
static inline void streamHyperedgesEquallyIntoNumaHypergraphs(const std::vector<std::string>& hyperedge_lines,
                                                              std::vector<StreamingHyperGraph>& numa_hypergraphs,
                                                              const bool has_hyperedge_weights) {
  // In case of equal hyperedge distribution, we split the range of hyperedges into
  // equidistant ranges and assign each range to one numa node. We do this by performing
  // parallel for in the corresponding numa task arena over the resp. range.
  size_t num_numa_hypergraphs = numa_hypergraphs.size();
  size_t num_hyperedges_per_hypergraph = hyperedge_lines.size() / num_numa_hypergraphs;
  TBB::instance().execute_parallel_on_all_numa_nodes([&](const int node) {
          size_t start = node * num_hyperedges_per_hypergraph;
          size_t end = node != (int)num_numa_hypergraphs - 1 ?
                       (node + 1) * num_hyperedges_per_hypergraph : hyperedge_lines.size();
          tbb::parallel_for(tbb::blocked_range<HyperedgeID>(start, end),
                            [&](const tbb::blocked_range<HyperedgeID> range) {
            for (HyperedgeID id = range.begin(); id < range.end(); ++id) {
              Hyperedge hyperedge = readStringAsHyperedge(hyperedge_lines[id], has_hyperedge_weights);
              ASSERT(HwTopology::instance().numa_node_of_cpu(sched_getcpu()) == node);
              numa_hypergraphs[node].streamHyperedge(hyperedge.first, id, hyperedge.second);
            }
          });
        });
}

template <typename StreamingHyperGraph = StreamingHypergraph,
          typename TBB = TBBNumaArena,
          typename HwTopology = HardwareTopology>
static inline void streamHyperedgesRandomIntoNumaHypergraphs(const std::vector<std::string>& hyperedge_lines,
                                                             std::vector<StreamingHyperGraph>& numa_hypergraphs,
                                                             const bool has_hyperedge_weights) {
  // In case of random hyperedge distribution, we perform
  // a parallel for on the global thread pool of TBB
  HwTopology& topology = HwTopology::instance();
  tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, hyperedge_lines.size()),
                    [&](const tbb::blocked_range<HyperedgeID> range) {
          for (HyperedgeID id = range.begin(); id < range.end(); ++id) {
            Hyperedge hyperedge = readStringAsHyperedge(hyperedge_lines[id], has_hyperedge_weights);
            // Stream hyperedge into numa hypergraph on which this cpu
            // is part of
            int node = topology.numa_node_of_cpu(sched_getcpu());
            numa_hypergraphs[node].streamHyperedge(hyperedge.first, id, hyperedge.second);
          }
        });
}

template <typename StreamingHyperGraph = StreamingHypergraph,
          typename TBB = TBBNumaArena,
          typename HwTopology = HardwareTopology>
static inline void streamHyperedgesIntoOneNumaHypergraph(const std::vector<std::string>& hyperedge_lines,
                                                         std::vector<StreamingHyperGraph>& numa_hypergraphs,
                                                         const bool has_hyperedge_weights) {
  // In case all hyperedges should be assigned to one numa node, we perform
  // a parallel for in the numa task arena of numa node 0.
  TBB::instance().numa_task_arena(0).execute([&] {
          tbb::parallel_for(tbb::blocked_range<HyperedgeID>(0UL, hyperedge_lines.size()),
                            [&](const tbb::blocked_range<HyperedgeID> range) {
            for (HyperedgeID id = range.begin(); id < range.end(); ++id) {
              Hyperedge hyperedge = readStringAsHyperedge(hyperedge_lines[id], has_hyperedge_weights);
              ASSERT(HwTopology::instance().numa_node_of_cpu(sched_getcpu()) == 0);
              numa_hypergraphs[0].streamHyperedge(hyperedge.first, id, hyperedge.second);
            }
          });
        });
}

template <typename HyperGraph = Hypergraph,
          typename StreamingHyperGraph = StreamingHypergraph,
          typename TBB = TBBNumaArena,
          typename HwTopology = HardwareTopology>
static inline HyperGraph readHyperedges(std::ifstream& file,
                                        const HypernodeID num_hypernodes,
                                        const HyperedgeID num_hyperedges,
                                        const mt_kahypar::Type type,
                                        const InitialHyperedgeDistribution distribution,
                                        const PartitionID k) {
  // Allocate numa hypergraph on their corresponding numa nodes
  int used_numa_nodes = TBB::instance().num_used_numa_nodes();
  std::vector<StreamingHyperGraph> numa_hypergraphs;
  TBB::instance().execute_sequential_on_all_numa_nodes([&](const int node) {
          numa_hypergraphs.emplace_back(node, k);
        });

  // Read input file line by line
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::vector<std::string> lines;
  std::string he_line;
  bool has_hyperedge_weights = type == mt_kahypar::Type::EdgeWeights ||
                               type == mt_kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  while (lines.size() < num_hyperedges) {
    std::getline(file, he_line);
    lines.push_back(he_line);
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("sequential_read", "Sequential Line Reading",
                                                  "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());

  // Parallel for, for reading hyperedges
  start = std::chrono::high_resolution_clock::now();
  switch (distribution) {
    case InitialHyperedgeDistribution::equally:
      streamHyperedgesEquallyIntoNumaHypergraphs<StreamingHyperGraph, TBB, HwTopology>(
        lines, numa_hypergraphs, has_hyperedge_weights);
      break;
    case InitialHyperedgeDistribution::random:
      streamHyperedgesRandomIntoNumaHypergraphs<StreamingHyperGraph, TBB, HwTopology>(
        lines, numa_hypergraphs, has_hyperedge_weights);
      break;
    case InitialHyperedgeDistribution::all_on_one:
      streamHyperedgesIntoOneNumaHypergraph<StreamingHyperGraph, TBB, HwTopology>(
        lines, numa_hypergraphs, has_hyperedge_weights);
      break;
    default:
      LOG << "Unknown distribution strategy";
      exit(0);
  }
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("stream_hyperedges", "Stream hyperedges",
                                                  "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 1, std::chrono::duration<double>(end - start).count());

  // Initialize numa hypergraph
  // Involves to memcpy streamed hyperedges of each cpu into
  // global data structure
  start = std::chrono::high_resolution_clock::now();
  for (int node = 0; node < used_numa_nodes; ++node) {
    TBB::instance().numa_task_arena(node).execute([&] {
            TBB::instance().numa_task_group(node).run([&, node] {
              numa_hypergraphs[node].initializeHyperedges(num_hypernodes);
            });
          });
  }
  TBB::instance().wait();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_hyperedges", "Initialize Hyperedges",
                                                  "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 2, std::chrono::duration<double>(end - start).count());

  start = std::chrono::high_resolution_clock::now();
  HyperGraph hypergraph(num_hypernodes, std::move(numa_hypergraphs), k);
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_hypernodes", "Initialize Hypernodes",
                                                  "hypergraph_import", mt_kahypar::utils::Timer::Type::IMPORT, 3, std::chrono::duration<double>(end - start).count());
  return hypergraph;
}

template <typename HyperGraph = Hypergraph>
static inline void readHypernodeWeights(std::ifstream& file,
                                        HyperGraph& hypergraph,
                                        const HypernodeID num_hypernodes,
                                        const mt_kahypar::Type type) {
  bool has_hypernode_weights = type == mt_kahypar::Type::NodeWeights ||
                               type == mt_kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  if (has_hypernode_weights) {
    std::string line;
    for (HypernodeID hn = 0; hn < num_hypernodes; ++hn) {
      std::getline(file, line);
      // skip any comments
      while (line[0] == '%') {
        std::getline(file, line);
      }
      std::istringstream line_stream(line);
      HypernodeWeight weight;
      line_stream >> weight;
      hypergraph.setNodeWeight(hypergraph.globalNodeID(hn), weight);
    }
    hypergraph.updateTotalWeight();
  }
}

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
}  // namespace

template <typename HyperGraph = Hypergraph,
          typename StreamingHyperGraph = StreamingHypergraph,
          typename TBB = TBBNumaArena,
          typename HwTopology = HardwareTopology>
static inline HyperGraph readHypergraphFile(const std::string& filename,
                                            const PartitionID k,
                                            const InitialHyperedgeDistribution distribution) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::ifstream file(filename);
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  mt_kahypar::Type type = mt_kahypar::Type::Unweighted;
  HyperGraph hypergraph;
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, type);
    hypergraph = readHyperedges<HyperGraph, StreamingHyperGraph, TBB, HwTopology>(
      file, num_hypernodes, num_hyperedges, type, distribution, k);
    readHypernodeWeights<HyperGraph>(file, hypergraph, num_hypernodes, type);
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
    exit(EXIT_FAILURE);
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("hypergraph_import", "Reading Hypergraph File",
                                                  "", mt_kahypar::utils::Timer::Type::IMPORT, 0, std::chrono::duration<double>(end - start).count());
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

static inline void writePartitionFile(const Hypergraph& hypergraph, const std::string& filename) {
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
}  // namespace io
}  // namespace mt_kahypar
