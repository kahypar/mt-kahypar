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
#include <thread>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace io {
namespace {

using Hyperedge = parallel::scalable_vector<HypernodeID>;
using HyperedgeVector = parallel::scalable_vector<Hyperedge>;

static int open_file(const std::string& filename) {
  int fd = open(filename.c_str(), O_RDONLY);
  if ( fd == -1 ) {
    ERROR("Could not open:" << filename);
  }
  return fd;
}

static size_t file_size(int fd) {
  struct stat file_info;
  if ( fstat(fd, &file_info) == -1 ) {
    ERROR("Error while getting file stats");
  }
  return static_cast<size_t>(file_info.st_size);
}

static char* mmap_file(int fd, const size_t length) {
  char* mapped_file = (char*) mmap(0, length, PROT_READ, MAP_SHARED, fd, 0);
  if ( mapped_file == MAP_FAILED ) {
    close(fd);
    ERROR("Error while mapping file to memory");
  }
  return mapped_file;
}

static void munmap_file(char* mapped_file, int fd, const size_t length) {
  if ( munmap(mapped_file, length) == -1 ) {
    close(fd);
    ERROR("Error while unmapping file from memory");
  }
}

static inline void goto_next_line(char* mapped_file, size_t& pos, const size_t length) {
  for ( ; ; ++pos ) {
    if ( pos == length || mapped_file[pos] == '\n' ) {
      ++pos;
      break;
    }
  }
}

static inline int64_t read_number(char* mapped_file, size_t& pos, const size_t length) {
  int64_t number = 0;
  for ( ; pos < length; ++pos ) {
    if ( mapped_file[pos] == ' ' || mapped_file[pos] == '\n' ) {
      while ( mapped_file[pos] == ' ' || mapped_file[pos] == '\n' ) {
        ++pos;
      }
      break;
    }
    ASSERT(mapped_file[pos] >= '0' && mapped_file[pos] <= '9');
    number = number * 10 + (mapped_file[pos] - '0');
  }
  return number;
}

static void readHGRHeader(char* mapped_file,
                          size_t& pos,
                          const size_t length,
                          HyperedgeID& num_hyperedges,
                          HypernodeID& num_hypernodes,
                          mt_kahypar::Type& type) {
  // Skip comments
  while ( mapped_file[pos] == '%' ) {
    goto_next_line(mapped_file, pos, length);
  }

  num_hyperedges = read_number(mapped_file, pos, length);
  num_hypernodes = read_number(mapped_file, pos, length);
  if ( mapped_file[pos - 1] != '\n' ) {
    type = static_cast<mt_kahypar::Type>(read_number(mapped_file, pos, length));
  }
  ASSERT(mapped_file[pos - 1] == '\n');
}

struct HyperedgeRange {
  const size_t start;
  const size_t end;
  const HyperedgeID start_id;
  const HyperedgeID num_hyperedges;
};

static inline bool isSinglePinHyperedge(char* mapped_file,
                                        size_t pos,
                                        const size_t length,
                                        const bool has_hyperedge_weights) {
  size_t num_spaces = 0;
  for ( ; pos < length; ++pos ) {
    if ( mapped_file[pos] == '\n' ) {
      break;
    } else if ( mapped_file[pos] == ' ' ) {
      ++num_spaces;
    }

    if ( num_spaces == 2 ) {
      break;
    }
  }
  return has_hyperedge_weights ? num_spaces == 1 : num_spaces == 0;
}

static HyperedgeID readHyperedges(char* mapped_file,
                                  size_t& pos,
                                  const size_t length,
                                  const HyperedgeID num_hyperedges,
                                  const mt_kahypar::Type type,
                                  HyperedgeVector& hyperedges,
                                  parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight) {
  utils::Timer::instance().start_timer("read_hyperedges", "Read Hyperedges");
  HyperedgeID num_removed_single_pin_hyperedges = 0;
  const bool has_hyperedge_weights = type == mt_kahypar::Type::EdgeWeights ||
                                     type == mt_kahypar::Type::EdgeAndNodeWeights ?
                                     true : false;

  parallel::scalable_vector<HyperedgeRange> hyperedge_ranges;
  tbb::parallel_invoke([&] {
    utils::Timer::instance().start_timer("preprocess_hyperedges", "Preprocess Hyperedges", true);
    // Sequential pass over all hyperedges to determine ranges in the
    // input file that are read in parallel.
    size_t current_range_start = pos;
    HyperedgeID current_range_start_id = 0;
    HyperedgeID current_range_num_hyperedges = 0;
    HyperedgeID current_num_hyperedges = 0;
    const HyperedgeID num_hyperedges_per_range = std::max(
      (num_hyperedges / ( 2 * std::thread::hardware_concurrency())), ID(1));
    while ( current_num_hyperedges < num_hyperedges ) {
      // Skip Comments
      ASSERT(pos < length);
      while ( mapped_file[pos] == '%' ) {
        goto_next_line(mapped_file, pos, length);
        ASSERT(pos < length);
      }

      ASSERT(mapped_file[pos - 1] == '\n');
      if ( !isSinglePinHyperedge(mapped_file, pos, length, has_hyperedge_weights) ) {
        ++current_range_num_hyperedges;
      } else {
        ++num_removed_single_pin_hyperedges;
      }
      ++current_num_hyperedges;
      goto_next_line(mapped_file, pos, length);

      // If there are enough hyperedges in the current scanned range
      // we store that range, which will be later processed in parallel
      if ( current_range_num_hyperedges == num_hyperedges_per_range ) {
        hyperedge_ranges.push_back(HyperedgeRange {
          current_range_start, pos, current_range_start_id, current_range_num_hyperedges});
        current_range_start = pos;
        current_range_start_id += current_range_num_hyperedges;
        current_range_num_hyperedges = 0;
      }
    }
    if ( current_range_num_hyperedges > 0 ) {
      hyperedge_ranges.push_back(HyperedgeRange {
        current_range_start, pos, current_range_start_id, current_range_num_hyperedges});
    }
    utils::Timer::instance().stop_timer("preprocess_hyperedges");
  }, [&] {
    utils::Timer::instance().start_timer("allocate_hyperedges", "Allocate Hyperedges", true);
    hyperedges.resize(num_hyperedges);
    utils::Timer::instance().stop_timer("allocate_hyperedges");
  }, [&] {
    utils::Timer::instance().start_timer("allocate_hyperedge_weights", "Allocate Hyperedge Weights", true);
    hyperedges_weight.assign(num_hyperedges, 1);
    utils::Timer::instance().stop_timer("allocate_hyperedge_weights");
  });

  const HyperedgeID tmp_num_hyperedges = num_hyperedges - num_removed_single_pin_hyperedges;
  hyperedges.resize(tmp_num_hyperedges);
  hyperedges_weight.resize(tmp_num_hyperedges);

  utils::Timer::instance().start_timer("parse_hyperedges", "Parse Hyperedges");
  // Process all ranges in parallel and build hyperedge vector
  tbb::parallel_for(0UL, hyperedge_ranges.size(), [&](const size_t i) {
    HyperedgeRange& range = hyperedge_ranges[i];
    size_t current_pos = range.start;
    const size_t current_end = range.end;
    HyperedgeID current_id = range.start_id;
    const HyperedgeID last_id = current_id + range.num_hyperedges;

    while ( current_id < last_id ) {
      // Skip Comments
      ASSERT(current_pos < current_end);
      while ( mapped_file[pos] == '%' ) {
        goto_next_line(mapped_file, current_pos, current_end);
        ASSERT(current_pos < current_end);
      }

      if ( !isSinglePinHyperedge(mapped_file, current_pos, current_end, has_hyperedge_weights) ) {
        ASSERT(current_id < hyperedges.size());
        if ( has_hyperedge_weights ) {
          hyperedges_weight[current_id] = read_number(mapped_file, current_pos, current_end);
        }

        Hyperedge& hyperedge = hyperedges[current_id];
        // Note, a hyperedge line must contain at least one pin
        HypernodeID pin = read_number(mapped_file, current_pos, current_end);
        ASSERT(pin > 0, V(current_id));
        hyperedge.push_back(pin - 1);
        while ( mapped_file[current_pos - 1] != '\n' ) {
          pin = read_number(mapped_file, current_pos, current_end);
          ASSERT(pin > 0, V(current_id));
          hyperedge.push_back(pin - 1);
        }
        ASSERT(hyperedge.size() >= 2);
        ++current_id;
      } else {
        goto_next_line(mapped_file, current_pos, current_end);
      }
    }
  });
  utils::Timer::instance().stop_timer("parse_hyperedges");

  utils::Timer::instance().stop_timer("read_hyperedges");
  return num_removed_single_pin_hyperedges;
}

static void readHypernodeWeights(char* mapped_file,
                                 size_t& pos,
                                 const size_t length,
                                 const HypernodeID num_hypernodes,
                                 const mt_kahypar::Type type,
                                 parallel::scalable_vector<HypernodeWeight>& hypernodes_weight) {
  utils::Timer::instance().start_timer("parse_hypernode_weights", "Parse Hypernode Weights");
  bool has_hypernode_weights = type == mt_kahypar::Type::NodeWeights ||
                               type == mt_kahypar::Type::EdgeAndNodeWeights ?
                               true : false;
  hypernodes_weight.assign(num_hypernodes, 1);
  if ( has_hypernode_weights ) {
    for ( HypernodeID hn = 0; hn < num_hypernodes; ++hn ) {
      ASSERT(pos > 0 && pos < length);
      ASSERT(mapped_file[pos - 1] == '\n');
      hypernodes_weight[hn] = read_number(mapped_file, pos, length);
    }
  }
  utils::Timer::instance().stop_timer("parse_hypernode_weights");
}

}  // namespace

static inline void readHypergraphFile(const std::string& filename,
                                      HyperedgeID& num_hyperedges,
                                      HypernodeID& num_hypernodes,
                                      HyperedgeID& num_removed_single_pin_hyperedges,
                                      HyperedgeVector& hyperedges,
                                      parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                                      parallel::scalable_vector<HypernodeWeight>& hypernodes_weight) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  mt_kahypar::utils::Timer::instance().start_timer(
    "construct_hypergraph_from_file", "Construct Hypergraph from File");
  int fd = open_file(filename);
  const size_t length = file_size(fd);
  char* mapped_file = mmap_file(fd, length);
  size_t pos = 0;

  // Read Hypergraph Header
  mt_kahypar::Type type = mt_kahypar::Type::Unweighted;
  readHGRHeader(mapped_file, pos, length, num_hyperedges, num_hypernodes, type);

  // Read Hyperedges
  num_removed_single_pin_hyperedges =
    readHyperedges(mapped_file, pos, length, num_hyperedges, type, hyperedges, hyperedges_weight);
  num_hyperedges -= num_removed_single_pin_hyperedges;

  // Read Hypernode Weights
  readHypernodeWeights(mapped_file, pos, length, num_hypernodes, type, hypernodes_weight);
  ASSERT(pos == length);

  munmap_file(mapped_file, fd, length);
  close(fd);
}

static inline Hypergraph readHypergraphFile(const std::string& filename,
                                            const TaskGroupID task_group_id,
                                            const bool stable_construction_of_incident_edges = false) {
  // Read Hypergraph File
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  HyperedgeID num_removed_single_pin_hyperedges = 0;
  HyperedgeVector hyperedges;
  parallel::scalable_vector<HyperedgeWeight> hyperedges_weight;
  parallel::scalable_vector<HypernodeWeight> hypernodes_weight;
  readHypergraphFile(filename, num_hyperedges, num_hypernodes,
    num_removed_single_pin_hyperedges, hyperedges, hyperedges_weight, hypernodes_weight);

  // Construct Hypergraph
  utils::Timer::instance().start_timer("construct_hypergraph", "Construct Hypergraph");
  Hypergraph hypergraph = HypergraphFactory::construct(
          task_group_id, num_hypernodes, num_hyperedges,
          hyperedges, hyperedges_weight.data(), hypernodes_weight.data(),
          stable_construction_of_incident_edges);
  hypergraph.setNumRemovedHyperedges(num_removed_single_pin_hyperedges);
  utils::Timer::instance().stop_timer("construct_hypergraph");
  mt_kahypar::utils::Timer::instance().stop_timer("construct_hypergraph_from_file");
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
      ASSERT(hn < partition.size());
      partition[hn] = hypergraph.partID(hn);
    }
    for (const PartitionID& part : partition) {
      out_stream << part << std::endl;
    }
    out_stream.close();
  }
}
}  // namespace io
}  // namespace mt_kahypar
