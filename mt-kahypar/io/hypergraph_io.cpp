/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "hypergraph_io.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <thread>
#include <memory>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "tbb/parallel_for.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::io {

  int open_file(const std::string& filename) {
    int fd = open(filename.c_str(), O_RDONLY);
    if ( fd == -1 ) {
      ERROR("Could not open:" << filename);
    }
    return fd;
  }

  size_t file_size(int fd) {
    struct stat file_info;
    if ( fstat(fd, &file_info) == -1 ) {
      ERROR("Error while getting file stats");
    }
    return static_cast<size_t>(file_info.st_size);
  }

  char* mmap_file(int fd, const size_t length) {
    char* mapped_file = (char*) mmap(0, length, PROT_READ, MAP_SHARED, fd, 0);
    if ( mapped_file == MAP_FAILED ) {
      close(fd);
      ERROR("Error while mapping file to memory");
    }
    return mapped_file;
  }

  void munmap_file(char* mapped_file, int fd, const size_t length) {
    if ( munmap(mapped_file, length) == -1 ) {
      close(fd);
      ERROR("Error while unmapping file from memory");
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void goto_next_line(char* mapped_file, size_t& pos, const size_t length) {
    for ( ; ; ++pos ) {
      if ( pos == length || mapped_file[pos] == '\n' ) {
        ++pos;
        break;
      }
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool is_line_ending(char* mapped_file, size_t& pos) {
    return mapped_file[pos] == '\n' || mapped_file[pos] == '\0';
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void do_line_ending(char* mapped_file, size_t& pos) {
    ASSERT(is_line_ending(mapped_file, pos));
    if (mapped_file[pos] != '\0') {
      ++pos;
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  int64_t read_number(char* mapped_file, size_t& pos, const size_t length) {
    int64_t number = 0;
    while ( mapped_file[pos] == ' ' ) {
      ++pos;
    }
    for ( ; pos < length; ++pos ) {
      if ( mapped_file[pos] == ' ' || is_line_ending(mapped_file, pos) ) {
        while ( mapped_file[pos] == ' ' ) {
          ++pos;
        }
        break;
      }
      ASSERT(mapped_file[pos] >= '0' && mapped_file[pos] <= '9');
      number = number * 10 + (mapped_file[pos] - '0');
    }
    return number;
  }

  void readHGRHeader(char* mapped_file,
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
    if ( mapped_file[pos] != '\n' ) {
      type = static_cast<mt_kahypar::Type>(read_number(mapped_file, pos, length));
    }
    do_line_ending(mapped_file, pos);
  }

  struct HyperedgeRange {
    const size_t start;
    const size_t end;
    const HyperedgeID start_id;
    const HyperedgeID num_hyperedges;
  };

  inline bool isSinglePinHyperedge(char* mapped_file,
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

  HyperedgeID readHyperedges(char* mapped_file,
                             size_t& pos,
                             const size_t length,
                             const HyperedgeID num_hyperedges,
                             const mt_kahypar::Type type,
                             HyperedgeVector& hyperedges,
                             parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                             const bool remove_single_pin_hes) {
    HyperedgeID num_removed_single_pin_hyperedges = 0;
    const bool has_hyperedge_weights = type == mt_kahypar::Type::EdgeWeights ||
                                       type == mt_kahypar::Type::EdgeAndNodeWeights ?
                                       true : false;

    parallel::scalable_vector<HyperedgeRange> hyperedge_ranges;
    tbb::parallel_invoke([&] {
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
        if ( !remove_single_pin_hes || !isSinglePinHyperedge(mapped_file, pos, length, has_hyperedge_weights) ) {
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
    }, [&] {
      hyperedges.resize(num_hyperedges);
    }, [&] {
      if ( has_hyperedge_weights ) {
        hyperedges_weight.resize(num_hyperedges);
      }
    });

    const HyperedgeID tmp_num_hyperedges = num_hyperedges - num_removed_single_pin_hyperedges;
    hyperedges.resize(tmp_num_hyperedges);
    if ( has_hyperedge_weights ) {
      hyperedges_weight.resize(tmp_num_hyperedges);
    }

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

        if ( !remove_single_pin_hes || !isSinglePinHyperedge(mapped_file, current_pos, current_end, has_hyperedge_weights) ) {
          ASSERT(current_id < hyperedges.size());
          if ( has_hyperedge_weights ) {
            hyperedges_weight[current_id] = read_number(mapped_file, current_pos, current_end);
          }

          Hyperedge& hyperedge = hyperedges[current_id];
          // Note, a hyperedge line must contain at least one pin
          HypernodeID pin = read_number(mapped_file, current_pos, current_end);
          ASSERT(pin > 0, V(current_id));
          hyperedge.push_back(pin - 1);
          while ( mapped_file[current_pos] != '\n' ) {
            pin = read_number(mapped_file, current_pos, current_end);
            ASSERT(pin > 0, V(current_id));
            hyperedge.push_back(pin - 1);
          }
          do_line_ending(mapped_file, current_pos);
          ASSERT(hyperedge.size() >= 2);
          ++current_id;
        } else {
          goto_next_line(mapped_file, current_pos, current_end);
        }
      }
    });
    return num_removed_single_pin_hyperedges;
  }

  void readHypernodeWeights(char* mapped_file,
                            size_t& pos,
                            const size_t length,
                            const HypernodeID num_hypernodes,
                            const mt_kahypar::Type type,
                            parallel::scalable_vector<HypernodeWeight>& hypernodes_weight) {
    bool has_hypernode_weights = type == mt_kahypar::Type::NodeWeights ||
                                 type == mt_kahypar::Type::EdgeAndNodeWeights ?
                                 true : false;
    if ( has_hypernode_weights ) {
      hypernodes_weight.resize(num_hypernodes);
      for ( HypernodeID hn = 0; hn < num_hypernodes; ++hn ) {
        ASSERT(pos > 0 && pos < length);
        ASSERT(mapped_file[pos - 1] == '\n');
        hypernodes_weight[hn] = read_number(mapped_file, pos, length);
        do_line_ending(mapped_file, pos);
      }
    }
  }


  void readHypergraphFile(const std::string& filename,
                          HyperedgeID& num_hyperedges,
                          HypernodeID& num_hypernodes,
                          HyperedgeID& num_removed_single_pin_hyperedges,
                          HyperedgeVector& hyperedges,
                          parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                          parallel::scalable_vector<HypernodeWeight>& hypernodes_weight,
                          const bool remove_single_pin_hes) {
    ASSERT(!filename.empty(), "No filename for hypergraph file specified");
    int fd = open_file(filename);
    const size_t length = file_size(fd);
    char* mapped_file = mmap_file(fd, length);
    size_t pos = 0;

    // Read Hypergraph Header
    mt_kahypar::Type type = mt_kahypar::Type::Unweighted;
    readHGRHeader(mapped_file, pos, length, num_hyperedges, num_hypernodes, type);

    // Read Hyperedges
    num_removed_single_pin_hyperedges =
            readHyperedges(mapped_file, pos, length, num_hyperedges,
              type, hyperedges, hyperedges_weight, remove_single_pin_hes);
    num_hyperedges -= num_removed_single_pin_hyperedges;

    // Read Hypernode Weights
    readHypernodeWeights(mapped_file, pos, length, num_hypernodes, type, hypernodes_weight);
    ASSERT(pos == length);

    munmap_file(mapped_file, fd, length);
    close(fd);
  }

  Hypergraph readHypergraphFile(const std::string& filename,
                                const bool stable_construction_of_incident_edges,
                                const bool remove_single_pin_hes) {
    // Read Hypergraph File
    HyperedgeID num_hyperedges = 0;
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_removed_single_pin_hyperedges = 0;
    HyperedgeVector hyperedges;
    parallel::scalable_vector<HyperedgeWeight> hyperedges_weight;
    parallel::scalable_vector<HypernodeWeight> hypernodes_weight;
    mt_kahypar::utils::Timer::instance().start_timer("read_input", "Read Hypergraph File");
    readHypergraphFile(filename, num_hyperedges, num_hypernodes,
                       num_removed_single_pin_hyperedges, hyperedges,
                       hyperedges_weight, hypernodes_weight, remove_single_pin_hes);
    mt_kahypar::utils::Timer::instance().stop_timer("read_input");

    // Construct Hypergraph
    utils::Timer::instance().start_timer("construct_hypergraph", "Construct Hypergraph");
    Hypergraph hypergraph = HypergraphFactory::construct(
            num_hypernodes, num_hyperedges,
            hyperedges, hyperedges_weight.data(), hypernodes_weight.data(),
            stable_construction_of_incident_edges);
    hypergraph.setNumRemovedHyperedges(num_removed_single_pin_hyperedges);
    utils::Timer::instance().stop_timer("construct_hypergraph");
    return hypergraph;
  }

  void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition) {
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

  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename) {
    if (filename.empty()) {
      LOG << "No filename for partition file specified";
    } else {
      std::ofstream out_stream(filename.c_str());
      std::vector<PartitionID> partition(phg.initialNumNodes(), -1);
      for (const HypernodeID& hn : phg.nodes()) {
        ASSERT(hn < partition.size());
        partition[hn] = phg.partID(hn);
      }
      for (const PartitionID& part : partition) {
        out_stream << part << std::endl;
      }
      out_stream.close();
    }
  }

  void readMetisHeader(char* mapped_file,
                       size_t& pos,
                       const size_t length,
                       HyperedgeID& num_edges,
                       HypernodeID& num_vertices,
                       bool& has_edge_weights,
                       bool& has_vertex_weights) {
    // Skip comments
    while ( mapped_file[pos] == '%' ) {
      goto_next_line(mapped_file, pos, length);
    }

    num_vertices = read_number(mapped_file, pos, length);
    num_edges = read_number(mapped_file, pos, length);

    if ( mapped_file[pos] != '\n' ) {
      // read the (up to) three 0/1 format digits
      uint32_t format_num = read_number(mapped_file, pos, length);
      ASSERT(format_num < 100, "Vertex sizes in input file are not supported.");
      ASSERT(format_num / 10 == 0 || format_num / 10 == 1);
      has_vertex_weights = (format_num / 10 == 1);
      ASSERT(format_num % 10 == 0 || format_num % 10 == 1);
      has_edge_weights = (format_num % 10 == 1);
    }
    do_line_ending(mapped_file, pos);
  }

  struct VertexRange {
    const size_t start;
    const size_t end;
    const HypernodeID vertex_start_id;
    const HypernodeID num_vertices;
    const HyperedgeID edge_start_id;
  };

  void readVertices(char* mapped_file,
                    size_t& pos,
                    const size_t length,
                    const HyperedgeID num_edges,
                    const HypernodeID num_vertices,
                    const bool has_edge_weights,
                    const bool has_vertex_weights,
                    EdgeVector& edges,
                    parallel::scalable_vector<HyperedgeWeight>& edges_weight,
                    parallel::scalable_vector<HypernodeWeight>& vertices_weight) {
    parallel::scalable_vector<VertexRange> vertex_ranges;
    tbb::parallel_invoke([&] {
      // Sequential pass over all vertices to determine ranges in the
      // input file that are read in parallel.
      // Additionally, we need to sum the vertex degrees to determine edge indices.
      size_t current_range_start = pos;
      HypernodeID current_range_vertex_id = 0;
      HypernodeID current_range_num_vertices = 0;
      HyperedgeID current_range_edge_id = 0;
      HyperedgeID current_range_num_edges = 0;
      const HypernodeID num_vertices_per_range = std::max(
              (num_vertices / ( 2 * std::thread::hardware_concurrency())), ID(1));
      while ( current_range_vertex_id + current_range_num_vertices < num_vertices ) {
        // Skip Comments
        ASSERT(pos < length);
        while ( mapped_file[pos] == '%' ) {
          goto_next_line(mapped_file, pos, length);
          ASSERT(pos < length);
        }

        ASSERT(mapped_file[pos - 1] == '\n');
        ++current_range_num_vertices;

        // Count the forward edges, ignore backward edges.
        // This is necessary because we can only calculate unique edge ids
        // efficiently if the edges are deduplicated.
        if ( has_vertex_weights ) {
          read_number(mapped_file, pos, length);
        }
        HyperedgeID vertex_degree = 0;
        while (mapped_file[pos] != '\n' && pos < length) {
          const HypernodeID source = current_range_vertex_id + current_range_num_vertices;
          const HypernodeID target = read_number(mapped_file, pos, length);
          ASSERT(source != target);
          if ( source < target ) {
            ++vertex_degree;
          }
          if ( has_edge_weights ) {
            read_number(mapped_file, pos, length);
          }
        }
        do_line_ending(mapped_file, pos);
        current_range_num_edges += vertex_degree;

        // If there are enough vertices in the current scanned range
        // we store that range, which will be processed in parallel later
        if ( current_range_num_vertices == num_vertices_per_range ) {
          vertex_ranges.push_back(VertexRange {
                  current_range_start, pos, current_range_vertex_id, current_range_num_vertices, current_range_edge_id});
          current_range_start = pos;
          current_range_vertex_id += current_range_num_vertices;
          current_range_num_vertices = 0;
          current_range_edge_id += current_range_num_edges;
          current_range_num_edges = 0;
        }
      }
      if ( current_range_num_vertices > 0 ) {
        vertex_ranges.push_back(VertexRange {
                current_range_start, pos, current_range_vertex_id, current_range_num_vertices, current_range_edge_id});
        current_range_vertex_id += current_range_num_vertices;
        current_range_edge_id += current_range_num_edges;
      }
      ASSERT(current_range_vertex_id == num_vertices);
      ASSERT(current_range_edge_id == num_edges);
    }, [&] {
      edges.resize(num_edges);
    }, [&] {
      if ( has_edge_weights ) {
        edges_weight.resize(num_edges);
      }
    }, [&] {
      if ( has_vertex_weights ) {
        vertices_weight.resize(num_vertices);
      }
    });

    ASSERT([&]() {
        HyperedgeID last_end = 0;
        for(const auto& range: vertex_ranges) {
          if (last_end > range.start) {
            return false;
          }
          last_end = range.end;
        }
        return true;
      }()
    );

    // Process all ranges in parallel, build edge vector and assign weights
    tbb::parallel_for(0UL, vertex_ranges.size(), [&](const size_t i) {
      const VertexRange& range = vertex_ranges[i];
      size_t current_pos = range.start;
      const size_t current_end = range.end;
      HypernodeID current_vertex_id = range.vertex_start_id;
      const HypernodeID last_vertex_id = current_vertex_id + range.num_vertices;
      HyperedgeID current_edge_id = range.edge_start_id;

      while ( current_vertex_id < last_vertex_id ) {
        // Skip Comments
        ASSERT(current_pos < current_end);
        while ( mapped_file[pos] == '%' ) {
          goto_next_line(mapped_file, current_pos, current_end);
          ASSERT(current_pos < current_end);
        }

        if ( has_vertex_weights ) {
          ASSERT(current_vertex_id < vertices_weight.size());
          vertices_weight[current_vertex_id] = read_number(mapped_file, current_pos, current_end);
        }

        while ( !is_line_ending(mapped_file, current_pos) ) {
          const HypernodeID target = read_number(mapped_file, current_pos, current_end);
          ASSERT(target > 0 && (target - 1) < num_vertices, V(target));

          // process forward edges, ignore backward edges
          if ( current_vertex_id < (target - 1) ) {
            ASSERT(current_edge_id < edges.size());
            // At this point, some magic is involved:
            // In case of the graph partitioner, the right handed expression is considered a pair.
            // In case of the hypergraph partitioner, the right handed expression is considered  a vector.
            edges[current_edge_id] = {current_vertex_id, target - 1};

            if ( has_edge_weights ) {
              edges_weight[current_edge_id] = read_number(mapped_file, current_pos, current_end);
            }
            ++current_edge_id;
          } else if ( has_edge_weights ) {
            read_number(mapped_file, current_pos, current_end);
          }
        }
        do_line_ending(mapped_file, current_pos);
        ++current_vertex_id;
      }
    });
  }

  void readMetisFile(const std::string& filename,
                     HyperedgeID& num_edges,
                     HypernodeID& num_vertices,
                     EdgeVector& edges,
                     parallel::scalable_vector<HyperedgeWeight>& edges_weight,
                     parallel::scalable_vector<HypernodeWeight>& vertices_weight) {
    ASSERT(!filename.empty(), "No filename for metis file specified");
    int fd = open_file(filename);
    const size_t length = file_size(fd);
    char* mapped_file = mmap_file(fd, length);
    size_t pos = 0;

    // Read Metis Header
    bool has_edge_weights = false;
    bool has_vertex_weights = false;
    readMetisHeader(mapped_file, pos, length, num_edges,
      num_vertices, has_edge_weights, has_vertex_weights);

    // Read Vertices
    readVertices(mapped_file, pos, length, num_edges, num_vertices,
      has_edge_weights, has_vertex_weights, edges, edges_weight, vertices_weight);
    ASSERT(pos == length);

    munmap_file(mapped_file, fd, length);
    close(fd);
  }

  Hypergraph readMetisFile(const std::string& filename,
                           const bool stable_construction_of_incident_edges) {
    // Read Metis File
    HyperedgeID num_edges = 0;
    HypernodeID num_vertices = 0;
    EdgeVector edges;
    parallel::scalable_vector<HyperedgeWeight> edges_weight;
    parallel::scalable_vector<HypernodeWeight> nodes_weight;
    mt_kahypar::utils::Timer::instance().start_timer("read_input", "Read Hypergraph File");
    readMetisFile(filename, num_edges, num_vertices, edges, edges_weight, nodes_weight);
    mt_kahypar::utils::Timer::instance().stop_timer("read_input");

    // Construct Graph
    utils::Timer::instance().start_timer("construct_hypergraph", "Construct Hypergraph");
    #ifdef USE_GRAPH_PARTITIONER
    Hypergraph graph = HypergraphFactory::construct_from_graph_edges(
            num_vertices, num_edges, edges,
            edges_weight.data(), nodes_weight.data(),
            stable_construction_of_incident_edges);
    #else
    Hypergraph graph = HypergraphFactory::construct(
            num_vertices, num_edges,
            edges, edges_weight.data(), nodes_weight.data(),
            stable_construction_of_incident_edges);
    #endif
    utils::Timer::instance().stop_timer("construct_hypergraph");
    return graph;
  }

  Hypergraph readInputFile(const std::string& filename,
                           const FileFormat format,
                           const bool stable_construction_of_incident_edges,
                           const bool remove_single_pin_hes) {
    switch (format) {
      case FileFormat::hMetis: return readHypergraphFile(filename, stable_construction_of_incident_edges, remove_single_pin_hes);
      case FileFormat::Metis: return readMetisFile(filename, stable_construction_of_incident_edges);
        // omit default case to trigger compiler warning for missing cases
    }
    return Hypergraph();
  }

} // namespace