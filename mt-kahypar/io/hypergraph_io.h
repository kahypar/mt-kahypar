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

#include <string>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace io {

  using Hyperedge = parallel::scalable_vector<HypernodeID>;
  using HyperedgeVector = parallel::scalable_vector<Hyperedge>;
  #ifdef USE_GRAPH_PARTITIONER
  using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
  #else
  using EdgeVector = HyperedgeVector;
  #endif

  void readHypergraphFile(const std::string& filename,
                          HyperedgeID& num_hyperedges,
                          HypernodeID& num_hypernodes,
                          HyperedgeID& num_removed_single_pin_hyperedges,
                          HyperedgeVector& hyperedges,
                          parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                          parallel::scalable_vector<HypernodeWeight>& hypernodes_weight);

  Hypergraph readHypergraphFile(const std::string& filename,
                                const bool stable_construction_of_incident_edges = false);

  Hypergraph readMetisFile(const std::string& filename,
                           const bool stable_construction_of_incident_edges = false);

  Hypergraph readInputFile(const std::string& filename,
                           const FileFormat format,
                           const bool stable_construction_of_incident_edges = false);

  void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition);
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename);

  namespace impl {
    int open_file(const std::string& filename);

    size_t file_size(int fd);

    char* mmap_file(int fd, const size_t length);

    void munmap_file(char* mapped_file, int fd, const size_t length);

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
    void line_ending(char* mapped_file, size_t& pos) {
      unused(mapped_file);
      ASSERT(mapped_file[pos] == '\n');
      ++pos;
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    int64_t read_number(char* mapped_file, size_t& pos, const size_t length) {
      int64_t number = 0;
      for ( ; pos < length; ++pos ) {
        if ( mapped_file[pos] == ' ' || mapped_file[pos] == '\n' ) {
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

    void readMetisHeader(char* mapped_file,
                         size_t& pos,
                         const size_t length,
                         HyperedgeID& num_edges,
                         HypernodeID& num_vertices,
                         bool& has_edge_weights,
                         bool& has_vertex_weights);

    struct VertexRange {
      const size_t start;
      const size_t end;
      const HypernodeID vertex_start_id;
      const HypernodeID num_vertices;
      const HyperedgeID edge_start_id;
    };
  } // namespace

  template<typename EdgeT>
  void readVertices(char* mapped_file,
                    size_t& pos,
                    const size_t length,
                    const HyperedgeID num_edges,
                    const HypernodeID num_vertices,
                    const bool has_edge_weights,
                    const bool has_vertex_weights,
                    parallel::scalable_vector<EdgeT>& edges,
                    parallel::scalable_vector<HyperedgeWeight>& edges_weight,
                    parallel::scalable_vector<HypernodeWeight>& vertices_weight) {
    using namespace impl;

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
        line_ending(mapped_file, pos);
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

        while ( mapped_file[current_pos] != '\n' ) {
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
        line_ending(mapped_file, current_pos);
        ++current_vertex_id;
      }
    });
  }

  // EdgeT can be either parallel::scalable_vector<HypernodeID> or std::pair<HypernodeID>
  template<typename EdgeT>
  void readMetisFile(const std::string& filename,
                     HyperedgeID& num_edges,
                     HypernodeID& num_vertices,
                     parallel::scalable_vector<EdgeT>& edges,
                     parallel::scalable_vector<HyperedgeWeight>& edges_weight,
                     parallel::scalable_vector<HypernodeWeight>& vertices_weight) {
    using namespace impl;

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

}  // namespace io
}  // namespace mt_kahypar
