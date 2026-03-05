/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "hypergraph_io.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <thread>
#include <memory>
#include <vector>

#include <tbb/parallel_for.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/io_utils.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::io {

  void readFormatInfo(char* mapped_file,
                      size_t& pos,
                      const size_t length,
                      bool& has_edge_weights,
                      bool& has_vertex_weights) {
    if (!is_line_ending(mapped_file, pos)) {
      // read the (up to) three 0/1 format digits
      uint32_t format_num = read_number(mapped_file, pos, length);
      ASSERT(format_num < 100, "Vertex sizes in input file are not supported.");
      ASSERT(format_num / 10 == 0 || format_num / 10 == 1);
      has_vertex_weights = (format_num / 10 == 1);
      ASSERT(format_num % 10 == 0 || format_num % 10 == 1);
      has_edge_weights = (format_num % 10 == 1);
    }
  }

  void readHGRHeader(char* mapped_file,
                     size_t& pos,
                     const size_t length,
                     HyperedgeID& num_hyperedges,
                     HypernodeID& num_hypernodes,
                     bool& has_hyperedge_weights,
                     bool& has_vertex_weights) {
    // Skip comments
    while ( mapped_file[pos] == '%' ) {
      goto_next_line(mapped_file, pos, length);
    }

    num_hyperedges = read_number(mapped_file, pos, length);
    num_hypernodes = read_number(mapped_file, pos, length);
    readFormatInfo(mapped_file, pos, length, has_hyperedge_weights, has_vertex_weights);
    do_line_ending(mapped_file, pos);
  }

  struct HyperedgeRange {
    const size_t start;
    const size_t end;
    const HyperedgeID start_id;
    const HyperedgeID num_hyperedges;
  };

  struct HyperedgeReadResult {
    HyperedgeReadResult() :
      num_removed_single_pin_hyperedges(0),
      num_duplicated_pins(0),
      num_hes_with_duplicated_pins(0) { }

    size_t num_removed_single_pin_hyperedges;
    size_t num_duplicated_pins;
    size_t num_hes_with_duplicated_pins;
  };

  HyperedgeReadResult readHyperedges(char* mapped_file,
                                     size_t& pos,
                                     const size_t length,
                                     const HyperedgeID num_hyperedges,
                                     const bool has_hyperedge_weights,
                                     HyperedgeVector& hyperedges,
                                     vec<HyperedgeWeight>& hyperedges_weight,
                                     const bool remove_single_pin_hes) {
    HyperedgeReadResult res;
    vec<LineRange> line_ranges;
    vec<HyperedgeVector> local_edges;
    vec<vec<HyperedgeWeight>> local_edges_weight;
    tbb::parallel_invoke([&] {
      // Sequential pass over all lines to determine ranges in the
      // input file that are read in parallel.
      line_ranges = split_lines(mapped_file, pos, length, num_hyperedges);
    }, [&] {
      hyperedges.resize(num_hyperedges);
    }, [&] {
      local_edges.resize(num_parallel_ranges());
      local_edges_weight.resize(num_parallel_ranges());
    }, [&] {
      if ( has_hyperedge_weights ) {
        hyperedges_weight.resize(num_hyperedges);
      }
    });

    // Process all ranges in parallel and build hyperedge vector
    tbb::parallel_for(UL(0), line_ranges.size(), [&](const size_t i) {
      HyperedgeVector& my_edges = local_edges[i];
      vec<HyperedgeWeight>& my_edges_weight = local_edges_weight[i];

      const LineRange& range = line_ranges[i];
      size_t current_pos = range.start;
      const size_t current_end = range.end;
      HyperedgeID current_hyperedge_id = range.line_start_index;
      const HyperedgeID last_hyperedge_id = current_hyperedge_id + range.num_lines;

      while ( current_hyperedge_id < last_hyperedge_id ) {
        // Skip Comments
        ASSERT(current_pos < current_end);
        while ( mapped_file[current_pos] == '%' ) {
          goto_next_line(mapped_file, current_pos, current_end);
          ASSERT(current_pos < current_end);
        }

        if ( has_hyperedge_weights ) {
          my_edges_weight.push_back(read_number(mapped_file, current_pos, current_end));
        }

        Hyperedge hyperedge;
        // Note, a hyperedge line must contain at least one pin
        HypernodeID pin = read_number(mapped_file, current_pos, current_end);
        ASSERT(pin > 0, V(current_hyperedge_id));
        hyperedge.push_back(pin - 1);
        while ( !is_line_ending(mapped_file, current_pos) ) {
          pin = read_number(mapped_file, current_pos, current_end);
          ASSERT(pin > 0, V(current_hyperedge_id));
          hyperedge.push_back(pin - 1);
        }
        do_line_ending(mapped_file, current_pos);
        ++current_hyperedge_id;

        if ( remove_single_pin_hes ) {
          // Detect duplicated pins
          std::sort(hyperedge.begin(), hyperedge.end());
          size_t j = 1;
          for ( size_t i = 1; i < hyperedge.size(); ++i ) {
            if ( hyperedge[j - 1] != hyperedge[i] ) {
              std::swap(hyperedge[i], hyperedge[j++]);
            }
          }
          if ( j < hyperedge.size() ) {
            // Remove duplicated pins
            __atomic_fetch_add(&res.num_hes_with_duplicated_pins, 1, __ATOMIC_RELAXED);
            __atomic_fetch_add(&res.num_duplicated_pins, hyperedge.size() - j, __ATOMIC_RELAXED);
            hyperedge.resize(j);
          }
        }

        if ( !remove_single_pin_hes || hyperedge.size() > 1 ) {
          my_edges.emplace_back(std::move(hyperedge));
        } else {
          __atomic_fetch_add(&res.num_removed_single_pin_hyperedges, 1, __ATOMIC_RELAXED);
          my_edges_weight.pop_back();  // also remove the previously added weight
        };
      }
    });

    // update lengths for removed hyperedges
    const HyperedgeID tmp_num_hyperedges = num_hyperedges - res.num_removed_single_pin_hyperedges;
    hyperedges.resize(tmp_num_hyperedges);
    if ( has_hyperedge_weights ) {
      hyperedges_weight.resize(tmp_num_hyperedges);
    }

    // finally, we copy the local edge lists to the global list
    copy_to_global_list(local_edges, hyperedges);
    if ( has_hyperedge_weights ) {
      copy_to_global_list(local_edges_weight, hyperedges_weight);
      ASSERT(hyperedges.size() == hyperedges_weight.size());
    }

    // parallel free
    tbb::parallel_invoke([&] {
      parallel::free(line_ranges);
      if ( has_hyperedge_weights ) {
        parallel::parallel_free(local_edges_weight);
      }
    }, [&] {
      parallel::parallel_free(local_edges);
    });

    return res;
  }

  void readHypernodeWeights(char* mapped_file,
                            size_t& pos,
                            const size_t length,
                            const HypernodeID num_hypernodes,
                            const bool has_hypernode_weights,
                            vec<HypernodeWeight>& hypernodes_weight) {
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
                          vec<HyperedgeWeight>& hyperedges_weight,
                          vec<HypernodeWeight>& hypernodes_weight,
                          const bool remove_single_pin_hes,
                          const bool print_warnings) {
    ASSERT(!filename.empty(), "No filename for hypergraph file specified");
    FileHandle handle = mmap_file(filename);
    size_t pos = 0;

    // Read Hypergraph Header
    bool has_hyperedge_weights = false;
    bool has_vertex_weights = false;
    readHGRHeader(handle.mapped_file, pos, handle.length, num_hyperedges,
      num_hypernodes, has_hyperedge_weights, has_vertex_weights);

    // Read Hyperedges
    HyperedgeReadResult res =
            readHyperedges(handle.mapped_file, pos, handle.length, num_hyperedges,
              has_hyperedge_weights, hyperedges, hyperedges_weight, remove_single_pin_hes);
    num_hyperedges -= res.num_removed_single_pin_hyperedges;
    num_removed_single_pin_hyperedges = res.num_removed_single_pin_hyperedges;

    if ( print_warnings && res.num_hes_with_duplicated_pins > 0 ) {
      WARNING("Removed " << res.num_duplicated_pins << " duplicated pins in "
        << res.num_hes_with_duplicated_pins << " hyperedges!");
    }

    // Read Hypernode Weights
    readHypernodeWeights(handle.mapped_file, pos, handle.length, num_hypernodes,
      has_vertex_weights, hypernodes_weight);

    // Check the end of the file
    while ( handle.mapped_file[pos] == '%' ) {
      goto_next_line(handle.mapped_file, pos, handle.length);
    }
    ASSERT(pos == handle.length);

    munmap_file(handle);
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
    readFormatInfo(mapped_file, pos, length, has_edge_weights, has_vertex_weights);
    do_line_ending(mapped_file, pos);
  }

  struct VertexRange {
    const size_t start;
    const size_t end;
    const HypernodeID vertex_start_id;
    const HypernodeID num_vertices;
    const HyperedgeID edge_start_id;
  };

  template<typename EdgeT>
  void readVertices(char* mapped_file,
                    size_t& pos,
                    const size_t length,
                    const HyperedgeID num_edges,
                    const HypernodeID num_vertices,
                    const bool has_edge_weights,
                    const bool has_vertex_weights,
                    vec<EdgeT>& edges,
                    vec<HyperedgeWeight>& edges_weight,
                    vec<HypernodeWeight>& vertices_weight) {
    vec<LineRange> line_ranges;
    vec<vec<EdgeT>> local_edges;
    vec<vec<HyperedgeWeight>> local_edges_weight;
    tbb::parallel_invoke([&] {
      // Sequential pass over all lines to determine ranges in the
      // input file that are read in parallel.
      line_ranges = split_lines(mapped_file, pos, length, num_vertices);
    }, [&] {
      edges.resize(num_edges);
    }, [&] {
      local_edges.resize(num_parallel_ranges());
      local_edges_weight.resize(num_parallel_ranges());
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
    tbb::parallel_for(UL(0), line_ranges.size(), [&](const size_t i) {
      vec<EdgeT>& my_edges = local_edges[i];
      vec<HyperedgeWeight>& my_edges_weight = local_edges_weight[i];

      const LineRange& range = line_ranges[i];
      size_t current_pos = range.start;
      const size_t current_end = range.end;
      HypernodeID current_vertex_id = range.line_start_index;
      const HypernodeID last_vertex_id = current_vertex_id + range.num_lines;

      while ( current_vertex_id < last_vertex_id ) {
        // Skip Comments
        ASSERT(current_pos < current_end);
        while ( mapped_file[current_pos] == '%' ) {
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
            // note: the edge type can be either std::pair<HypernodeID, HypernodeID> or vec<HypernodeID>,
            // the initializer list {u, v} is a valid constructor call for both
            my_edges.push_back({current_vertex_id, target - 1});

            if ( has_edge_weights ) {
              my_edges_weight.push_back(read_number(mapped_file, current_pos, current_end));
            }
          } else if ( has_edge_weights ) {
            read_number(mapped_file, current_pos, current_end);
          }
        }
        do_line_ending(mapped_file, current_pos);
        ++current_vertex_id;
      }
    });

    // finally, we copy the local edge lists to the global list
    copy_to_global_list(local_edges, edges);
    if ( has_edge_weights ) {
      copy_to_global_list(local_edges_weight, edges_weight);
      ASSERT(edges.size() == edges_weight.size());
    }

    // parallel free
    tbb::parallel_invoke([&] {
      parallel::free(line_ranges);
      if ( has_edge_weights ) {
        parallel::parallel_free(local_edges_weight);
      }
    }, [&] {
      parallel::parallel_free(local_edges);
    });
  }

  template<typename EdgeT>
  void readGraphFile(const std::string& filename,
                     HyperedgeID& num_edges,
                     HypernodeID& num_vertices,
                     vec<EdgeT>& edges,
                     vec<HyperedgeWeight>& edges_weight,
                     vec<HypernodeWeight>& vertices_weight) {
    ASSERT(!filename.empty(), "No filename for metis file specified");
    FileHandle handle = mmap_file(filename);
    size_t pos = 0;

    // Read Metis Header
    bool has_edge_weights = false;
    bool has_vertex_weights = false;
    readMetisHeader(handle.mapped_file, pos, handle.length, num_edges,
      num_vertices, has_edge_weights, has_vertex_weights);

    // Read Vertices
    readVertices(handle.mapped_file, pos, handle.length, num_edges, num_vertices,
      has_edge_weights, has_vertex_weights, edges, edges_weight, vertices_weight);

    // Check the end of the file
    while ( handle.mapped_file[pos] == '%' ) {
      goto_next_line(handle.mapped_file, pos, handle.length);
    }
    ASSERT(pos == handle.length);

    munmap_file(handle);
  }

  template<typename InitFunc>
  void readPartitionFileImpl(const std::string& filename, HypernodeID num_nodes, InitFunc init_func) {
    ASSERT(!filename.empty(), "No filename for partition file specified");
    std::ifstream file(filename);
    if (file) {
      PartitionID* partition = init_func();
      int part;
      HypernodeID hn = 0;
      while (file >> part) {
        if (hn >= num_nodes) {
          throw InvalidInputException(std::string("Input file has more entries than the number of nodes: ") + filename);
        }
        partition[hn++] = part;
      }
      file.close();
      if (hn < num_nodes) {
        throw InvalidInputException(std::string("Input file has less entries than the number of nodes: ") + filename);
      }
    } else {
      throw InvalidInputException(std::string("File not found: ") + filename);
    }
  }

  void readPartitionFile(const std::string& filename, HypernodeID num_nodes, std::vector<PartitionID>& partition) {
    readPartitionFileImpl(filename, num_nodes, [&]{
      partition.clear();
      partition.resize(num_nodes);
      return partition.data();
    });
  }

  void readPartitionFile(const std::string& filename, HypernodeID num_nodes, PartitionID* partition) {
    readPartitionFileImpl(filename, num_nodes, [=]{ return partition; });
  }

  template<typename PartitionedHypergraph>
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename) {
    if (filename.empty()) {
      throw InvalidInputException("No filename for output partition file specified");
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

  namespace {
  #define WRITE_PARTITION_FILE(X) void writePartitionFile(const X& phg, const std::string& filename)
  }

  INSTANTIATE_FUNC_WITH_PARTITIONED_HG(WRITE_PARTITION_FILE)

  template void readGraphFile(const std::string& filename,
                              HyperedgeID& num_hyperedges,
                              HypernodeID& num_hypernodes,
                              HyperedgeVector& hyperedges,
                              vec<HyperedgeWeight>& hyperedges_weight,
                              vec<HypernodeWeight>& hypernodes_weight);

  template void readGraphFile(const std::string& filename,
                              HyperedgeID& num_hyperedges,
                              HypernodeID& num_hypernodes,
                              EdgeVector& hyperedges,
                              vec<HyperedgeWeight>& hyperedges_weight,
                              vec<HypernodeWeight>& hypernodes_weight);

} // namespace
