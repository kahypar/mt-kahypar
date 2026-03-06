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
#include <chrono>
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
#include "mt-kahypar/utils/deduplicate.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::io {

  using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

  void readFormatInfo(char* mapped_file,
                      size_t& pos,
                      const size_t current_line,
                      const size_t length,
                      bool& has_edge_weights,
                      bool& has_vertex_weights) {
    if (!is_line_ending(mapped_file, pos)) {
      // read the (up to) three 0/1 format digits
      uint32_t format_num = read_number(mapped_file, pos, current_line, length);
      ASSERT(format_num < 100, "Vertex sizes in input file are not supported.");
      ASSERT(format_num / 10 == 0 || format_num / 10 == 1);
      has_vertex_weights = (format_num / 10 == 1);
      ASSERT(format_num % 10 == 0 || format_num % 10 == 1);
      has_edge_weights = (format_num % 10 == 1);
    }
  }

  void readHGRHeader(char* mapped_file,
                     size_t& pos,
                     size_t& current_line,
                     const size_t length,
                     HyperedgeID& num_hyperedges,
                     HypernodeID& num_hypernodes,
                     bool& has_hyperedge_weights,
                     bool& has_vertex_weights) {
    // Skip comments
    while ( mapped_file[pos] == '%' ) {
      goto_next_line(mapped_file, pos, current_line, length);
    }

    num_hyperedges = read_number(mapped_file, pos, current_line, length);
    num_hypernodes = read_number(mapped_file, pos, current_line, length);
    readFormatInfo(mapped_file, pos, current_line, length, has_hyperedge_weights, has_vertex_weights);
    do_line_ending(mapped_file, pos, current_line);
  }

  struct HyperedgeRange {
    const size_t start;
    const size_t end;
    const HyperedgeID start_id;
    const HyperedgeID num_hyperedges;
  };

  bool isSinglePinHyperedge(char* mapped_file,
                            size_t pos,
                            const size_t length,
                            const bool has_hyperedge_weights) {
    ASSERT(!is_line_ending(mapped_file, pos), "Empty line where hyperedge was expected");
    ASSERT(mapped_file[pos] != ' ', "Line may not start with space");
    ++pos;

    size_t num_spaces = 0;
    for ( ; pos < length; ++pos ) {
      if (is_line_ending(mapped_file, pos)) {
        break;
        // shift the check by -1 so that we can tolerate a space at the end of the line
      } else if ( mapped_file[pos - 1] == ' ' ) {
        ++num_spaces;
      }

      if ( num_spaces == 2 ) {
        break;
      }
    }
    return has_hyperedge_weights ? num_spaces == 1 : num_spaces == 0;
  }

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
                                     size_t& current_line,
                                     const size_t length,
                                     const HyperedgeID num_hyperedges,
                                     const bool has_hyperedge_weights,
                                     HyperedgeVector& hyperedges,
                                     vec<HyperedgeWeight>& hyperedges_weight,
                                     const bool remove_single_pin_hes) {
    HyperedgeReadResult res;
    vec<HyperedgeRange> hyperedge_ranges;
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
          goto_next_line(mapped_file, pos, current_line, length);
          ASSERT(pos < length);
        }

        // This check is fine even with windows line endings!
        ASSERT(mapped_file[pos - 1] == '\n');
        if ( !remove_single_pin_hes || !isSinglePinHyperedge(mapped_file, pos, length, has_hyperedge_weights) ) {
          ++current_range_num_hyperedges;
        } else {
          ++res.num_removed_single_pin_hyperedges;
        }
        ++current_num_hyperedges;
        goto_next_line(mapped_file, pos, current_line, length);

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

    const HyperedgeID tmp_num_hyperedges = num_hyperedges - res.num_removed_single_pin_hyperedges;
    hyperedges.resize(tmp_num_hyperedges);
    if ( has_hyperedge_weights ) {
      hyperedges_weight.resize(tmp_num_hyperedges);
    }

    // Process all ranges in parallel and build hyperedge vector
    tbb::parallel_for(UL(0), hyperedge_ranges.size(), [&](const size_t i) {
      HyperedgeRange& range = hyperedge_ranges[i];
      size_t current_pos = range.start;
      size_t current_line = 0;
      const size_t current_end = range.end;
      HyperedgeID current_id = range.start_id;
      const HyperedgeID last_id = current_id + range.num_hyperedges;

      while ( current_id < last_id ) {
        // Skip Comments
        ASSERT(current_pos < current_end);
        while ( mapped_file[current_pos] == '%' ) {
          goto_next_line(mapped_file, current_pos, current_line, current_end);
          ASSERT(current_pos < current_end);
        }

        if ( !remove_single_pin_hes || !isSinglePinHyperedge(mapped_file, current_pos, current_end, has_hyperedge_weights) ) {
          ASSERT(current_id < hyperedges.size());
          if ( has_hyperedge_weights ) {
            hyperedges_weight[current_id] = read_number(mapped_file, current_pos, current_line, current_end);
          }

          Hyperedge& hyperedge = hyperedges[current_id];
          // Note, a hyperedge line must contain at least one pin
          HypernodeID pin = read_number(mapped_file, current_pos, current_line, current_end);
          ASSERT(pin > 0, V(current_id));
          hyperedge.push_back(pin - 1);
          while ( !is_line_ending(mapped_file, current_pos) ) {
            pin = read_number(mapped_file, current_pos, current_line, current_end);
            ASSERT(pin > 0, V(current_id));
            hyperedge.push_back(pin - 1);
          }
          do_line_ending(mapped_file, current_pos, current_line);

          utils::deduplicateHyperedgePins(hyperedge, res.num_duplicated_pins, res.num_hes_with_duplicated_pins);
          ASSERT(hyperedge.size() >= 2);
          ++current_id;
        } else {
          goto_next_line(mapped_file, current_pos, current_line, current_end);
        }
      }
    });
    return res;
  }

  void readHypernodeWeights(char* mapped_file,
                            size_t& pos,
                            size_t& current_line,
                            const size_t length,
                            const HypernodeID num_hypernodes,
                            const bool has_hypernode_weights,
                            vec<HypernodeWeight>& hypernodes_weight) {
    if ( has_hypernode_weights ) {
      hypernodes_weight.resize(num_hypernodes);
      for ( HypernodeID hn = 0; hn < num_hypernodes; ++hn ) {
        ASSERT(pos > 0 && pos < length);
        ASSERT(mapped_file[pos - 1] == '\n');
        hypernodes_weight[hn] = read_number(mapped_file, pos, current_line, length);
        do_line_ending(mapped_file, pos, current_line);
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
    size_t current_line = 0;

    // Read Hypergraph Header
    bool has_hyperedge_weights = false;
    bool has_vertex_weights = false;
    readHGRHeader(handle.mapped_file, pos, current_line, handle.length,
      num_hyperedges, num_hypernodes, has_hyperedge_weights, has_vertex_weights);

    // Read Hyperedges
    HyperedgeReadResult res =
            readHyperedges(handle.mapped_file, pos, current_line, handle.length, num_hyperedges,
              has_hyperedge_weights, hyperedges, hyperedges_weight, remove_single_pin_hes);
    num_hyperedges -= res.num_removed_single_pin_hyperedges;
    num_removed_single_pin_hyperedges = res.num_removed_single_pin_hyperedges;

    if ( print_warnings && res.num_hes_with_duplicated_pins > 0 ) {
      WARNING("Removed" << res.num_duplicated_pins << "duplicated pins in"
        << res.num_hes_with_duplicated_pins << "hyperedges!");
    }

    // Read Hypernode Weights
    readHypernodeWeights(handle.mapped_file, pos, current_line, handle.length, num_hypernodes,
      has_vertex_weights, hypernodes_weight);

    // Check the end of the file
    while ( handle.mapped_file[pos] == '%' ) {
      goto_next_line(handle.mapped_file, pos, current_line, handle.length);
    }
    if (pos != handle.length) {
      parsing_exception(current_line, "unexpected content; hypergraph data is already complete");
    }

    munmap_file(handle);
  }

  void readMetisHeader(char* mapped_file,
                       size_t& pos,
                       size_t& current_line,
                       const size_t length,
                       HyperedgeID& num_edges,
                       HypernodeID& num_vertices,
                       bool& has_edge_weights,
                       bool& has_vertex_weights) {
    // Skip comments
    while ( mapped_file[pos] == '%' ) {
      goto_next_line(mapped_file, pos, current_line, length);
    }

    num_vertices = read_number(mapped_file, pos, current_line, length);
    num_edges = read_number(mapped_file, pos, current_line, length);
    readFormatInfo(mapped_file, pos, current_line, length, has_edge_weights, has_vertex_weights);
    do_line_ending(mapped_file, pos, current_line);
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
                    size_t& current_line,
                    const size_t length,
                    const HyperedgeID num_edges,
                    const HypernodeID num_vertices,
                    const bool has_edge_weights,
                    const bool has_vertex_weights,
                    vec<EdgeT>& edges,
                    vec<HyperedgeWeight>& edges_weight,
                    vec<HypernodeWeight>& vertices_weight) {
    HighResClockTimepoint first_pass_start = std::chrono::high_resolution_clock::now();
    vec<LineRange> line_ranges;
    vec<vec<EdgeT>> local_edges;
    vec<vec<HyperedgeWeight>> local_edges_weight;
    tbb::parallel_invoke([&] {
      HighResClockTimepoint split_start = std::chrono::high_resolution_clock::now();

      // Sequential pass over all lines to determine ranges in the
      // input file that are read in parallel.
      line_ranges = split_lines(mapped_file, pos, current_line, length, num_vertices);

      HighResClockTimepoint split_end = std::chrono::high_resolution_clock::now();
      double split_time = std::chrono::duration<double>(split_end - split_start).count();
      LOG << V(split_time);
    }, [&] {
      HighResClockTimepoint edge_alloc_start = std::chrono::high_resolution_clock::now();
      edges.resize(num_edges);
      HighResClockTimepoint edge_alloc_end = std::chrono::high_resolution_clock::now();
      double edge_alloc_time = std::chrono::duration<double>(edge_alloc_end - edge_alloc_start).count();
      LOG << V(edge_alloc_time);
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

    HighResClockTimepoint second_pass_start = std::chrono::high_resolution_clock::now();
    double first_pass_time = std::chrono::duration<double>(second_pass_start - first_pass_start).count();
    LOG << V(first_pass_time);

    // Process all ranges in parallel, build edge vector and assign weights
    tbb::parallel_for(UL(0), line_ranges.size(), [&](const size_t i) {
      vec<EdgeT>& my_edges = local_edges[i];
      vec<HyperedgeWeight>& my_edges_weight = local_edges_weight[i];

      const LineRange& range = line_ranges[i];
      size_t current_pos = range.start;
      size_t current_line = range.first_line_number_in_file;
      const size_t current_end = range.end;
      HypernodeID current_vertex_id = range.line_start_index;
      const HypernodeID last_vertex_id = current_vertex_id + range.num_lines;

      while ( current_vertex_id < last_vertex_id ) {
        // Skip Comments
        ASSERT(current_pos < current_end);
        while ( mapped_file[current_pos] == '%' ) {
          goto_next_line(mapped_file, current_pos, current_line, current_end);
          ASSERT(current_pos < current_end);
        }

        if ( has_vertex_weights ) {
          ASSERT(current_vertex_id < vertices_weight.size());
          vertices_weight[current_vertex_id] = read_number(mapped_file, current_pos, current_line, current_end);
        }

        while ( !is_line_ending(mapped_file, current_pos) ) {
          const HypernodeID target = read_number(mapped_file, current_pos, current_line, current_end);
          ASSERT(target > 0 && (target - 1) < num_vertices, V(target));

          // process forward edges, ignore backward edges
          if ( current_vertex_id < (target - 1) ) {
            my_edges.push_back({current_vertex_id, target - 1});

            if ( has_edge_weights ) {
              my_edges_weight.push_back(read_number(mapped_file, current_pos, current_line, current_end));
            }
          } else if ( has_edge_weights ) {
            read_number(mapped_file, current_pos, current_line, current_end);
          }
        }
        do_line_ending(mapped_file, current_pos, current_line);
        ++current_vertex_id;
      }
    });

    HighResClockTimepoint second_pass_end = std::chrono::high_resolution_clock::now();
    double second_pass_time = std::chrono::duration<double>(second_pass_end - second_pass_start).count();
    LOG << V(second_pass_time);

    // finally, we copy the local edge lists to the global list
    copy_to_global_list(local_edges, edges);
    if ( has_edge_weights ) {
      copy_to_global_list(local_edges_weight, edges_weight);
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

    HighResClockTimepoint copy_end = std::chrono::high_resolution_clock::now();
    double copy_time = std::chrono::duration<double>(copy_end - second_pass_end).count();
    LOG << V(copy_time);
    double total_vertices_time = std::chrono::duration<double>(copy_end - first_pass_start).count();
    LOG << V(total_vertices_time);
  }

  template<typename EdgeT>
  void readGraphFile(const std::string& filename,
                     HyperedgeID& num_edges,
                     HypernodeID& num_vertices,
                     vec<EdgeT>& edges,
                     vec<HyperedgeWeight>& edges_weight,
                     vec<HypernodeWeight>& vertices_weight) {
    ASSERT(!filename.empty(), "No filename for metis file specified");
    HighResClockTimepoint mmap_start = std::chrono::high_resolution_clock::now();
    FileHandle handle = mmap_file(filename);
    HighResClockTimepoint mmap_end = std::chrono::high_resolution_clock::now();
    double mmap_time = std::chrono::duration<double>(mmap_end - mmap_start).count();
    LOG << V(mmap_time);

    size_t pos = 0;
    size_t current_line = 1;

    // Read Metis Header
    bool has_edge_weights = false;
    bool has_vertex_weights = false;
    readMetisHeader(handle.mapped_file, pos, current_line, handle.length, num_edges,
      num_vertices, has_edge_weights, has_vertex_weights);

    // Read Vertices
    readVertices(handle.mapped_file, pos, current_line, handle.length, num_edges, num_vertices,
      has_edge_weights, has_vertex_weights, edges, edges_weight, vertices_weight);

    // Check the end of the file
    while ( handle.mapped_file[pos] == '%' ) {
      goto_next_line(handle.mapped_file, pos, current_line, handle.length);
      ++current_line;
    }
    if (pos != handle.length) {
      parsing_exception(current_line, "unexpected content; graph data is already complete");
    }

    HighResClockTimepoint unmap_start = std::chrono::high_resolution_clock::now();
    munmap_file(handle);
    HighResClockTimepoint unmap_end = std::chrono::high_resolution_clock::now();
    double unmap_time = std::chrono::duration<double>(unmap_end - unmap_start).count();
    LOG << V(unmap_time);
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
