/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "libmtkahypar.h"
#include "libmtkahyparhgp.h"
#include "libmtkahypargp.h"

#include "mt-kahypar/parallel/tbb_initializer.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/macros.h"

namespace {
  template<typename T>
  using vec = mt_kahypar::parallel::scalable_vector<T>;
}


mt_kahypar_context_t* mt_kahypar_context_new() {
  return reinterpret_cast<mt_kahypar_context_t*>(new mt_kahypar::Context(false));
}

void mt_kahypar_free_context(mt_kahypar_context_t* context) {
  if (context == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::Context*>(context);
}

void mt_kahypar_configure_context_from_file(mt_kahypar_context_t* kahypar_context,
                                            const char* ini_file_name) {
  mt_kahypar::parseIniToContext(
    *reinterpret_cast<mt_kahypar::Context*>(kahypar_context), ini_file_name);
}

void mt_kahypar_load_preset(mt_kahypar_context_t* context,
                                           const mt_kahypar_preset_type_t preset) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  switch(preset) {
    case DETERMINISTIC:
      c.load_deterministic_preset();
      break;
    case SPEED:
      c.load_default_preset();
      break;
    case HIGH_QUALITY:
      c.load_default_flow_preset();
      break;
  }
}

int mt_kahypar_set_context_parameter(mt_kahypar_context_t* context,
                                     const mt_kahypar_context_parameter_type_t type,
                                     const char* value) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  switch(type) {
    case NUM_BLOCKS:
      c.partition.k = atoi(value);
      if ( c.partition.k > 0 ) return 0; /** success **/
      else return 2; /** integer conversion error **/
    case EPSILON:
      c.partition.epsilon = atof(value);
      return 0;
    case OBJECTIVE:
      {
        std::string objective(value);
        if ( objective == "km1" ) {
          c.partition.objective = mt_kahypar::Objective::km1;
          return 0;
        } else if ( objective == "cut" ) {
          c.partition.objective = mt_kahypar::Objective::cut;
          return 0;
        }
        return 3;
      }
    case SEED:
      c.partition.seed = atoi(value);
      return 0;
    case NUM_VCYCLES:
      c.partition.num_vcycles = atoi(value);
      return 0;
    case VERBOSE:
      c.partition.verbose_output = atoi(value);
      return 0;
  }
  return 1; /** no valid parameter type **/
}

void mt_kahypar_set_partitioning_parameters(mt_kahypar_context_t* context,
                                            const mt_kahypar_partition_id_t num_blocks,
                                            const double epsilon,
                                            const mt_kahypar_objective_t objective,
                                            const size_t seed) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.k = num_blocks;
  c.partition.epsilon = epsilon;
  c.partition.objective = objective == KM1 ?
    mt_kahypar::Objective::km1 : mt_kahypar::Objective::cut;
  c.partition.seed = seed;
}

void mt_kahypar_set_individual_target_block_weights(mt_kahypar_context_t* context,
                                                    const mt_kahypar_partition_id_t num_blocks,
                                                    const mt_kahypar_hypernode_weight_t* block_weights) {
  mt_kahypar::Context& c = *reinterpret_cast<mt_kahypar::Context*>(context);
  c.partition.use_individual_part_weights = true;
  c.partition.max_part_weights.assign(num_blocks, 0);
  for ( mt_kahypar_partition_id_t i = 0; i < num_blocks; ++i ) {
    c.partition.max_part_weights[i] = block_weights[i];
  }
}


void mt_kahypar_initialize_thread_pool(const size_t num_threads,
                                       const bool interleaved_allocations) {
  size_t P = num_threads;
  size_t num_available_cpus = mt_kahypar::HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
      << "Setting number of threads from" << num_threads
      << "to" << num_available_cpus);
    P = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  mt_kahypar::TBBInitializer::instance(P);

  if ( interleaved_allocations ) {
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
    mt_kahypar::parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
    hwloc_bitmap_free(cpuset);
  }
}

mt_kahypar_hypergraph_t* mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                              const mt_kahypar_context_t* context,
                                                              const mt_kahypar_file_format_type_t file_format) {
  return hgp::mt_kahypar_read_hypergraph_from_file(file_name, context, file_format);
}

mt_kahypar_graph_t* mt_kahypar_read_graph_from_file(const char* file_name,
                                                    const mt_kahypar_context_t* context,
                                                    const mt_kahypar_file_format_type_t file_format) {
  return gp::mt_kahypar_read_graph_from_file(file_name, context, file_format);
}

mt_kahypar_hypergraph_t* mt_kahypar_create_hypergraph(const mt_kahypar_hypernode_id_t num_vertices,
                                                      const mt_kahypar_hyperedge_id_t num_hyperedges,
                                                      const size_t* hyperedge_indices,
                                                      const mt_kahypar_hyperedge_id_t* hyperedges,
                                                      const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                                      const mt_kahypar_hypernode_weight_t* vertex_weights) {
  return hgp::mt_kahypar_create_hypergraph(num_vertices, num_hyperedges,
    hyperedge_indices, hyperedges, hyperedge_weights, vertex_weights);
}

mt_kahypar_graph_t* mt_kahypar_create_graph(const mt_kahypar_hypernode_id_t num_vertices,
                                            const mt_kahypar_hyperedge_id_t num_edges,
                                            const mt_kahypar_hypernode_id_t* edges,
                                            const mt_kahypar_hyperedge_weight_t* edge_weights,
                                            const mt_kahypar_hypernode_weight_t* vertex_weights) {
  return gp::mt_kahypar_create_graph(num_vertices, num_edges, edges, edge_weights, vertex_weights);
}

void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t* hypergraph) {
  hgp::mt_kahypar_free_hypergraph(hypergraph);
}

void mt_kahypar_free_graph(mt_kahypar_graph_t* graph) {
  gp::mt_kahypar_free_graph(graph);
}

mt_kahypar_hypernode_id_t mt_kahypar_num_hypernodes(mt_kahypar_hypergraph_t* hypergraph) {
  return hgp::mt_kahypar_num_nodes(hypergraph);
}

mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_graph_t* graph) {
  return gp::mt_kahypar_num_nodes(graph);
}

mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t* hypergraph) {
  return hgp::mt_kahypar_num_hyperedges(hypergraph);
}

mt_kahypar_hyperedge_id_t mt_kahypar_num_edges(mt_kahypar_graph_t* graph) {
  return gp::mt_kahypar_num_edges(graph);
}

mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t* hypergraph) {
  return hgp::mt_kahypar_num_pins(hypergraph);
}

mt_kahypar_hypernode_id_t mt_kahypar_hypergraph_weight(mt_kahypar_hypergraph_t* hypergraph) {
  return hgp::mt_kahypar_total_weight(hypergraph);
}

mt_kahypar_hypernode_id_t mt_kahypar_graph_weight(mt_kahypar_graph_t* graph) {
  return gp::mt_kahypar_total_weight(graph);
}

void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  hgp::mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}

void mt_kahypar_free_partitioned_graph(mt_kahypar_partitioned_graph_t* partitioned_graph) {
  gp::mt_kahypar_free_partitioned_graph(partitioned_graph);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                     mt_kahypar_context_t* context) {
  return hgp::mt_kahypar_partition(hypergraph, context);
}

mt_kahypar_partitioned_graph_t* mt_kahypar_partition_graph(mt_kahypar_graph_t* graph,
                                                           mt_kahypar_context_t* context) {
  return gp::mt_kahypar_partition(graph, context);
}


void mt_kahypar_improve_hypergraph_partition(mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                             mt_kahypar_context_t* context,
                                             const size_t num_vcycles) {
  hgp::mt_kahypar_improve_partition(partitioned_hg, context, num_vcycles);
}

void mt_kahypar_improve_graph_partition(mt_kahypar_partitioned_graph_t* partitioned_graph,
                                        mt_kahypar_context_t* context,
                                        const size_t num_vcycles) {
  gp::mt_kahypar_improve_partition(partitioned_graph, context, num_vcycles);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                              const mt_kahypar_partition_id_t num_blocks,
                                                                              const mt_kahypar_partition_id_t* partition) {
  return hgp::mt_kahypar_create_partitioned_hypergraph(hypergraph, num_blocks, partition);
}

mt_kahypar_partitioned_graph_t* mt_kahypar_create_partitioned_graph(mt_kahypar_graph_t* graph,
                                                                    const mt_kahypar_partition_id_t num_blocks,
                                                                    const mt_kahypar_partition_id_t* partition) {
  return gp::mt_kahypar_create_partitioned_graph(graph, num_blocks, partition);
}

mt_kahypar_partitioned_hypergraph_t* mt_kahypar_read_hypergraph_partition_from_file(mt_kahypar_hypergraph_t* hypergraph,
                                                                                    const mt_kahypar_partition_id_t num_blocks,
                                                                                    const char* partition_file) {
  return hgp::mt_kahypar_read_partition_from_file(hypergraph, num_blocks, partition_file);
}

mt_kahypar_partitioned_graph_t* mt_kahypar_read_graph_partition_from_file(mt_kahypar_graph_t* graph,
                                                                          const mt_kahypar_partition_id_t num_blocks,
                                                                          const char* partition_file) {
  return gp::mt_kahypar_read_partition_from_file(graph, num_blocks, partition_file);
}

void mt_kahypar_write_hypergraph_partition_to_file(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                   const char* partition_file) {
  hgp::mt_kahypar_write_partition_to_file(partitioned_hg, partition_file);
}

void mt_kahypar_write_graph_partition_to_file(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                              const char* partition_file) {
  gp::mt_kahypar_write_partition_to_file(partitioned_graph, partition_file);
}

void mt_kahypar_get_hypergraph_partition(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                         mt_kahypar_partition_id_t* partition) {
  hgp::mt_kahypar_get_partition(partitioned_hg, partition);
}

void mt_kahypar_get_graph_partition(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                    mt_kahypar_partition_id_t* partition) {
  gp::mt_kahypar_get_partition(partitioned_graph, partition);
}

void mt_kahypar_get_hypergraph_block_weights(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                             mt_kahypar_hypernode_weight_t* block_weights) {
  hgp::mt_kahypar_get_block_weights(partitioned_hg, block_weights);
}

void mt_kahypar_get_graph_block_weights(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                        mt_kahypar_hypernode_weight_t* block_weights) {
  gp::mt_kahypar_get_block_weights(partitioned_graph, block_weights);
}

double mt_kahypar_hypergraph_imbalance(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                       const mt_kahypar_context_t* context) {
  return hgp::mt_kahypar_imbalance(partitioned_hg, context);
}

double mt_kahypar_graph_imbalance(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                            const mt_kahypar_context_t* context) {
  return gp::mt_kahypar_imbalance(partitioned_graph, context);
}

mt_kahypar_hyperedge_weight_t mt_kahypar_hypergraph_cut(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return hgp::mt_kahypar_cut(partitioned_hg);
}

mt_kahypar_hyperedge_weight_t mt_kahypar_graph_cut(const mt_kahypar_partitioned_graph_t* partitioned_graph) {
  return gp::mt_kahypar_cut(partitioned_graph);
}

mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return hgp::mt_kahypar_km1(partitioned_hg);
}

mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg) {
  return hgp::mt_kahypar_soed(partitioned_hg);
}