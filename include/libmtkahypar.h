/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#ifndef LIBMTKAHYPAR_H
#define LIBMTKAHYPAR_H

#include <stddef.h>

#include "libmtkahypartypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Configurable parameters of the partitioning context.
 */
typedef enum {
  // number of blocks of the partition
  NUM_BLOCKS,
  // imbalance factor
  EPSILON,
  // objective function (either 'cut' or 'km1')
  OBJECTIVE,
  // seed for randomization
  SEED,
  // number of V-cycles
  NUM_VCYCLES,
  // disables or enables logging
  VERBOSE
} mt_kahypar_context_parameter_type_t;

/**
 * Supported objective functions.
 */
typedef enum {
  CUT, // TODO: add cut tests
  KM1
} mt_kahypar_objective_t;

/**
 * Preset types for partitioning context.
 */
typedef enum {
  // deterministic partitioning mode (corresponds to Mt-KaHyPar-SDet)
  DETERMINISTIC,
  // computes good partitions very fast (corresponds to Mt-KaHyPar-D)
  SPEED,
  // extends speed preset with flow-based refinement
  // -> computes high-quality partitions (corresponds to Mt-KaHyPar-D-F)
  HIGH_QUALITY
} mt_kahypar_preset_type_t;


// ####################### Setup Context #######################

/**
 * Creates a new empty partitioning context object.
 */
MT_KAHYPAR_API mt_kahypar_context_t* mt_kahypar_context_new();

/**
 * Deletes the partitioning context object.
 */
MT_KAHYPAR_API void mt_kahypar_free_context(mt_kahypar_context_t* context);

/**
 * Loads a partitioning context from a configuration file.
 */
MT_KAHYPAR_API void mt_kahypar_configure_context_from_file(mt_kahypar_context_t* context,
                                                           const char* ini_file_name);

/**
 * Loads a partitioning context of a predefined preset type.
 * Possible preset types are DETERMINISTIC (corresponds to Mt-KaHyPar-SDet),
 * SPEED (corresponds to Mt-KaHyPar-D) and HIGH_QUALITY (corresponds to Mt-KaHyPar-D-F)
 */
MT_KAHYPAR_API void mt_kahypar_load_preset(mt_kahypar_context_t* context,
                                           const mt_kahypar_preset_type_t preset);

/**
 * Sets a new value for a context parameter.
 *
 * Usage:
 * mt_kahypar_set_context_parameter(context, OBJECTIVE, "km1") // sets the objective function to the connectivity metric
 *
 * \return exit code zero if the corresponding parameter is successfully set to the value. Otherwise, it returns
 * 1 for an unknown parameter type, 2 for an integer conversion error or 3 for an unknown value type.
 */
MT_KAHYPAR_API int mt_kahypar_set_context_parameter(mt_kahypar_context_t* context,
                                                    const mt_kahypar_context_parameter_type_t type,
                                                    const char* value);

/**
 * Sets all required parameters for a partitioning call.
 */
MT_KAHYPAR_API void mt_kahypar_set_partitioning_parameters(mt_kahypar_context_t* context,
                                                           const mt_kahypar_partition_id_t num_blocks,
                                                           const double epsilon,
                                                           const mt_kahypar_objective_t objective,
                                                           const size_t seed);

/**
 * Sets individual target block weights for each block of the partition.
 * A balanced partition then satisfies that the weight of each block is smaller or equal than the
 * defined target block weight for the corresponding block.
 */
MT_KAHYPAR_API void mt_kahypar_set_individual_target_block_weights(mt_kahypar_context_t* context,
                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                   const mt_kahypar_hypernode_weight_t* block_weights);


// ####################### Thread Pool Initialization #######################

MT_KAHYPAR_API void mt_kahypar_initialize_thread_pool(const size_t num_threads,
                                                      const bool interleaved_allocations);

// ####################### Load/Construct Hypergraph #######################

/**
 * Reads a (hyper)graph from a file. The file can be either in hMetis or Metis file format.
 *
 * \note Note that for deterministic partitioning, you must call
 *       mt_kahypar_load_preset(context, DETERMINISTIC) before reading the hypergraph.
 * \note When reading a hMetis file with mt_kahypar_read_graph_from_file(...), make sure that
 *       the file represents graph. Otherwise, your program terminates.
 */
MT_KAHYPAR_API mt_kahypar_hypergraph_t* mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                                             const mt_kahypar_context_t* context,
                                                                             const mt_kahypar_file_format_type_t file_format);
MT_KAHYPAR_API mt_kahypar_graph_t* mt_kahypar_read_graph_from_file(const char* file_name,
                                                                   const mt_kahypar_context_t* context,
                                                                   const mt_kahypar_file_format_type_t file_format);

/**
 * Constructs a hypergraph from a given adjacency array that specifies the hyperedges.
 *
 * For example:
 * hyperedge_indices: | 0   | 2       | 6     | 9     | 12
 * hyperedges:        | 0 2 | 0 1 3 4 | 3 4 6 | 2 5 6 |
 * Defines a hypergraph with four hyperedges, e.g., e_0 = {0,2}, e_1 = {0,1,3,4}, ...
 *
 * \note For unweighted hypergraphs, you can pass nullptr to either hyperedge_weights or vertex_weights.
 * \note After construction, the arguments of this function are no longer needed and can be deleted.
 */
MT_KAHYPAR_API mt_kahypar_hypergraph_t* mt_kahypar_create_hypergraph(const mt_kahypar_hypernode_id_t num_vertices,
                                                                     const mt_kahypar_hyperedge_id_t num_hyperedges,
                                                                     const size_t* hyperedge_indices,
                                                                     const mt_kahypar_hyperedge_id_t* hyperedges,
                                                                     const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                                                     const mt_kahypar_hypernode_weight_t* vertex_weights);

/**
 * Constructs a graph from a given edge list vector.
 *
 * Example:
 * edges:        | 0 2 | 0 1 | 2 3 | 1 3 |
 * Defines a graph with four edges -> e_0 = {0,2}, e_1 = {0,1}, e_2 = {2,3}, e_3 = {1,3}
 *
 * \note For unweighted graphs, you can pass nullptr to either hyperedge_weights or vertex_weights.
 * \note After construction, the arguments of this function are no longer needed and can be deleted.
 */
MT_KAHYPAR_API mt_kahypar_graph_t* mt_kahypar_create_graph(const mt_kahypar_hypernode_id_t num_vertices,
                                                           const mt_kahypar_hyperedge_id_t num_edges,
                                                           const mt_kahypar_hypernode_id_t* edges,
                                                           const mt_kahypar_hyperedge_weight_t* edge_weights,
                                                           const mt_kahypar_hypernode_weight_t* vertex_weights);


/**
 * Deletes the (hyper)graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API void mt_kahypar_free_graph(mt_kahypar_graph_t* graph);

/**
 * Returns the number of nodes of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_hypernodes(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_graph_t* graph);

/**
 * Returns the number of (hyper)edges of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_edges(mt_kahypar_graph_t* graph);

/**
 * Returns the number of pins of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t* hypergraph);

/**
 * Returns the sum of all node weights of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_hypergraph_weight(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_graph_weight(mt_kahypar_graph_t* graph);

// ####################### Partition #######################

/**
 * Partitions a (hyper)graph according to the parameters specified in the partitioning context.
 *
 * \note Before partitioning, the number of blocks, imbalance parameter and objective function must be
 *       set in the partitioning context. This can be done either via mt_kahypar_set_context_parameter(...)
 *       or mt_kahypar_set_partitioning_parameters(...).
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                                    mt_kahypar_context_t* context);
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_partition_graph(mt_kahypar_graph_t* graph,
                                                                          mt_kahypar_context_t* context);

/**
 * Improves a given partition (using the V-cycle technique).
 *
 * \note The number of blocks specified in the partitioning context must be equal to the
 *       number of blocks of the given partition.
 * \note There is no guarantee that this call will find an improvement.
 */
MT_KAHYPAR_API void mt_kahypar_improve_hypergraph_partition(mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                            mt_kahypar_context_t* context,
                                                            const size_t num_vcycles);
MT_KAHYPAR_API void mt_kahypar_improve_graph_partition(mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                       mt_kahypar_context_t* context,
                                                       const size_t num_vcycles);

/**
 * Constructs a partitioned (hyper)graph out of the given partition.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                                             const mt_kahypar_partition_id_t num_blocks,
                                                                                             const mt_kahypar_partition_id_t* partition);
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_create_partitioned_graph(mt_kahypar_graph_t* graph,
                                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                                   const mt_kahypar_partition_id_t* partition);

/**
 * Constructs a partitioned (hyper)graph from a given partition file.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_read_hypergraph_partition_from_file(mt_kahypar_hypergraph_t* hypergraph,
                                                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                                                   const char* partition_file);
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_read_graph_partition_from_file(mt_kahypar_graph_t* graph,
                                                                                         const mt_kahypar_partition_id_t num_blocks,
                                                                                         const char* partition_file);

/**
 * Writes a partition to a file.
 */
MT_KAHYPAR_API void mt_kahypar_write_hypergraph_partition_to_file(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                                  const char* partition_file);
MT_KAHYPAR_API void mt_kahypar_write_graph_partition_to_file(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                             const char* partition_file);

/**
 * Extracts a partition from a partitioned (hyper)graph.
 */
MT_KAHYPAR_API void mt_kahypar_get_hypergraph_partition(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                        mt_kahypar_partition_id_t* partition);
MT_KAHYPAR_API void mt_kahypar_get_graph_partition(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                   mt_kahypar_partition_id_t* partition);

/**
 * Extracts the weight of each block from a partition.
 */
MT_KAHYPAR_API void mt_kahypar_get_hypergraph_block_weights(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                            mt_kahypar_hypernode_weight_t* block_weights);
MT_KAHYPAR_API void mt_kahypar_get_graph_block_weights(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                       mt_kahypar_hypernode_weight_t* block_weights);

/**
 * Computes the imbalance of the partition.
 */
MT_KAHYPAR_API double mt_kahypar_hypergraph_imbalance(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                      const mt_kahypar_context_t* context);
MT_KAHYPAR_API double mt_kahypar_graph_imbalance(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                 const mt_kahypar_context_t* context);

/**
 * Computes the cut metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_hypergraph_cut(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_graph_cut(const mt_kahypar_partitioned_graph_t* partitioned_graph);

/**
 * Computes the connectivity metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);

/**
 * Computes the sum-of-external-degree metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);


/**
 * Deletes the partitioned (hyper)graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg);
MT_KAHYPAR_API void mt_kahypar_free_partitioned_graph(mt_kahypar_partitioned_graph_t* partitioned_graph);

#ifdef __cplusplus
}
#endif

#endif    // LIBMTKAHYPAR_H