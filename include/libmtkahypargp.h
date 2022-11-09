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

#ifndef LIBMTKAHYPARGP_H
#define LIBMTKAHYPARGP_H

#include <stddef.h>

#include "include/libmtkahypartypes.h"

namespace gp {

// ####################### Load/Construct Graph #######################

/**
 * Reads a graph from a file. The file can be either in hMetis or Metis file format.
 *
 * \note Note that for deterministic partitioning, you must call
 *       mt_kahypar_load_preset(context, DETERMINISTIC) before reading the graph.
 */
MT_KAHYPAR_API mt_kahypar_graph_t* mt_kahypar_read_graph_from_file(const char* file_name,
                                                                   const mt_kahypar_context_t* context,
                                                                   const mt_kahypar_file_format_type_t file_format);

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
 * Deletes the graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_graph(mt_kahypar_graph_t* graph);

/**
 * Returns the number of nodes of the graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_graph_t* graph);

/**
 * Returns the number of edges of the graph.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_edges(mt_kahypar_graph_t* graph);

/**
 * Returns the sum of all node weights of the graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_total_weight(mt_kahypar_graph_t* graph);

// ####################### Partition #######################

/**
 * Partitions a graph according to the parameters specified in the partitioning context.
 *
 * \note Before partitioning, the number of blocks, imbalance parameter and objective function must be
 *       set in the partitioning context. This can be done either via mt_kahypar_set_context_parameter(...)
 *       or mt_kahypar_set_partitioning_parameters(...).
 */
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_partition(mt_kahypar_graph_t* graph,
                                                                    mt_kahypar_context_t* context);

/**
 * Improves a given partition (using the V-cycle technique).
 *
 * \note The number of blocks specified in the partitioning context must be equal to the
 *       number of blocks of the given partition.
 * \note There is no guarantee that this call will find an improvement.
 */
MT_KAHYPAR_API void mt_kahypar_improve_partition(mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                 mt_kahypar_context_t* context,
                                                 const size_t num_vcycles);

/**
 * Constructs a partitioned graph out of the given partition.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_create_partitioned_graph(mt_kahypar_graph_t* graph,
                                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                                   const mt_kahypar_partition_id_t* partition);

/**
 * Constructs a partitioned graph from a given partition file.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_graph_t* mt_kahypar_read_partition_from_file(mt_kahypar_graph_t* graph,
                                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                                   const char* partition_file);

/**
 * Writes a partition to a file.
 */
MT_KAHYPAR_API void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                       const char* partition_file);

/**
 * Extracts a partition from a partitioned graph.
 */
MT_KAHYPAR_API void mt_kahypar_get_partition(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                             mt_kahypar_partition_id_t* partition);

/**
 * Extracts the weight of each block from a partition.
 */
MT_KAHYPAR_API void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                                 mt_kahypar_hypernode_weight_t* block_weights);

/**
 * Computes the imbalance of the partition.
 */
MT_KAHYPAR_API double mt_kahypar_imbalance(const mt_kahypar_partitioned_graph_t* partitioned_graph,
                                           const mt_kahypar_context_t* context);

/**
 * Computes the cut metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_graph_t* partitioned_graph);


/**
 * Deletes the partitioned graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_partitioned_graph(mt_kahypar_partitioned_graph_t* partitioned_graph);


} // namespace gp

#endif  // LIBMTKAHYPARGP_H