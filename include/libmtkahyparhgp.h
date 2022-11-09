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

#ifndef LIBMTKAHYPARHGP_H
#define LIBMTKAHYPARHGP_H

#include <stddef.h>

#include "include/libmtkahypartypes.h"

namespace hgp {

// ####################### Load/Construct Hypergraph #######################

/**
 * Reads a hypergraph from a file. The file can be either in hMetis or Metis file format.
 *
 * \note Note that for deterministic partitioning, you must call
 *       mt_kahypar_load_preset(context, DETERMINISTIC) before reading the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypergraph_t* mt_kahypar_read_hypergraph_from_file(const char* file_name,
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
 * Deletes the hypergraph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t* hypergraph);

/**
 * Returns the number of nodes of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_hypergraph_t* hypergraph);

/**
 * Returns the number of hyperedges of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t* hypergraph);

/**
 * Returns the number of pins of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t* hypergraph);

/**
 * Returns the sum of all node weights of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_total_weight(mt_kahypar_hypergraph_t* hypergraph);

// ####################### Partition #######################

/**
 * Partitions a hypergraph according to the parameters specified in the partitioning context.
 *
 * \note Before partitioning, the number of blocks, imbalance parameter and objective function must be
 *       set in the partitioning context. This can be done either via mt_kahypar_set_context_parameter(...)
 *       or mt_kahypar_set_partitioning_parameters(...).
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition(mt_kahypar_hypergraph_t* hypergraph,
                                                                         mt_kahypar_context_t* context);

/**
 * Improves a given partition (using the V-cycle technique).
 *
 * \note The number of blocks specified in the partitioning context must be equal to the
 *       number of blocks of the given partition.
 * \note There is no guarantee that this call will find an improvement.
 */
MT_KAHYPAR_API void mt_kahypar_improve_partition(mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                 mt_kahypar_context_t* context,
                                                 const size_t num_vcycles);

/**
 * Constructs a partitioned hypergraph out of the given partition.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                                             const mt_kahypar_partition_id_t num_blocks,
                                                                                             const mt_kahypar_partition_id_t* partition);

/**
 * Constructs a partitioned hypergraph from a given partition file.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_read_partition_from_file(mt_kahypar_hypergraph_t* hypergraph,
                                                                                        const mt_kahypar_partition_id_t num_blocks,
                                                                                        const char* partition_file);

/**
 * Writes a partition to a file.
 */
MT_KAHYPAR_API void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                       const char* partition_file);

/**
 * Extracts a partition from a partitioned hypergraph.
 */
MT_KAHYPAR_API void mt_kahypar_get_partition(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                             mt_kahypar_partition_id_t* partition);

/**
 * Extracts the weight of each block from a partition.
 */
MT_KAHYPAR_API void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                 mt_kahypar_hypernode_weight_t* block_weights);

/**
 * Computes the imbalance of the partition.
 */
MT_KAHYPAR_API double mt_kahypar_imbalance(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                           const mt_kahypar_context_t* context);

/**
 * Computes the cut metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);

/**
 * Computes the connectivity metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);

/**
 * Computes the sum-of-external-degree metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);


/**
 * Deletes the partitioned hypergraph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg);


} // namespace hgp

#endif  // LIBMTKAHYPARHGP_H