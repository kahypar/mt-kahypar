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

#include "mtkahypartypes.h"

#ifdef __cplusplus
extern "C" {
#endif



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
                                                           const mt_kahypar_objective_t objective);


/**
 * Initializes the random number generator with the given seed value (not thread-safe).
 */
MT_KAHYPAR_API void mt_kahypar_set_seed(const size_t seed);

/**
 * Sets individual target block weights for each block of the partition.
 * A balanced partition then satisfies that the weight of each block is smaller or equal than the
 * defined target block weight for the corresponding block.
 */
MT_KAHYPAR_API void mt_kahypar_set_individual_target_block_weights(mt_kahypar_context_t* context,
                                                                   const mt_kahypar_partition_id_t num_blocks,
                                                                   const mt_kahypar_hypernode_weight_t* block_weights);


// ####################### Thread Pool Initialization #######################

MT_KAHYPAR_API void mt_kahypar_initialize(const size_t num_threads, const bool interleaved_allocations);

// ####################### Load/Construct Hypergraph #######################

/**
 * Reads a (hyper)graph from a file for a given configuration (preset).
 * The file can be either in hMetis or Metis file format.
 *
 * \note Note that we use different (hyper)graph data structures for different configurations.
 * Make sure that you partition the hypergraph with the same configuration as it is loaded.
 */
MT_KAHYPAR_API mt_kahypar_hypergraph_t mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                                            const mt_kahypar_preset_type_t preset,
                                                                            const mt_kahypar_file_format_type_t file_format);

/**
 * Reads a target graph in Metis file format. The target graph can be used in the
 * 'mt_kahypar_map' function to map a (hyper)graph onto it.
 */
MT_KAHYPAR_API mt_kahypar_target_graph_t* mt_kahypar_read_target_graph_from_file(const char* file_name);

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
MT_KAHYPAR_API mt_kahypar_hypergraph_t mt_kahypar_create_hypergraph(mt_kahypar_preset_type_t preset,
                                                                    const mt_kahypar_hypernode_id_t num_vertices,
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
MT_KAHYPAR_API mt_kahypar_hypergraph_t mt_kahypar_create_graph(const mt_kahypar_preset_type_t preset,
                                                               const mt_kahypar_hypernode_id_t num_vertices,
                                                               const mt_kahypar_hyperedge_id_t num_edges,
                                                               const mt_kahypar_hypernode_id_t* edges,
                                                               const mt_kahypar_hyperedge_weight_t* edge_weights,
                                                               const mt_kahypar_hypernode_weight_t* vertex_weights);
/**
 * Constructs a target graph from a given edge list vector. The target graph can be used in the
 * 'mt_kahypar_map' function to map a (hyper)graph onto it.
 *
 * Example:
 * edges:        | 0 2 | 0 1 | 2 3 | 1 3 |
 * Defines a graph with four edges -> e_0 = {0,2}, e_1 = {0,1}, e_2 = {2,3}, e_3 = {1,3}
 *
 * \note For unweighted graphs, you can pass nullptr to either hyperedge_weights.
 * \note After construction, the arguments of this function are no longer needed and can be deleted.
 */
MT_KAHYPAR_API mt_kahypar_target_graph_t* mt_kahypar_create_target_graph(const mt_kahypar_hypernode_id_t num_vertices,
                                                                         const mt_kahypar_hyperedge_id_t num_edges,
                                                                         const mt_kahypar_hypernode_id_t* edges,
                                                                         const mt_kahypar_hyperedge_weight_t* edge_weights);

/**
 * Deletes the (hyper)graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t hypergraph);

/**
 * Deletes a target graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_target_graph(mt_kahypar_target_graph_t* target_graph);

/**
 * Returns the number of nodes of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_hypernodes(mt_kahypar_hypergraph_t hypergraph);

/**
 * Returns the number of (hyper)edges of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t hypergraph);

/**
 * Returns the number of pins of the hypergraph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t hypergraph);

/**
 * Returns the sum of all node weights of the (hyper)graph.
 */
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_hypergraph_weight(mt_kahypar_hypergraph_t hypergraph);

// ####################### Fixed Vertices #######################

/**
 * Adds fixed vertices to the (hyper)graph as specified in the array to which 'fixed_vertices' points to.
 * The array should contain n entries (n = number of nodes). Each entry contains either the fixed vertex
 * block ID of the corresponding node or -1 if the node is not fixed to a block.
 */
MT_KAHYPAR_API void mt_kahypar_add_fixed_vertices(mt_kahypar_hypergraph_t hypergraph,
                                                  mt_kahypar_partition_id_t* fixed_vertices,
                                                  mt_kahypar_partition_id_t num_blocks);

/**
 * Reads fixed vertices from a file and stores them in the array to which 'fixed_vertices' points to.
 */
MT_KAHYPAR_API void mt_kahypar_read_fixed_vertices_from_file(const char* file_name,
                                                             mt_kahypar_partition_id_t* fixed_vertices);

/**
 * Adds fixed vertices to the (hyper)graph as specified in the fixed vertex file (expected in hMetis fix file format).
 * The file should contain n lines (n = number of nodes). Each line contains either the fixed vertex
 * block ID of the corresponding node or -1 if the node is not fixed to a block.
 */
MT_KAHYPAR_API void mt_kahypar_add_fixed_vertices_from_file(mt_kahypar_hypergraph_t hypergraph,
                                                            const char* file_name,
                                                            mt_kahypar_partition_id_t num_blocks);

/**
 * Removes all fixed vertices from the hypergraph.
 */
MT_KAHYPAR_API void mt_kahypar_remove_fixed_vertices(mt_kahypar_hypergraph_t hypergraph);

// ####################### Partition #######################

/**
 * Checks whether or not the given hypergraph can be partitioned with the corresponding preset.
 */
MT_KAHYPAR_API bool mt_kahypar_check_compatibility(mt_kahypar_hypergraph_t hypergraph,
                                                   mt_kahypar_preset_type_t preset);

/**
 * Partitions a (hyper)graph with the configuration specified in the partitioning context.
 *
 * \note Before partitioning, the number of blocks, imbalance parameter and objective function must be
 *       set in the partitioning context. This can be done either via mt_kahypar_set_context_parameter(...)
 *       or mt_kahypar_set_partitioning_parameters(...).
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t mt_kahypar_partition(mt_kahypar_hypergraph_t hypergraph,
                                                                        mt_kahypar_context_t* context);

/**
 * Maps a (hyper)graph onto a target graph with the configuration specified in the partitioning context.
 * The number of blocks of the output mapping/partition is the same as the number of nodes in the target graph
 * (each node of the target graph represents a block). The objective is to minimize the total weight of
 * all Steiner trees spanned by the (hyper)edges on the target graph. A Steiner tree is a tree with minimal weight
 * that spans a subset of the nodes (in our case the hyperedges) on the target graph. This objective function
 * is able to acurately model wire-lengths in VLSI design or communication costs in a distributed system where some
 * processors do not communicate directly with each other or different speeds.
 *
 * \note Since computing Steiner trees is an NP-hard problem, we currently restrict the size of the target graph
 * to at most 64 nodes. If you want to map hypergraphs onto larger target graphs, you can use recursive multisectioning.
 * For example, if the target graph has 4096 nodes, you can first map the hypergraph onto a coarser approximation of the
 * target graph with 64 nodes, and subsequently map each block of the mapping to the corresponding subgraph of the
 * target graph each having 64 nodes.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t mt_kahypar_map(mt_kahypar_hypergraph_t hypergraph,
                                                                  mt_kahypar_target_graph_t* target_graph,
                                                                  mt_kahypar_context_t* context);

/**
 * Checks whether or not the given partitioned hypergraph can
 * be improved with the corresponding preset.
 */
MT_KAHYPAR_API bool mt_kahypar_check_partition_compatibility(mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                                             mt_kahypar_preset_type_t preset);

/**
 * Improves a given partition (using the V-cycle technique).
 *
 * \note The number of blocks specified in the partitioning context must be equal to the
 *       number of blocks of the given partition.
 * \note There is no guarantee that this call will find an improvement.
 */
MT_KAHYPAR_API void mt_kahypar_improve_partition(mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                                 mt_kahypar_context_t* context,
                                                 const size_t num_vcycles);

/**
 * Improves a given mapping (using the V-cycle technique).
 *
 * \note The number of nodes of the target graph must be equal to the
 *       number of blocks of the given partition.
 * \note There is no guarantee that this call will find an improvement.
 */
MT_KAHYPAR_API void mt_kahypar_improve_mapping(mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                               mt_kahypar_target_graph_t* target_graph,
                                               mt_kahypar_context_t* context,
                                               const size_t num_vcycles);

/**
 * Constructs a partitioned (hyper)graph out of the given partition.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t hypergraph,
                                                                                            const mt_kahypar_preset_type_t preset,
                                                                                            const mt_kahypar_partition_id_t num_blocks,
                                                                                            const mt_kahypar_partition_id_t* partition);

/**
 * Constructs a partitioned (hyper)graph from a given partition file.
 */
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t mt_kahypar_read_partition_from_file(mt_kahypar_hypergraph_t hypergraph,
                                                                                       const mt_kahypar_preset_type_t preset,
                                                                                       const mt_kahypar_partition_id_t num_blocks,
                                                                                       const char* partition_file);

/**
 * Writes a partition to a file.
 */
MT_KAHYPAR_API void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                                       const char* partition_file);

/**
 * Extracts a partition from a partitioned (hyper)graph.
 */
MT_KAHYPAR_API void mt_kahypar_get_partition(const mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                             mt_kahypar_partition_id_t* partition);

/**
 * Extracts the weight of each block from a partition.
 */
MT_KAHYPAR_API void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                                 mt_kahypar_hypernode_weight_t* block_weights);
/**
 * Computes the imbalance of the partition.
 */
MT_KAHYPAR_API double mt_kahypar_imbalance(const mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                           const mt_kahypar_context_t* context);

/**
 * Computes the cut metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_hypergraph_t partitioned_hg);

/**
 * Computes the connectivity metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t partitioned_hg);

/**
 * Computes the sum-of-external-degree metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t partitioned_hg);

/**
 * Computes the steiner tree metric.
 */
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_steiner_tree(const mt_kahypar_partitioned_hypergraph_t partitioned_hg,
                                                                     mt_kahypar_target_graph_t* target_graph);


/**
 * Deletes the partitioned (hyper)graph object.
 */
MT_KAHYPAR_API void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t partitioned_hg);

#ifdef __cplusplus
}
#endif

#endif    // LIBMTKAHYPAR_H