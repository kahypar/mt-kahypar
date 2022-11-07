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

#ifndef LIBKAHYPAR_H
#define LIBKAHYPAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MT_KAHYPAR_API
#   if __GNUC__ >= 4
#       define MT_KAHYPAR_API __attribute__ ((visibility("default")))
#   else
#       define MT_KAHYPAR_API
#   endif
#endif

struct mt_kahypar_context_s;
typedef struct mt_kahypar_context_s mt_kahypar_context_t;
typedef struct mt_kahypar_hypergraph_s mt_kahypar_hypergraph_t;
typedef struct mt_kahypar_partitioned_hypergraph_s mt_kahypar_partitioned_hypergraph_t;

typedef unsigned long int mt_kahypar_hypernode_id_t;
typedef unsigned long int mt_kahypar_hyperedge_id_t;
typedef int mt_kahypar_hypernode_weight_t;
typedef int mt_kahypar_hyperedge_weight_t;
typedef int mt_kahypar_partition_id_t;

typedef enum {
  NUM_BLOCKS,
  EPSILON,
  OBJECTIVE, /** either 'cut' or 'km1' **/
  SEED,
  NUM_VCYCLES,
  VERBOSE
} mt_kahypar_context_parameter_type_t;

typedef enum {
  CUT, // TODO: add cut tests
  KM1
} mt_kahypar_objective_t;

typedef enum {
  DETERMINISTIC,
  SPEED,
  HIGH_QUALITY
} mt_kahypar_preset_type_t;

typedef enum {
  METIS, /** for graph files **/
  HMETIS /** for hypergraph files **/
} mt_kahypar_file_format_type_t;


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
 * Returns exit code zero if the corresponding parameter is successfully set to the value. Otherwise, it returns
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

MT_KAHYPAR_API mt_kahypar_hypergraph_t* mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                                             const mt_kahypar_context_t* context,
                                                                             const mt_kahypar_file_format_type_t file_format);
MT_KAHYPAR_API mt_kahypar_hypergraph_t* mt_kahypar_create_hypergraph(const mt_kahypar_hypernode_id_t num_vertices,
                                                                     const mt_kahypar_hyperedge_id_t num_hyperedges,
                                                                     const size_t* hyperedge_indices,
                                                                     const mt_kahypar_hyperedge_id_t* hyperedges,
                                                                     const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                                                     const mt_kahypar_hypernode_weight_t* vertex_weights);
MT_KAHYPAR_API void mt_kahypar_free_hypergraph(mt_kahypar_hypergraph_t* hypergraph);


MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_nodes(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hyperedge_id_t mt_kahypar_num_hyperedges(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_num_pins(mt_kahypar_hypergraph_t* hypergraph);
MT_KAHYPAR_API mt_kahypar_hypernode_id_t mt_kahypar_total_weight(mt_kahypar_hypergraph_t* hypergraph);

// ####################### Partition #######################

MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_partition(mt_kahypar_hypergraph_t* hypergraph,
                                                                         mt_kahypar_context_t* context);

MT_KAHYPAR_API void mt_kahypar_improve_partition(mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                 mt_kahypar_context_t* context,
                                                 const size_t num_vcycles);

MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_create_partitioned_hypergraph(mt_kahypar_hypergraph_t* hypergraph,
                                                                                             const mt_kahypar_partition_id_t num_blocks,
                                                                                             const mt_kahypar_partition_id_t* partition);
MT_KAHYPAR_API mt_kahypar_partitioned_hypergraph_t* mt_kahypar_read_partition_from_file(mt_kahypar_hypergraph_t* hypergraph,
                                                                                        const mt_kahypar_partition_id_t num_blocks,
                                                                                        const char* partition_file);
MT_KAHYPAR_API void mt_kahypar_write_partition_to_file(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                       const char* partition_file);

MT_KAHYPAR_API void mt_kahypar_get_partition(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                             mt_kahypar_partition_id_t* partition);
MT_KAHYPAR_API void mt_kahypar_get_block_weights(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                                 mt_kahypar_hypernode_weight_t* block_weights);
MT_KAHYPAR_API double mt_kahypar_imbalance(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg,
                                           const mt_kahypar_context_t* context);
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_cut(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_km1(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);
MT_KAHYPAR_API mt_kahypar_hyperedge_weight_t mt_kahypar_soed(const mt_kahypar_partitioned_hypergraph_t* partitioned_hg);


MT_KAHYPAR_API void mt_kahypar_free_partitioned_hypergraph(mt_kahypar_partitioned_hypergraph_t* partitioned_hg);


#ifdef __cplusplus
}
#endif

#endif    // LIBKAHYPAR_H