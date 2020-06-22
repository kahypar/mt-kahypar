/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#ifndef LIBKAHYPAR_H
#define LIBKAHYPAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef KAHYPAR_API
#  ifdef _WIN32
#     if defined(KAHYPAR_BUILD_SHARED)  /* build dll */
#         define KAHYPAR_API __declspec(dllexport)
#     elif !defined(KAHYPAR_BUILD_STATIC)  /* use dll */
#         define KAHYPAR_API __declspec(dllimport)
#     else  /* static library */
#         define KAHYPAR_API
#     endif
#  else
#     if __GNUC__ >= 4
#         define KAHYPAR_API __attribute__ ((visibility("default")))
#     else
#         define KAHYPAR_API
#     endif
#  endif
#endif

struct mt_kahypar_context_s;
typedef struct mt_kahypar_context_s mt_kahypar_context_t;

typedef unsigned long int mt_kahypar_hypernode_id_t;
typedef unsigned long int mt_kahypar_hyperedge_id_t;
typedef int mt_kahypar_hypernode_weight_t;
typedef int mt_kahypar_hyperedge_weight_t;
typedef unsigned int mt_kahypar_partition_id_t;

KAHYPAR_API mt_kahypar_context_t* mt_kahypar_context_new();
KAHYPAR_API void mt_kahypar_context_free(mt_kahypar_context_t* kahypar_context);
KAHYPAR_API void mt_kahypar_configure_context_from_file(mt_kahypar_context_t* kahypar_context,
                                                        const char* ini_file_name);

KAHYPAR_API void mt_kahypar_initialize_thread_pool(const size_t num_threads,
                                                   const bool interleaved_allocations);

KAHYPAR_API void mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                                      mt_kahypar_hypernode_id_t* num_vertices,
                                                      mt_kahypar_hyperedge_id_t* num_hyperedges,
                                                      size_t** hyperedge_indices,
                                                      mt_kahypar_hyperedge_id_t** hyperedges,
                                                      mt_kahypar_hyperedge_weight_t** hyperedge_weights,
                                                      mt_kahypar_hypernode_weight_t** vertex_weights);

KAHYPAR_API void mt_kahypar_partition(const mt_kahypar_hypernode_id_t num_vertices,
                                      const mt_kahypar_hyperedge_id_t num_hyperedges,
                                      const double epsilon,
                                      const mt_kahypar_partition_id_t num_blocks,
                                      const int seed,
                                      const mt_kahypar_hypernode_weight_t* vertex_weights,
                                      const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                                      const size_t* hyperedge_indices,
                                      const mt_kahypar_hyperedge_id_t* hyperedges,
                                      mt_kahypar_hyperedge_weight_t* objective,
                                      mt_kahypar_context_t* kahypar_context,
                                      mt_kahypar_partition_id_t* partition,
                                      const bool verbose = false);


#ifdef __cplusplus
}
#endif

#endif    // LIBKAHYPAR_H