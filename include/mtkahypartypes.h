#ifndef MTKAHYPAR_TYPEDEFS_H
#define MTKAHYPAR_TYPEDEFS_H

#include <stddef.h>

typedef enum {
  STATIC_GRAPH,
  DYNAMIC_GRAPH,
  STATIC_HYPERGRAPH,
  DYNAMIC_HYPERGRAPH,
  NULLPTR_HYPERGRAPH
} mt_kahypar_hypergraph_type_t;

typedef enum {
  MULTILEVEL_GRAPH_PARTITIONING,
  N_LEVEL_GRAPH_PARTITIONING,
  MULTILEVEL_HYPERGRAPH_PARTITIONING,
  N_LEVEL_HYPERGRAPH_PARTITIONING,
  LARGE_K_PARTITIONING,
  NULLPTR_PARTITION
} mt_kahypar_partition_type_t;

/**
 * Either success or type of the error.
 */
typedef enum {
  SUCCESS = 0,
  // input files are not found or an input is syntactically or semantically invalid
  INVALID_INPUT,
  // an algorithm parameter has an invalid value
  INVALID_PARAMETER,
  // the attempted operation is incompatible with the config,
  // e.g. the hypergraph type is incompatible with the preset type
  UNSUPPORTED_OPERATION,
  // errors originating from the OS, e.g. out of memory or failed mmap
  SYSTEM_ERROR,
  OTHER_ERROR
} mt_kahypar_status_t;

/**
 * Indicates whether an error occured.
 *
 * If an error occurs, mt_kahypar_free_error_content needs to be called
 * afterwards to free the space allocated for the error message.
 */
typedef struct {
  // null-terminated error message
  const char* msg;
  // length of error message without null terminator
  size_t msg_len;
  // either success or the error type
  mt_kahypar_status_t status;
} mt_kahypar_error_t;

struct mt_kahypar_context_s;
typedef struct mt_kahypar_context_s mt_kahypar_context_t;
struct mt_kahypar_target_graph_s;
typedef struct mt_kahypar_target_graph_s mt_kahypar_target_graph_t;

struct mt_kahypar_hypergraph_s;
typedef struct {
  mt_kahypar_hypergraph_s* hypergraph;
  mt_kahypar_hypergraph_type_t type;
} mt_kahypar_hypergraph_t;

typedef struct {
  const mt_kahypar_hypergraph_s* hypergraph;
  mt_kahypar_hypergraph_type_t type;
} mt_kahypar_hypergraph_const_t;

struct mt_kahypar_partitioned_hypergraph_s;
typedef struct {
  mt_kahypar_partitioned_hypergraph_s* partitioned_hg;
  mt_kahypar_partition_type_t type;
} mt_kahypar_partitioned_hypergraph_t;

typedef struct {
  const mt_kahypar_partitioned_hypergraph_s* partitioned_hg;
  mt_kahypar_partition_type_t type;
} mt_kahypar_partitioned_hypergraph_const_t;

typedef unsigned long int mt_kahypar_hypernode_id_t;
typedef unsigned long int mt_kahypar_hyperedge_id_t;
typedef int mt_kahypar_hypernode_weight_t;
typedef int mt_kahypar_hyperedge_weight_t;
typedef int mt_kahypar_partition_id_t;

/**
 * Configurable parameters of the partitioning context.
 */
typedef enum {
  // number of blocks of the partition (integer)
  NUM_BLOCKS,
  // imbalance factor (float)
  EPSILON,
  // objective function (either 'cut', 'km1' or 'soed')
  OBJECTIVE,
  // number of V-cycles (integer)
  NUM_VCYCLES,
  // enables logging (bool: 1/0)
  VERBOSE
} mt_kahypar_context_parameter_type_t;

/**
 * Supported objective functions.
 */
typedef enum {
  CUT,
  KM1,
  SOED
} mt_kahypar_objective_t;

/**
 * Preset types for partitioning context.
 */
typedef enum {
  // deterministic partitioning mode (corresponds to Mt-KaHyPar-SDet)
  DETERMINISTIC,
  // partitioning mode for partitioning a (hyper)graph into a large number of blocks
  LARGE_K,
  // computes good partitions very fast (corresponds to Mt-KaHyPar-D)
  DEFAULT,
  // extends default preset with flow-based refinement
  // -> computes high-quality partitions (corresponds to Mt-KaHyPar-D-F)
  QUALITY,
  // n-level code with flow-based refinement
  // => highest quality configuration (corresponds to Mt-KaHyPar-Q-F)
  HIGHEST_QUALITY
} mt_kahypar_preset_type_t;

/**
 * Supported (hyper)graph file formats.
 */
typedef enum {
  // Standard file format for graphs
  METIS,
  // Standard file format for hypergraphs
  HMETIS
} mt_kahypar_file_format_type_t;

#ifndef MT_KAHYPAR_API
#   if __GNUC__ >= 4
#       define MT_KAHYPAR_API __attribute__ ((visibility("default")))
#   else
#       define MT_KAHYPAR_API
#   endif
#endif

#endif // MTKAHYPAR_TYPEDEFS_H
