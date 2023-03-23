#ifndef TYPEDEFS_H
#define TYPEDEFS_H

typedef enum {
  STATIC_GRAPH,
  STATIC_HYPERGRAPH,
  DYNAMIC_GRAPH,
  DYNAMIC_HYPERGRAPH,
  NULLPTR_HYPERGRAPH
} mt_kahypar_hypergraph_type_t;

typedef enum {
  STATIC_PARTITIONED_GRAPH,
  STATIC_PARTITIONED_HYPERGRAPH,
  STATIC_SPARSE_PARTITIONED_HYPERGRAPH,
  DYNAMIC_PARTITIONED_GRAPH,
  DYNAMIC_PARTITIONED_HYPERGRAPH,
  NULLPTR_PARTITION
} mt_kahypar_partition_type_t;

struct mt_kahypar_context_s;
typedef struct mt_kahypar_context_s mt_kahypar_context_t;

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


typedef struct mt_kahypar_graph_s mt_kahypar_graph_t;
typedef struct mt_kahypar_partitioned_graph_s mt_kahypar_partitioned_graph_t;

typedef unsigned long int mt_kahypar_hypernode_id_t;
typedef unsigned long int mt_kahypar_hyperedge_id_t;
typedef int mt_kahypar_hypernode_weight_t;
typedef int mt_kahypar_hyperedge_weight_t;
typedef int mt_kahypar_partition_id_t;

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

#endif // TYPEDEFS_H