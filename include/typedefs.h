#ifndef TYPEDEFS_H
#define TYPEDEFS_H

struct mt_kahypar_context_s;
typedef struct mt_kahypar_context_s mt_kahypar_context_t;
typedef struct mt_kahypar_hypergraph_s mt_kahypar_hypergraph_t;
typedef struct mt_kahypar_graph_s mt_kahypar_graph_t;
typedef struct mt_kahypar_partitioned_hypergraph_s mt_kahypar_partitioned_hypergraph_t;
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