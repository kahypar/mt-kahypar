#include <stdio.h>
#include <assert.h>

#include <mtkahypar.h>

// Install library interface via 'sudo make install.mtkahypar' in build folder
// Compile with: gcc -std=c99 -DNDEBUG -O3 partition_c_example.c -o example -lmtkahypar
int main(int argc, char* argv[]) {
  mt_kahypar_error_t error = {0};

  // Initialize
  mt_kahypar_initialize(
    4, /* use 4 threads */
    true /* activate interleaved NUMA allocation policy */ );

  // Setup partitioning context
  mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
  // In the following, we partition a graph into two blocks
  // with an allowed imbalance of 3% and optimize the edge cut (CUT)
  mt_kahypar_set_partitioning_parameters(context,
    2 /* number of blocks */, 0.03 /* imbalance parameter */,
    CUT /* objective function */);
  mt_kahypar_set_seed(42 /* seed */);
  // Enable logging
  mt_kahypar_status_t status =
    mt_kahypar_set_context_parameter(context, VERBOSE, "1", &error);
  assert(status == SUCCESS);

  // Read Graph
  mt_kahypar_hypergraph_t graph = mt_kahypar_read_hypergraph_from_file(
    "delaunay_n15.graph", context, METIS /* file format */, &error);
  if (graph.hypergraph == NULL) {
    printf("%s\n", error.msg); return 1;
  }

  // Partition Graph
  mt_kahypar_partitioned_hypergraph_t partitioned_graph =
    mt_kahypar_partition(graph, context, &error);
  if (partitioned_graph.partitioned_hg == NULL) {
    printf("%s\n", error.msg); return 1;
  }

  // Compute Metrics
  const double imbalance = mt_kahypar_imbalance(partitioned_graph, context);
  const int cut = mt_kahypar_cut(partitioned_graph);

  // Output Results
  printf("Partitioning Results:\n");
  printf("Imbalance         = %f\n", imbalance);
  printf("Cut               = %d\n", cut);

  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(graph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_graph);
}
