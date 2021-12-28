#include <memory>
#include <vector>
#include <iostream>

#include <libmtkahypar.h>

// Install library interface via 'sudo make install.mtkahypar' in build folder
// Compile with: g++ -std=c++14 -DNDEBUG -O3 partition_hypergraph.cc -o example_1 -lmtkahypard
// or flag "-lmtkahyparq", if you want to use our strong hypergraph partitioner
int main(int argc, char* argv[]) {

  // Initialize thread pool with 8 threads and NUMA allocation policy INTERLEAVED
  mt_kahypar_initialize_thread_pool(8, true /* activate interleaved NUMA allocation policy */ );

  // Load context from file
  mt_kahypar_context_t* context = mt_kahypar_context_new();
  // Use "../config/quality_preset.ini", if compiled with flag "-lmtkahyparq"
  mt_kahypar_configure_context_from_file(context, "../config/default_preset.ini");

  // Setup Hypergraph
  const mt_kahypar_hypernode_id_t num_vertices = 7;
  const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

  std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(4);

  // force the cut to contain hyperedge 0 and 2
  hyperedge_weights[0] = 1;  hyperedge_weights[1] = 1000;
  hyperedge_weights[2] = 1;  hyperedge_weights[3] = 1000;

  std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);

  hyperedge_indices[0] = 0; hyperedge_indices[1] = 2;
  hyperedge_indices[2] = 6; hyperedge_indices[3] = 9;
  hyperedge_indices[4] = 12;

  std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);

  // hypergraph from hMetis manual page 14
  hyperedges[0] = 0;  hyperedges[1] = 2;
  hyperedges[2] = 0;  hyperedges[3] = 1;
  hyperedges[4] = 3;  hyperedges[5] = 4;
  hyperedges[6] = 3;  hyperedges[7] = 4;
  hyperedges[8] = 6;  hyperedges[9] = 2;
  hyperedges[10] = 5; hyperedges[11] = 6;

  const double imbalance = 0.03;
  const mt_kahypar_partition_id_t k = 2;

  mt_kahypar_hyperedge_weight_t objective = 0;

  std::vector<mt_kahypar_partition_id_t> partition(num_vertices, -1);

  // Partition Hypergraph
  mt_kahypar_partition(num_vertices, num_hyperedges,
       	               imbalance, k, 0 /* seed */,
               	       nullptr /* unit vertex_weights */, hyperedge_weights.get(),
               	       hyperedge_indices.get(), hyperedges.get(),
       	               &objective, context, partition.data(),
                       true /* verbose output */ );

  // Print objective and block of each vertex
  std::cout << "Objective: " << objective << std::endl;
  for ( int i = 0; i != num_vertices; ++i ) {
    std::cout << "Vertex " << i << " = " << partition[i] << std::endl;
  }

  mt_kahypar_context_free(context);
}