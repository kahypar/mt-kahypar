#include <memory>
#include <vector>
#include <iostream>
#include <thread>

#include <mtkahypar.h>

// Install library interface via 'sudo make install.mtkahypar' in build folder
// Compile with: g++ -std=c++14 -DNDEBUG -O3 construct_and_read_hypergraph.cc -o example -lmtkahypar
int main(int argc, char* argv[]) {
  mt_kahypar_error_t error{};

  // Initialize thread pool
  mt_kahypar_initialize(
    std::thread::hardware_concurrency() /* use all available cores */,
    true /* activate interleaved NUMA allocation policy */ );

  // In the following, we construct a hypergraph with 7 nodes and 4 hyperedges
  const mt_kahypar_hypernode_id_t num_nodes = 7;
  const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

  // The hyperedge indices points to the hyperedge vector and defines the
  // the ranges containing the pins of each hyperedge
  std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);
  hyperedge_indices[0] = 0; hyperedge_indices[1] = 2; hyperedge_indices[2] = 6;
  hyperedge_indices[3] = 9; hyperedge_indices[4] = 12;

  std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges =
    std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);
  // First hyperedge contains nodes with ID 0 and 2
  hyperedges[0] = 0;  hyperedges[1] = 2;
  // Second hyperedge contains nodes with ID 0, 1, 3 and 4
  hyperedges[2] = 0;  hyperedges[3] = 1; hyperedges[4] = 3;  hyperedges[5] = 4;
  // Third hyperedge contains nodes with ID 3, 4 and 6
  hyperedges[6] = 3;  hyperedges[7] = 4; hyperedges[8] = 6;
  // Fourth hyperedge contains nodes with ID 2, 5 and 6
  hyperedges[9] = 2; hyperedges[10] = 5; hyperedges[11] = 6;

  // Define node weights
  std::unique_ptr<mt_kahypar_hypernode_weight_t[]> node_weights =
    std::make_unique<mt_kahypar_hypernode_weight_t[]>(7);
  node_weights[0] = 2; node_weights[1] = 1; node_weights[2] = 2; node_weights[3] = 4;
  node_weights[4] = 1; node_weights[5] = 3; node_weights[6] = 3;

  // Define hyperedge weights
  std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> hyperedge_weights =
    std::make_unique<mt_kahypar_hyperedge_weight_t[]>(4);
  hyperedge_weights[0] = 1; hyperedge_weights[1] = 10;
  hyperedge_weights[2] = 1; hyperedge_weights[3] = 10;

  // Construct hypergraph
  mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar_create_hypergraph(context, num_nodes, num_hyperedges,
      hyperedge_indices.get(), hyperedges.get(),
      hyperedge_weights.get(), node_weights.get(), &error);
  if (hypergraph.hypergraph == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  std::cout << "Number of Nodes            = " << mt_kahypar_num_hypernodes(hypergraph) << std::endl;
  std::cout << "Number of Hyperedges       = " << mt_kahypar_num_hyperedges(hypergraph) << std::endl;
  std::cout << "Number of Pins             = " << mt_kahypar_num_pins(hypergraph) << std::endl;
  std::cout << "Total Weight of Hypergraph = " << mt_kahypar_hypergraph_weight(hypergraph) << std::endl;

  // Iterate over the nodes
  std::vector<mt_kahypar_hyperedge_id_t> incident_edges_buffer;
  for (mt_kahypar_hypernode_id_t node = 0; node < mt_kahypar_num_hypernodes(hypergraph); ++node) {
    std::cout << "Node " << node << std::endl;
    std::cout << "Node Degree = " << mt_kahypar_hypernode_degree(hypergraph, node) << std::endl;
    std::cout << "Node Weight = " << mt_kahypar_hypernode_weight(hypergraph, node) << std::endl;

    // Get incident hyperedges of the node by writing them to the buffer
    incident_edges_buffer.resize(mt_kahypar_hypernode_degree(hypergraph, node));
    mt_kahypar_get_incident_hyperedges(hypergraph, node, incident_edges_buffer.data());

    // Iterate over all incident hyperedges
    std::cout << "Incident Hyperedges: ";
    for (mt_kahypar_hyperedge_id_t edge: incident_edges_buffer) {
      std::cout << edge << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Iterate over the hyperedges
  std::vector<mt_kahypar_hypernode_id_t> pins_buffer;
  for (mt_kahypar_hyperedge_id_t edge = 0; edge < mt_kahypar_num_hyperedges(hypergraph); ++edge) {
    std::cout << "Edge " << edge << std::endl;
    std::cout << "Edge Size = " << mt_kahypar_hyperedge_size(hypergraph, edge) << std::endl;
    std::cout << "Edge Weight = " << mt_kahypar_hyperedge_weight(hypergraph, edge) << std::endl;

    // Get pins of the hyperedge by writing them to the buffer
    pins_buffer.resize(mt_kahypar_hyperedge_size(hypergraph, edge));
    mt_kahypar_get_hyperedge_pins(hypergraph, edge, pins_buffer.data());

    // Iterate over all pins of the hyperedge
    std::cout << "Pins: ";
    for (mt_kahypar_hypernode_id_t pin: pins_buffer) {
      std::cout << pin << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Create a partition of the graph
  const mt_kahypar_partition_id_t number_of_blocks = 3;
  const std::vector<mt_kahypar_partition_id_t> partition = {0, 0 ,0, 1, 1, 2, 2};
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_create_partitioned_hypergraph(hypergraph, context, number_of_blocks, partition.data(), &error);
  if (partitioned_hg.partitioned_hg == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  // Iterate over all nodes
  for (mt_kahypar_hypernode_id_t node = 0; node < mt_kahypar_num_hypernodes(hypergraph); ++node) {
    std::cout << "Node " << node << std::endl;
    std::cout << "Block ID                     = "
      << mt_kahypar_block_id(partitioned_hg, node) << std::endl;
    std::cout << "Number of Incident Cut Edges = "
      << mt_kahypar_num_incident_cut_hyperedges(partitioned_hg, node) << std::endl;
  }
  std::cout << std::endl;

  // Iterate over all edges
  for (mt_kahypar_hyperedge_id_t edge = 0; edge < mt_kahypar_num_hyperedges(hypergraph); ++edge) {
    std::cout << "Edge " << edge << std::endl;
    // Print number of blocks connected by edge
    std::cout << "Connectivity = " << mt_kahypar_connectivity(partitioned_hg, edge) << std::endl;
  }
  std::cout << std::endl;

  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(hypergraph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}
