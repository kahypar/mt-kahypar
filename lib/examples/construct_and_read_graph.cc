#include <memory>
#include <vector>
#include <iostream>
#include <thread>

#include <mtkahypar.h>

// Install library interface via 'sudo make install.mtkahypar' in build folder
// Compile with: g++ -std=c++14 -DNDEBUG -O3 construct_and_read_graph.cc -o example -lmtkahypar
int main(int argc, char* argv[]) {
  mt_kahypar_error_t error{};

  // Initialize thread pool
  mt_kahypar_initialize(
    std::thread::hardware_concurrency() /* use all available cores */,
    true /* activate interleaved NUMA allocation policy */ );

  // In the following, we construct a graph with 5 nodes and 6 edges
  const mt_kahypar_hypernode_id_t num_nodes = 5;
  const mt_kahypar_hyperedge_id_t num_edges = 6;

  // We represent the edges of the graph as edge list vector.
  // Two consecutive node IDs in the edge list vector form
  // an undirected edge in the graph.
  std::unique_ptr<mt_kahypar_hypernode_id_t[]> edges =
    std::make_unique<mt_kahypar_hypernode_id_t[]>(12);
  edges[0] = 0;  edges[1] = 1;  // first edge connects node 0 and 1
  edges[2] = 0;  edges[3] = 2;  // second edge connects node 0 and 2
  edges[4] = 1;  edges[5] = 2;  // third edge connects node 1 and 2
  edges[6] = 1;  edges[7] = 3;  // fourth edge connects node 1 and 3
  edges[8] = 2;  edges[9] = 3;  // fifth edge connects node 2 and 3
  edges[10] = 3; edges[11] = 4; // sixth edge connects node 3 and 4

  // Define node weights
  std::unique_ptr<mt_kahypar_hypernode_weight_t[]> node_weights =
    std::make_unique<mt_kahypar_hypernode_weight_t[]>(5);
  node_weights[0] = 2; node_weights[1] = 1; node_weights[2] = 2;
  node_weights[3] = 4; node_weights[4] = 1;

  // Define edge weights
  std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edge_weights =
    std::make_unique<mt_kahypar_hyperedge_weight_t[]>(6);
  edge_weights[0] = 1; edge_weights[1] = 10;
  edge_weights[2] = 1; edge_weights[3] = 10;
  edge_weights[4] = 1; edge_weights[5] = 10;

  // Construct graph
  mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
  mt_kahypar_hypergraph_t graph =
    mt_kahypar_create_graph(context, num_nodes, num_edges,
      edges.get(), edge_weights.get(), node_weights.get(), &error);
  if (graph.hypergraph == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  std::cout << "Number of Nodes       = " << mt_kahypar_num_hypernodes(graph) << std::endl;
  std::cout << "Number of Edges       = " << mt_kahypar_num_hyperedges(graph) << std::endl;
  std::cout << "Total Weight of Graph = " << mt_kahypar_hypergraph_weight(graph) << std::endl;

  // Iterate over the nodes
  std::vector<mt_kahypar_hyperedge_id_t> incident_edges_buffer;
  for (mt_kahypar_hypernode_id_t node = 0; node < mt_kahypar_num_hypernodes(graph); ++node) {
    std::cout << "Node " << node << std::endl;
    std::cout << "Node Degree = " << mt_kahypar_hypernode_degree(graph, node) << std::endl;
    std::cout << "Node Weight = " << mt_kahypar_hypernode_weight(graph, node) << std::endl;

    // Get incident edges of the node by writing them to the buffer
    incident_edges_buffer.resize(mt_kahypar_hypernode_degree(graph, node));
    mt_kahypar_get_incident_hyperedges(graph, node, incident_edges_buffer.data());

    // Iterate over all neighbors
    std::cout << "Neighbors: ";
    for (mt_kahypar_hyperedge_id_t edge: incident_edges_buffer) {
      mt_kahypar_hypernode_id_t neighbor = mt_kahypar_edge_target(graph, edge);
      std::cout << neighbor << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Iterate over the edges
  // Note that the graph has six undirected edges, but internally we represent
  // the graph as a directed graph. Thus, we iterate over 12 directed edges here.
  // For an directed edge (u,v), we call u the source and v the target node.
  for (mt_kahypar_hyperedge_id_t edge = 0; edge < mt_kahypar_num_hyperedges(graph); ++edge) {
    std::cout << "Edge " << edge << std::endl;
    std::cout << "Edge Weight = " << mt_kahypar_hyperedge_weight(graph, edge) << std::endl;
    std::cout << "Source Node = " << mt_kahypar_edge_source(graph, edge) << std::endl;
    std::cout << "Target Node = " << mt_kahypar_edge_target(graph, edge) << std::endl;
  }
  std::cout << std::endl;

  // Create a partition of the graph
  const mt_kahypar_partition_id_t number_of_blocks = 3;
  const std::vector<mt_kahypar_partition_id_t> partition = {0, 0, 1, 1, 2};
  mt_kahypar_partitioned_hypergraph_t partitioned_graph =
    mt_kahypar_create_partitioned_hypergraph(graph, context, number_of_blocks, partition.data(), &error);
  if (partitioned_graph.partitioned_hg == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  // Iterate over all nodes
  for (mt_kahypar_hypernode_id_t node = 0; node < mt_kahypar_num_hypernodes(graph); ++node) {
    std::cout << "Node " << node << std::endl;
    std::cout << "Block ID                     = "
      << mt_kahypar_block_id(partitioned_graph, node) << std::endl;
    std::cout << "Number of Incident Cut Edges = "
      << mt_kahypar_num_incident_cut_hyperedges(partitioned_graph, node) << std::endl;
  }
  std::cout << std::endl;

  // Iterate over all edges
  for (mt_kahypar_hyperedge_id_t edge = 0; edge < mt_kahypar_num_hyperedges(graph); ++edge) {
    std::cout << "Edge " << edge << std::endl;
    // Print number of blocks connected by edge
    std::cout << "Connectivity = " << mt_kahypar_connectivity(partitioned_graph, edge) << std::endl;
  }
  std::cout << std::endl;

  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(graph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_graph);
}
