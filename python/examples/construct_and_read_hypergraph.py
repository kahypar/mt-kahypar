# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.

import mtkahypar
import multiprocessing

# Initialize
mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Create context
context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)

# Creates a weighted hypergraph
hypergraph = mtk.create_hypergraph(
  context,
  7,               # with seven nodes
  4,               # and four hyperedges
  [[0,2],          # The first hyperedge contains node 0 and 2
   [0,1,3,4],      # The second hyperedge contains node 0, 1, 3 and 4
   [3,4,6],        # The third hyperedge contains node 3, 4 and 6
   [2,5,6]],       # The fourth hyperedge contains node 2, 5 and 6
  [4,2,3,4,2,2,1], # node weights
  [1,2,1,3])       # hyperedge weights

# Output statistics of hypergraph
print("Hypergraph Stats:")
print("Number of Nodes      = " + str(hypergraph.num_nodes()))
print("Number of Hyperedges = " + str(hypergraph.num_edges()))
print("Number of Pins       = " + str(hypergraph.num_pins()))
print("Weight of Hypergraph = " + str(hypergraph.total_weight()))
print()

# Iterate over the nodes
for node in hypergraph.nodes():
  print("Node " + str(node))
  print("Node Degree = " + str(hypergraph.node_degree(node)))
  print("Node Weight = " + str(hypergraph.node_weight(node)))
  # Iterate over all incident hyperedges of the node
  print("Incident Hyperedges:")
  for edge in hypergraph.incident_edges(node):
    print(edge, end = " ")
  print()
print()

# Iterate over the hyperedges
for edge in hypergraph.edges():
  print("Hyperedge " + str(edge))
  print("Edge Size   = " + str(hypergraph.edge_size(edge)))
  print("Edge Weight = " + str(hypergraph.edge_weight(edge)))
  # Iterate over all pins of the edge
  print("Pins:")
  for pin in hypergraph.pins(edge):
    print(pin, end = " ")
  print()
print()

# Create a partition of the hypergraph
number_of_blocks = 3
partitioned_hg = hypergraph.create_partitioned_hypergraph(
  context,
  number_of_blocks, # partition hypergraph into three blocks
  [0,0,0,1,1,2,2])  # block IDs for each node

# Iterate over all nodes
for node in hypergraph.nodes():
  print("Node " + str(node))
  print("Block ID                     = " + str(partitioned_hg.block_id(node)))
  print("Number of Incident Cut Edges = " + str(partitioned_hg.num_incident_cut_edges(node)))
print()

# Iterate over all hyperedges
for edge in hypergraph.edges():
  print("Hyperedge " + str(edge))
  print("Connectivity = " + str(partitioned_hg.connectivity(edge))) # number of blocks contained in edge
  print("Edge Size    = " + str(hypergraph.edge_size(edge)))
  for block in partitioned_hg.connectivity_set(edge):
    print("Number of pins contained in block " + str(block) + " = "
      + str(partitioned_hg.num_pins_in_block(edge, block)))
