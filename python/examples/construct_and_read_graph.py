# Please follow the instructions in the README to install the python library interface
# and copy mtkahypar.so to the folder containing the examples.

import mtkahypar
import multiprocessing

# Initialize
mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Create context
context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)

# Creates a weighted graph
graph = mtk.create_graph(
  context,
  5,             # with 5 nodes
  6,             # and 6 undirected edges
  # Define six undirected edges
  [(0,1),        # first edge connects node 0 and 1
   (0,2),        # second edge connects node 0 and 2
   (1,2),        # third edge connects node 1 and 2
   (1,3),        # fourth edge connects node 1 and 3
   (2,3),        # fifth edge connects node 2 and 3
   (3,4)],       # sixth edge connects node 3 and 4
  [4,2,4,2,1],   # node weights
  [1,2,1,3,3,2]) # edge weights

# Output statistics of graph
print("Graph Stats:")
print("Number of Nodes            = " + str(graph.num_nodes()))
print("Number of Undirected Edges = " + str(graph.num_edges()))
print("Number of Directed Edges   = " + str(graph.num_directed_edges()))
print("Weight of Graph            = " + str(graph.total_weight()))
print()

# Iterate over the nodes
for node in graph.nodes():
  print("Node " + str(node))
  print("Node Degree = " + str(graph.node_degree(node)))
  print("Node Weight = " + str(graph.node_weight(node)))
  # Iterate over all neighbors of the node
  print("Neighbors:")
  for edge in graph.incident_edges(node):
    print(graph.edge_target(edge), end = " ")
  print()
print()

# Iterate over the edges
# Note that the graph has six undirected edges, but internally we represent
# the graph as a directed graph. Thus, we iterate over 12 directed edges here.
# For an directed edge (u,v), we call u the source and v the target node.
for edge in graph.edges():
  print("Edge " + str(edge))
  print("Edge Weight = " + str(graph.edge_weight(edge)))
  print("Source Node = " + str(graph.edge_source(edge)))
  print("Target Node = " + str(graph.edge_target(edge)))
  print()
print()

# Create a partition of the graph
number_of_blocks = 3
partitioned_graph = graph.create_partitioned_hypergraph(
  context,
  number_of_blocks, # partition graph into three blocks
  [0,0,1,1,2])      # block IDs for each node

# Iterate over all nodes
for node in graph.nodes():
  print("Node " + str(node))
  print("Block ID                     = " + str(partitioned_graph.block_id(node)))
  print("Number of Incident Cut Edges = " + str(partitioned_graph.num_incident_cut_edges(node)))
print()

# Iterate over all edges
for edge in graph.edges():
  print("Edge " + str(edge))
  print("Connectivity = " + str(partitioned_graph.connectivity(edge))) # number of blocks contained in edge
