# Please follow the instructions in the README to install the python library interface
# and copy mtkahypar.so to the folder containing the examples.

import mtkahypar

# Creates a weighted hypergraph
graph = mtkahypar.Graph(
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
  [1,2,1,3,3,2]) # hyperedge weights

# Output statistics of graph
print("Graph Stats:")
print("Number of Nodes            = " + str(graph.numNodes()))
print("Number of Undirected Edges = " + str(graph.numEdges()))
print("Number of Directed Edges   = " + str(graph.numDirectedEdges()))
print("Weight of Graph            = " + str(graph.totalWeight()))
print()

# Iterate over the nodes
graph.doForAllNodes(lambda node : (
  print("Node " + str(node)),
  print("Node Degree = " + str(graph.nodeDegree(node))),
  print("Node Weight = " + str(graph.nodeWeight(node))),
  # Iterate over all neighbors of the node
  print("Neighbors:"),
  graph.doForAllNeighbors(node, lambda neighbor : print(neighbor, end = " ")),
  print()
))
print()

# Iterate over the edges
# Note that the graph has six undirected edges, but internally we represent
# the graph as a directed graph. Thus, we iterate over 12 directed edges here.
# For an directed edge (u,v), we call u the source and v the target node.
graph.doForAllEdges(lambda edge : (
  print("Edge " + str(edge)),
  print("Edge Weight = " + str(graph.edgeWeight(edge))),
  print("Source Node = " + str(graph.source(edge))),
  print("Target Node = " + str(graph.target(edge))),
  print()
))
print()

# Create a partition of the graph
number_of_blocks = 3
partitioned_graph = mtkahypar.PartitionedGraph(
  graph,
  number_of_blocks, # partition hypergraph into three blocks
  [0,0,1,1,2])      # block IDs for each node

# Output Metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_graph.imbalance())) # deviation from average block weight in percentage
print("cut       = " + str(partitioned_graph.cut())) # edge-cut metric
print("Block Weights:")
print("Weight of Block 0 = " + str(partitioned_graph.blockWeight(0)))
print("Weight of Block 1 = " + str(partitioned_graph.blockWeight(1)))
print("Weight of Block 2 = " + str(partitioned_graph.blockWeight(2)))
print()

# Iterate over all nodes
graph.doForAllNodes(lambda node : (
  print("Node " + str(node)),
  print("Block ID                     = " + str(partitioned_graph.blockID(node))),
  print("Number of Incident Cut Edges = " + str(partitioned_graph.numIncidentCutEdges(node)))
))
print()

# Iterate over all edges
graph.doForAllEdges(lambda edge : (
  print("Hyperedge " + str(edge)),
  print("Connectivity = " + str(partitioned_graph.connectivity(edge))), # number of blocks contained in edge
))