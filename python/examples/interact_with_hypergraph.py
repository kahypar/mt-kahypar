# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.

import mtkahypar

# Creates a weighted hypergraph
hypergraph = mtkahypar.Hypergraph(
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
print("Number of Nodes      = " + str(hypergraph.numNodes()))
print("Number of Hyperedges = " + str(hypergraph.numEdges()))
print("Number of Pins       = " + str(hypergraph.numPins()))
print("Weight of Hypergraph = " + str(hypergraph.totalWeight()))
print()

# Iterate over the nodes
hypergraph.doForAllNodes(lambda node : (
  print("Node " + str(node)),
  print("Node Degree = " + str(hypergraph.nodeDegree(node))),
  print("Node Weight = " + str(hypergraph.nodeWeight(node))),
  # Iterate over all incident hyperedges of the node
  print("Incident Hyperedges:"),
  hypergraph.doForAllIncidentEdges(node, lambda edge : print(edge, end = " ")),
  print()
))
print()

# Iterate over the hyperedges
hypergraph.doForAllEdges(lambda edge : (
  print("Hyperedge " + str(edge)),
  print("Edge Size   = " + str(hypergraph.edgeSize(edge))),
  print("Edge Weight = " + str(hypergraph.edgeWeight(edge))),
  # Iterate over all pins of the edge
  print("Pins:"),
  hypergraph.doForAllPins(edge, lambda pin : print(pin, end = " ")),
  print()
))
print()

# Create a partition of the hypergraph
number_of_blocks = 3
partitioned_hg = mtkahypar.PartitionedHypergraph(
  hypergraph,
  number_of_blocks, # partition hypergraph into three blocks
  [0,0,0,1,1,2])    # block IDs for each node

# Output Metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance())) # deviation from average block weight in percentage
print("cut       = " + str(partitioned_hg.cut())) # cut-net metric
print("km1       = " + str(partitioned_hg.km1())) # connectivity metric
print("Block Weights:")
print("Weight of Block 0 = " + str(partitioned_hg.blockWeight(0)))
print("Weight of Block 1 = " + str(partitioned_hg.blockWeight(1)))
print("Weight of Block 2 = " + str(partitioned_hg.blockWeight(2)))
print()

# Iterate over all nodes
hypergraph.doForAllNodes(lambda node : (
  print("Node " + str(node)),
  print("Block ID                     = " + str(partitioned_hg.blockID(node))),
  print("Number of Incident Cut Edges = " + str(partitioned_hg.numIncidentCutEdges(node)))
))
print()

# Iterate over all hyperedges
hypergraph.doForAllEdges(lambda edge : (
  print("Hyperedge " + str(edge)),
  print("Connectivity = " + str(partitioned_hg.connectivity(edge))), # number of blocks contained in edge
  print("Edge Size    = " + str(hypergraph.edgeSize(edge))),
  partitioned_hg.doForAllBlocksInEdge(edge, lambda block : (
    print("Number of pins contained in block " + str(block) + " = "
      + str(partitioned_hg.numPinsInBlock(edge, block)))
  ))
))