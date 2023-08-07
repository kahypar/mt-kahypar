# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.

import os
import multiprocessing
import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize thread pool
mtkahypar.initializeThreadPool(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtkahypar.Context()
context.loadPreset(mtkahypar.PresetType.DEFAULT) # corresponds to Mt-KaHyPar-D
# In the following, we partition a hypergraph into two blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(
  8,                       # number of blocks - number of nodes of the target graph
  0.03,                    # imbalance parameter
  mtkahypar.Objective.KM1) # objective function - not relevant for mapping
mtkahypar.setSeed(42)      # seed
context.logging = True

# Load hypergraph from file
hypergraph = mtkahypar.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Load target graph from file
graph = mtkahypar.Graph(
  mydir + "/../tests/test_instances/target.graph", # graph file
  mtkahypar.FileFormat.METIS) # target graph is stored in Metis file format

# Map hypergraph onto graph (optimizes Steiner tree metric)
partitioned_hg = hypergraph.mapOntoGraph(graph, context)

# Output metrics
print("Partition Stats:")
print("Imbalance    = " + str(partitioned_hg.imbalance()))
print("steiner_tree = " + str(partitioned_hg.steiner_tree(graph)))
print("km1          = " + str(partitioned_hg.km1()))
print("soed         = " + str(partitioned_hg.soed()))
print("cut          = " + str(partitioned_hg.cut()))
print("Block Weights:")
for i in range(0,8):
  print("Weight of Block " + str(i) + " = " + str(partitioned_hg.blockWeight(i)))
