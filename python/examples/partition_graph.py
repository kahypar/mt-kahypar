# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.

import os
import multiprocessing
import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize thread pool
mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtkahypar.Context()
context.loadPreset(mtkahypar.PresetType.DEFAULT) # corresponds to Mt-KaHyPar-D
# In the following, we partition a graph into two blocks
# with an allowed imbalance of 3% and optimize the cut metric
context.setPartitioningParameters(
  2,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.CUT) # objective function
mtkahypar.setSeed(42)      # seed
context.logging = True

# Load graph from file
graph = mtkahypar.Graph(
  mydir + "/../tests/test_instances/delaunay_n15.graph", # graph file
  mtkahypar.FileFormat.METIS) # graph is stored in Metis file format

# Partition graph
partitioned_graph = graph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_graph.imbalance()))
print("cut       = " + str(partitioned_graph.cut()))
print("Block Weights:")
print("Weight of Block 0 = " + str(partitioned_graph.blockWeight(0)))
print("Weight of Block 1 = " + str(partitioned_graph.blockWeight(1)))
