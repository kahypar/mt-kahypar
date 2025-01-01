# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.

import os
import multiprocessing
import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize
mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
# In the following, we partition a graph into two blocks
# with an allowed imbalance of 3% and optimize the cut metric
context.set_partitioning_parameters(
  2,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.CUT) # objective function
mtkahypar.set_seed(42)      # seed
context.logging = True

# Load graph from file
graph = mtk.graph_from_file(
  mydir + "/../tests/test_instances/delaunay_n15.graph", # graph file
  context,
  mtkahypar.FileFormat.METIS) # graph is stored in Metis file format

# Partition graph
partitioned_graph = graph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_graph.imbalance(context)))
print("cut       = " + str(partitioned_graph.cut()))
print("Block Weights:")
for i in partitioned_graph.blocks():
  print("Weight of Block " + str(i) + " = " + str(partitioned_graph.block_weight(i)))
