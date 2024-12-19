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
# In the following, we partition a hypergraph into two blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.set_mapping_parameters(
  8,                       # number of blocks - number of nodes of the target graph
  0.03)                    # imbalance parameter
mtkahypar.set_seed(42)     # seed
context.logging = True

# Load hypergraph from file
hypergraph =  mtk.hypergraph_from_file(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  context,
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Load target graph from file
graph = mtk.target_graph_from_file(
  mydir + "/../tests/test_instances/target.graph", # graph file
  context,
  mtkahypar.FileFormat.METIS) # target graph is stored in Metis file format

# Map hypergraph onto graph (optimizes Steiner tree metric)
partitioned_hg = hypergraph.map_onto_graph(graph, context)

# Output metrics
print("Partition Stats:")
print("Imbalance    = " + str(partitioned_hg.imbalance(context)))
print("steiner_tree = " + str(partitioned_hg.steiner_tree(graph)))
print("km1          = " + str(partitioned_hg.km1()))
print("soed         = " + str(partitioned_hg.soed()))
print("cut          = " + str(partitioned_hg.cut()))
print("Block Weights:")
for i in partitioned_hg.blocks():
  print("Weight of Block " + str(i) + " = " + str(partitioned_hg.block_weight(i)))
