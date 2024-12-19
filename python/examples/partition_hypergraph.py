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
context.set_partitioning_parameters(
  2,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.KM1) # objective function
mtkahypar.set_seed(42)      # seed
context.logging = True

# Load hypergraph from file
hypergraph = mtk.hypergraph_from_file(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  context,
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance(context)))
print("km1       = " + str(partitioned_hg.km1()))
print("cut       = " + str(partitioned_hg.cut()))
print("Block Weights:")
for i in partitioned_hg.blocks():
  print("Weight of Block " + str(i) + " = " + str(partitioned_hg.block_weight(i)))
