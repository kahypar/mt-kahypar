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
# In the following, we partition a hypergraph into four blocks
# and optimize the connectivity metric
context.setPartitioningParameters(
  4,                 # number of blocks
  0.03,              # imbalance parameter not relevant for partitioning with individual block weights
  mtkahypar.Objective.KM1) # objective function
mtkahypar.setSeed(42)      # seed
# Set individual target block weights for each block
context.max_block_weights = [
  2131,  # The weight of the first block must be smaller or equal than 2131
  1213,  # The weight of the second block must be smaller or equal than 1213
  7287,  # The weight of the third block must be smaller or equal than 7287
  2501]  # The weight of the fourth block must be smaller or equal than 2501
context.logging = True

# Load hypergraph from file
hypergraph = mtkahypar.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Output metrics
print("Partition Stats:")
print("km1       = " + str(partitioned_hg.km1()))
print("Weight of Block 0 = " + str(partitioned_hg.blockWeight(0)) + " (<= 2131)")
print("Weight of Block 1 = " + str(partitioned_hg.blockWeight(1)) + " (<= 1213)")
print("Weight of Block 2 = " + str(partitioned_hg.blockWeight(2)) + " (<= 7287)")
print("Weight of Block 3 = " + str(partitioned_hg.blockWeight(3)) + " (<= 2501)")
