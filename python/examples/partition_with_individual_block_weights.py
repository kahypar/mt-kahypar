# Please follow the instructions in the README to install the python library interface and
# and copy mtkahyparhgp.so and mtkahypargp.so to this folder to run the examples.

import os
import multiprocessing
import mtkahyparhgp as hgp

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize thread pool
hgp.initializeThreadPool(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = hgp.Context()
context.loadPreset(hgp.PresetType.SPEED) # corresponds to Mt-KaHyPar-D
# In the following, we partition a hypergraph into four blocks
# and optimize the connectivity metric
context.setPartitioningParameters(
  4,                 # number of blocks
  0.03,              # imbalance parameter not relevant for partitioning with individual block weights
  hgp.Objective.KM1, # objective function
  42)                # seed
# Set individual target block weights for each block
context.setIndividualBlockWeights([
  2131,  # The weight of the first block must be smaller or equal than 2131
  1213,  # The weight of the second block must be smaller or equal than 1213
  7287,  # The weight of the third block must be smaller or equal than 7287
  2501]) # The weight of the fourth block must be smaller or equal than 2501
context.enableLogging(True)

# Load hypergraph from file
hypergraph = hgp.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  hgp.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hgp.partition(hypergraph, context)

# Output metrics
print("Partition Stats:")
print("km1       = " + str(partitioned_hg.km1()))
print("Weight of Block 0 = " + str(partitioned_hg.blockWeight(0)) + " (<= 2131)")
print("Weight of Block 1 = " + str(partitioned_hg.blockWeight(1)) + " (<= 1213)")
print("Weight of Block 2 = " + str(partitioned_hg.blockWeight(2)) + " (<= 7287)")
print("Weight of Block 3 = " + str(partitioned_hg.blockWeight(3)) + " (<= 2501)")
