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
# In the following, we partition a hypergraph into two blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(
  2,                 # number of blocks
  0.03,              # imbalance parameter
  hgp.Objective.KM1, # objective function
  42)                # seed
context.enableLogging(True)

# Load hypergraph from file
hypergraph = hgp.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  hgp.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hgp.partition(hypergraph, context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance()))
print("km1       = " + str(partitioned_hg.km1()))
print("Block Weights:")
print("Weight of Block 0 = " + str(partitioned_hg.blockWeight(0)))
print("Weight of Block 1 = " + str(partitioned_hg.blockWeight(1)))
