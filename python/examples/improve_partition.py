# Please follow the instructions in the README to install the python library interface and
# and copy mtkahyparmtkahypar.so and mtkahypargp.so to this folder to run the examples.

import os
import multiprocessing
import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize thread pool
mtkahypar.initializeThreadPool(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtkahypar.Context()
# In the following, we partition a hypergraph into eight blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(
  8,                  # number of blocks
  0.03,               # imbalance parameter
  mtkahypar.Objective.KM1)  # objective function
mtkahypar.setSeed(42) # seed
context.enableLogging(True)

# Load hypergraph from file
hypergraph = mtkahypar.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Read partition from file
partitioned_hg = mtkahypar.PartitionedHypergraph(
  hypergraph,
  8, # partition of the hypergraph contains eight blocks
  mydir + "/../tests/test_instances/ibm01.hgr.part8") # partition file

# Improve partition
km1_before = partitioned_hg.km1()
context.loadPreset(mtkahypar.PresetType.QUALITY) # use high quality preset for improvement
# We perform one multilevel improvement cycle (also called V-cycle)
mtkahypar.improvePartition(partitioned_hg, context, 1)
km1_after = partitioned_hg.km1()

# Output metrics
print("Partition Stats:")
print("Km1 before Improvement = " + str(km1_before))
print("Km1 after Improvement = " + str(km1_after))
