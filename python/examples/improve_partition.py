# Please follow the instructions in the README to install the python library interface and
# and copy mtkahyparmtkahypar.so and mtkahypargp.so to this folder to run the examples.

import os
import multiprocessing
import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize
mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
# In the following, we partition a hypergraph into eight blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.set_partitioning_parameters(
  8,                        # number of blocks
  0.03,                     # imbalance parameter
  mtkahypar.Objective.KM1)  # objective function
mtkahypar.set_seed(42)      # seed
context.logging = True

# Load hypergraph from file
hypergraph = mtk.hypergraph_from_file(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  context,
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Read partition from file
partitioned_hg = hypergraph.partitioned_hypergraph_from_file(
  context,
  8, # partition of the hypergraph contains eight blocks
  mydir + "/../tests/test_instances/ibm01.hgr.part8") # partition file

# Improve partition
km1_before = partitioned_hg.km1()
context = mtk.context_from_preset(mtkahypar.PresetType.QUALITY) # use high quality preset for improvement
context.set_partitioning_parameters(8, 0.03, mtkahypar.Objective.KM1)
context.logging = True

# We perform one multilevel improvement cycle (also called V-cycle)
partitioned_hg.improve_partition(context, 1)
km1_after = partitioned_hg.km1()

# Output metrics
print("Partition Stats:")
print("Km1 before Improvement = " + str(km1_before))
print("Km1 after Improvement = " + str(km1_after))
