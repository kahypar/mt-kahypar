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
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(
  4,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.KM1) # objective function
mtkahypar.setSeed(42)      # seed
context.logging = True

# Load hypergraph from file
hypergraph = mtkahypar.Hypergraph(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Add fixed vertices
hypergraph.addFixedVerticesFromFile(
  mydir + "/../tests/test_instances/ibm01.k4.p1.fix", 4)

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance()))
print("km1       = " + str(partitioned_hg.km1()))
print("cut       = " + str(partitioned_hg.cut()))

correct_assignment = True
for hn in range(0, hypergraph.numNodes()):
  if partitioned_hg.isFixed(hn):
    if partitioned_hg.fixedVertexBlock(hn) != partitioned_hg.blockID(hn):
      print("Node" + str(hn) + "is fixed to block" + str(partitioned_hg.fixedVertexBlock(hn)) +
            ", but is assigned to block" + str(partitioned_hg.blockID(hn)))
      correct_assignment = False

if correct_assignment:
  print("\033[1;92mFixed vertex assignment was successful :)\033[0m")

