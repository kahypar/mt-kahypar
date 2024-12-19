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
# In the following, we partition a hypergraph into four blocks
# and optimize the connectivity metric
context.set_partitioning_parameters(
  4,                 # number of blocks
  0.03,              # imbalance parameter not relevant for partitioning with individual block weights
  mtkahypar.Objective.KM1) # objective function
mtkahypar.set_seed(42)      # seed
# Set individual target block weights for each block
context.set_individual_target_block_weights([
  2131,  # The weight of the first block must be smaller or equal than 2131
  1213,  # The weight of the second block must be smaller or equal than 1213
  7287,  # The weight of the third block must be smaller or equal than 7287
  2501]) # The weight of the fourth block must be smaller or equal than 2501
context.logging = True

# Load hypergraph from file
hypergraph = mtk.hypergraph_from_file(
  mydir + "/../tests/test_instances/ibm01.hgr", # hypergraph file
  context,
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Get max block weights from context
max_block_weights = context.compute_max_block_weights(hypergraph.total_weight())

# Output metrics
print("Partition Stats:")
print("km1       = " + str(partitioned_hg.km1()))
print("Block Weights:")
for i in partitioned_hg.blocks():
  print(f"Weight of Block {i} = {partitioned_hg.block_weight(i)} (<= {max_block_weights[i]})")
