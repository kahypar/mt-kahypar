#******************************************************************************
# MIT License
#
# This file is part of Mt-KaHyPar.
#
# Copyright (C) 2022 Tobias Heuer <tobias.heuer@kit.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#*****************************************************************************/

import unittest
import os
import multiprocessing
import math

import mtkahyparhgp as hgp

mydir = os.path.dirname(os.path.realpath(__file__))
logging = False

class MainTest(unittest.TestCase):

  def setUp(self):
    hgp.initializeThreadPool(multiprocessing.cpu_count())

  def test_check_hypergraph_stats(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.numNodes(), 7)
    self.assertEqual(hypergraph.numEdges(), 4)
    self.assertEqual(hypergraph.numPins(), 12)
    self.assertEqual(hypergraph.totalWeight(), 7)

  def test_check_node_degrees(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.nodeDegree(0), 2)
    self.assertEqual(hypergraph.nodeDegree(1), 1)
    self.assertEqual(hypergraph.nodeDegree(2), 2)
    self.assertEqual(hypergraph.nodeDegree(3), 2)
    self.assertEqual(hypergraph.nodeDegree(4), 2)
    self.assertEqual(hypergraph.nodeDegree(5), 1)
    self.assertEqual(hypergraph.nodeDegree(6), 2)

  def test_check_edge_sizes(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.edgeSize(0), 2)
    self.assertEqual(hypergraph.edgeSize(1), 4)
    self.assertEqual(hypergraph.edgeSize(2), 3)
    self.assertEqual(hypergraph.edgeSize(3), 3)

  def test_check_node_weights(self):
    hypergraph = hgp.Hypergraph(
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,1,1,1])

    self.assertEqual(hypergraph.totalWeight(), 28)
    self.assertEqual(hypergraph.nodeWeight(0), 1)
    self.assertEqual(hypergraph.nodeWeight(1), 2)
    self.assertEqual(hypergraph.nodeWeight(2), 3)
    self.assertEqual(hypergraph.nodeWeight(3), 4)
    self.assertEqual(hypergraph.nodeWeight(4), 5)
    self.assertEqual(hypergraph.nodeWeight(5), 6)
    self.assertEqual(hypergraph.nodeWeight(6), 7)

  def test_check_edge_weights(self):
    hypergraph = hgp.Hypergraph(
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,2,3,4])

    self.assertEqual(hypergraph.edgeWeight(0), 1)
    self.assertEqual(hypergraph.edgeWeight(1), 2)
    self.assertEqual(hypergraph.edgeWeight(2), 3)
    self.assertEqual(hypergraph.edgeWeight(3), 4)

  def test_load_hypergraph_in_hmetis_file_format(self):
    hypergraph = hgp.Hypergraph(
      mydir + "/test_instances/ibm01.hgr", hgp.FileFormat.HMETIS)

    self.assertEqual(hypergraph.numNodes(), 12752)
    self.assertEqual(hypergraph.numEdges(), 14111)
    self.assertEqual(hypergraph.numPins(), 50566)
    self.assertEqual(hypergraph.totalWeight(), 12752)

  def test_load_hypergraph_in_metis_file_format(self):
    hypergraph = hgp.Hypergraph(
      mydir + "/test_instances/delaunay_n15.graph", hgp.FileFormat.METIS)

    self.assertEqual(hypergraph.numNodes(), 32768)
    self.assertEqual(hypergraph.numEdges(), 98274)
    self.assertEqual(hypergraph.totalWeight(), 32768)

  def test_all_nodes_in_correct_block(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.blockID(0), 0)
    self.assertEqual(partitioned_hg.blockID(1), 0)
    self.assertEqual(partitioned_hg.blockID(2), 0)
    self.assertEqual(partitioned_hg.blockID(3), 1)
    self.assertEqual(partitioned_hg.blockID(4), 1)
    self.assertEqual(partitioned_hg.blockID(5), 1)
    self.assertEqual(partitioned_hg.blockID(6), 2)

  def test_correct_block_weights(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.blockWeight(0), 3)
    self.assertEqual(partitioned_hg.blockWeight(1), 3)
    self.assertEqual(partitioned_hg.blockWeight(2), 1)

  def test_metrics(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.cut(), 3)
    self.assertEqual(partitioned_hg.km1(), 4)
    self.assertEqual(partitioned_hg.soed(), 7)

  def test_all_nodes_contains_correct_number_of_incident_cut_hyperedges(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.numIncidentCutEdges(0), 1)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(1), 1)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(2), 1)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(3), 2)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(4), 2)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(5), 1)
    self.assertEqual(partitioned_hg.numIncidentCutEdges(6), 2)

  def test_all_edges_have_correct_connectivity(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.connectivity(0), 1)
    self.assertEqual(partitioned_hg.connectivity(1), 2)
    self.assertEqual(partitioned_hg.connectivity(2), 2)
    self.assertEqual(partitioned_hg.connectivity(3), 3)

  def test_all_edges_have_correct_number_of_pins_in_blocks(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.numPinsInBlock(0,0), 2)
    self.assertEqual(partitioned_hg.numPinsInBlock(0,1), 0)
    self.assertEqual(partitioned_hg.numPinsInBlock(0,2), 0)
    self.assertEqual(partitioned_hg.numPinsInBlock(1,0), 2)
    self.assertEqual(partitioned_hg.numPinsInBlock(1,1), 2)
    self.assertEqual(partitioned_hg.numPinsInBlock(1,2), 0)
    self.assertEqual(partitioned_hg.numPinsInBlock(2,0), 0)
    self.assertEqual(partitioned_hg.numPinsInBlock(2,1), 2)
    self.assertEqual(partitioned_hg.numPinsInBlock(2,2), 1)
    self.assertEqual(partitioned_hg.numPinsInBlock(3,0), 1)
    self.assertEqual(partitioned_hg.numPinsInBlock(3,1), 1)
    self.assertEqual(partitioned_hg.numPinsInBlock(3,2), 1)

  def test_load_partition_from_file(self):
    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3,
      mydir + "/test_instances/test_partition.part3")

    self.assertEqual(partitioned_hg.blockID(0), 0)
    self.assertEqual(partitioned_hg.blockID(1), 0)
    self.assertEqual(partitioned_hg.blockID(2), 0)
    self.assertEqual(partitioned_hg.blockID(3), 1)
    self.assertEqual(partitioned_hg.blockID(4), 1)
    self.assertEqual(partitioned_hg.blockID(5), 1)
    self.assertEqual(partitioned_hg.blockID(6), 2)

  def test_write_partition_to_file(self):
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    hypergraph = hgp.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hgp.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,2,2])

    partitioned_hg.writePartitionToFile(mydir + "/test_partition.part3")
    partitioned_hg_2 = hgp.PartitionedHypergraph(hypergraph, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_hg.blockID(0), 0)
    self.assertEqual(partitioned_hg.blockID(1), 0)
    self.assertEqual(partitioned_hg.blockID(2), 0)
    self.assertEqual(partitioned_hg.blockID(3), 1)
    self.assertEqual(partitioned_hg.blockID(4), 1)
    self.assertEqual(partitioned_hg.blockID(5), 2)
    self.assertEqual(partitioned_hg.blockID(6), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  class HypergraphPartitioner(unittest.TestCase):

    def __init__(self, preset_type, num_blocks, epsilon, objective, force_logging):
      self.context = hgp.Context()
      self.context.loadPreset(preset_type)
      self.context.setPartitioningParameters(num_blocks, epsilon, objective, 42)
      self.context.enableLogging(logging or force_logging)
      self.hypergraph = hgp.Hypergraph(
        mydir + "/test_instances/ibm01.hgr", hgp.FileFormat.HMETIS)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.hypergraph.totalWeight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.setIndividualBlockWeights(individualBlockWeights)

    def partition(self):
      self.partitioned_hg = hgp.partition(self.hypergraph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      objective_before = self.partitioned_hg.km1()
      hgp.improvePartition(self.partitioned_hg, self.context, num_vcycles)
      objective_after = self.partitioned_hg.km1()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def __verifyPartition(self):
      if not self.useIndividualBlockWeights:
        # Check if imbalance is smaller than allowed imbalance
        self.assertLessEqual(self.partitioned_hg.imbalance(), self.epsilon)
        # Check if block weights are smaller than maximum allowed block weight
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_hg.blockWeight(block), self.maxAllowedBlockWeight)
      else:
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_hg.blockWeight(block), self.individualBlockWeights[block])

      # Verify block IDs of nodes
      self.hypergraph.doForAllNodes(lambda hn : (
        self.assertGreaterEqual(self.partitioned_hg.blockID(hn), 0),
        self.assertLess(self.partitioned_hg.blockID(hn), self.k)
      ))


  def test_partitions_a_hypergraph_with_speed_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 2, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_speed_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_high_quality_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.HIGH_QUALITY, 2, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_high_quality_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.HIGH_QUALITY, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.DETERMINISTIC, 2, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.DETERMINISTIC, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()

  def test_checks_if_deterministic_preset_produces_same_result_for_hypergraphs(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.DETERMINISTIC, 8, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()
    objective_1 = partitioner.partitioned_hg.km1()
    partitioner.partition()
    objective_2 = partitioner.partitioned_hg.km1()
    partitioner.partition()
    objective_3 = partitioner.partitioned_hg.km1()
    self.assertEqual(objective_1, objective_2)
    self.assertEqual(objective_1, objective_3)

  def test_improves_a_partition_with_one_vcycle(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_improves_a_partition_with_one_vcycle_and_different_preset_type(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()
    partitioner.context.loadPreset(hgp.PresetType.HIGH_QUALITY)
    partitioner.improvePartition(1)

  def test_improves_a_partition_with_three_vcycle(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(3)

  def test_partitions_a_hypergraph_with_individual_block_weights(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([2131,1213,7287,2501])
    partitioner.partition()

  def test_partitions_a_hypergraph_with_individual_block_weights_and_one_vcycle(self):
    partitioner = self.HypergraphPartitioner(hgp.PresetType.SPEED, 4, 0.03, hgp.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([2131,1213,7287,2501])
    partitioner.partition()
    partitioner.improvePartition(1)

if __name__ == '__main__':
  unittest.main()
