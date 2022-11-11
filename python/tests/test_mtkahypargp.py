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

import mtkahypargp as gp

mydir = os.path.dirname(os.path.realpath(__file__))
logging = False

class MainTest(unittest.TestCase):

  def setUp(self):
    gp.initializeThreadPool(multiprocessing.cpu_count())

  def test_check_graph_stats(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.numNodes(), 5)
    self.assertEqual(graph.numEdges(), 6)
    self.assertEqual(graph.numDirectedEdges(), 12)
    self.assertEqual(graph.totalWeight(), 5)

  def test_check_node_degrees(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.nodeDegree(0), 2)
    self.assertEqual(graph.nodeDegree(1), 3)
    self.assertEqual(graph.nodeDegree(2), 3)
    self.assertEqual(graph.nodeDegree(3), 3)
    self.assertEqual(graph.nodeDegree(4), 1)

  def test_check_node_weights(self):
    graph = gp.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,1,1,1,1,1])

    self.assertEqual(graph.totalWeight(), 15)
    self.assertEqual(graph.nodeWeight(0), 1)
    self.assertEqual(graph.nodeWeight(1), 2)
    self.assertEqual(graph.nodeWeight(2), 3)
    self.assertEqual(graph.nodeWeight(3), 4)
    self.assertEqual(graph.nodeWeight(4), 5)

  def test_check_edge_weights(self):
    graph = gp.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.edgeWeight(0),  1) # (0,1)
    self.assertEqual(graph.edgeWeight(1),  2) # (0,2)
    self.assertEqual(graph.edgeWeight(2),  1) # (1,0)
    self.assertEqual(graph.edgeWeight(3),  3) # (1,2)
    self.assertEqual(graph.edgeWeight(4),  4) # (1,3)
    self.assertEqual(graph.edgeWeight(5),  2) # (2,0)
    self.assertEqual(graph.edgeWeight(6),  3) # (2,1)
    self.assertEqual(graph.edgeWeight(7),  5) # (2,3)
    self.assertEqual(graph.edgeWeight(8),  4) # (3,1)
    self.assertEqual(graph.edgeWeight(9),  5) # (3,2)
    self.assertEqual(graph.edgeWeight(10), 6) # (3,4)
    self.assertEqual(graph.edgeWeight(11), 6) # (4,3)

  def test_source_nodes_of_edges(self):
    graph = gp.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.source(0),  0) # (0,1)
    self.assertEqual(graph.source(1),  0) # (0,2)
    self.assertEqual(graph.source(2),  1) # (1,0)
    self.assertEqual(graph.source(3),  1) # (1,2)
    self.assertEqual(graph.source(4),  1) # (1,3)
    self.assertEqual(graph.source(5),  2) # (2,0)
    self.assertEqual(graph.source(6),  2) # (2,1)
    self.assertEqual(graph.source(7),  2) # (2,3)
    self.assertEqual(graph.source(8),  3) # (3,1)
    self.assertEqual(graph.source(9),  3) # (3,2)
    self.assertEqual(graph.source(10), 3) # (3,4)
    self.assertEqual(graph.source(11), 4) # (4,3)

  def test_target_nodes_of_edges(self):
    graph = gp.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.target(0),  1) # (0,1)
    self.assertEqual(graph.target(1),  2) # (0,2)
    self.assertEqual(graph.target(2),  0) # (1,0)
    self.assertEqual(graph.target(3),  2) # (1,2)
    self.assertEqual(graph.target(4),  3) # (1,3)
    self.assertEqual(graph.target(5),  0) # (2,0)
    self.assertEqual(graph.target(6),  1) # (2,1)
    self.assertEqual(graph.target(7),  3) # (2,3)
    self.assertEqual(graph.target(8),  1) # (3,1)
    self.assertEqual(graph.target(9),  2) # (3,2)
    self.assertEqual(graph.target(10), 4) # (3,4)
    self.assertEqual(graph.target(11), 3) # (4,3)


  def test_load_graph_in_metis_file_format(self):
    graph = gp.Graph(
      mydir + "/test_instances/delaunay_n15.graph", gp.FileFormat.METIS)

    self.assertEqual(graph.numNodes(), 32768)
    self.assertEqual(graph.numEdges(), 98274)
    self.assertEqual(graph.numDirectedEdges(), 2 * 98274)
    self.assertEqual(graph.totalWeight(), 32768)

  def test_all_nodes_in_correct_block(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.blockID(0), 0)
    self.assertEqual(partitioned_graph.blockID(1), 1)
    self.assertEqual(partitioned_graph.blockID(2), 1)
    self.assertEqual(partitioned_graph.blockID(3), 2)
    self.assertEqual(partitioned_graph.blockID(4), 2)

  def test_correct_block_weights(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.blockWeight(0), 1)
    self.assertEqual(partitioned_graph.blockWeight(1), 2)
    self.assertEqual(partitioned_graph.blockWeight(2), 2)

  def test_cut_metric(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.cut(), 4)

  def test_all_nodes_contains_correct_number_of_incident_cut_edges(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.numIncidentCutEdges(0), 2)
    self.assertEqual(partitioned_graph.numIncidentCutEdges(1), 2)
    self.assertEqual(partitioned_graph.numIncidentCutEdges(2), 2)
    self.assertEqual(partitioned_graph.numIncidentCutEdges(3), 2)
    self.assertEqual(partitioned_graph.numIncidentCutEdges(4), 0)

  def test_all_edges_have_correct_connectivity(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.connectivity(0),  2) # (0,1)
    self.assertEqual(partitioned_graph.connectivity(1),  2) # (0,2)
    self.assertEqual(partitioned_graph.connectivity(2),  2) # (1,0)
    self.assertEqual(partitioned_graph.connectivity(3),  1) # (1,2)
    self.assertEqual(partitioned_graph.connectivity(4),  2) # (1,3)
    self.assertEqual(partitioned_graph.connectivity(5),  2) # (2,0)
    self.assertEqual(partitioned_graph.connectivity(6),  1) # (2,1)
    self.assertEqual(partitioned_graph.connectivity(7),  2) # (2,3)
    self.assertEqual(partitioned_graph.connectivity(8),  2) # (3,1)
    self.assertEqual(partitioned_graph.connectivity(9),  2) # (3,2)
    self.assertEqual(partitioned_graph.connectivity(10), 1) # (3,4)
    self.assertEqual(partitioned_graph.connectivity(11), 1) # (4,3)

  def test_load_partition_from_file(self):
    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3,
      mydir + "/test_instances/test_graph_partition.part3")

    self.assertEqual(partitioned_graph.blockID(0), 0)
    self.assertEqual(partitioned_graph.blockID(1), 0)
    self.assertEqual(partitioned_graph.blockID(2), 1)
    self.assertEqual(partitioned_graph.blockID(3), 1)
    self.assertEqual(partitioned_graph.blockID(4), 2)

  def test_write_partition_to_file(self):
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    graph = gp.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = gp.PartitionedGraph(graph, 3, [0,0,1,2,2])

    partitioned_graph.writePartitionToFile(mydir + "/test_partition.part3")
    partitioned_graph_2 = gp.PartitionedGraph(graph, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_graph_2.blockID(0), 0)
    self.assertEqual(partitioned_graph_2.blockID(1), 0)
    self.assertEqual(partitioned_graph_2.blockID(2), 1)
    self.assertEqual(partitioned_graph_2.blockID(3), 2)
    self.assertEqual(partitioned_graph_2.blockID(4), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  class GraphPartitioner(unittest.TestCase):

    def __init__(self, preset_type, num_blocks, epsilon, objective, force_logging):
      self.context = gp.Context()
      self.context.loadPreset(preset_type)
      self.context.setPartitioningParameters(num_blocks, epsilon, objective, 42)
      self.context.enableLogging(logging or force_logging)
      self.graph = gp.Graph(
        mydir + "/test_instances/delaunay_n15.graph", gp.FileFormat.METIS)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.graph.totalWeight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.setIndividualBlockWeights(individualBlockWeights)

    def partition(self):
      self.partitioned_graph = gp.partition(self.graph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      objective_before = self.partitioned_graph.cut()
      gp.improvePartition(self.partitioned_graph, self.context, num_vcycles)
      objective_after = self.partitioned_graph.cut()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def __verifyPartition(self):
      if not self.useIndividualBlockWeights:
        # Check if imbalance is smaller than allowed imbalance
        self.assertLessEqual(self.partitioned_graph.imbalance(), self.epsilon)
        # Check if block weights are smaller than maximum allowed block weight
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_graph.blockWeight(block), self.maxAllowedBlockWeight)
      else:
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_graph.blockWeight(block), self.individualBlockWeights[block])

      # Verify block IDs of nodes
      self.graph.doForAllNodes(lambda hn : (
        self.assertGreaterEqual(self.partitioned_graph.blockID(hn), 0),
        self.assertLess(self.partitioned_graph.blockID(hn), self.k)
      ))

  def test_partitions_a_graph_with_speed_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 2, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_graph_with_speed_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_graph_with_high_quality_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.HIGH_QUALITY, 2, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_graph_with_high_quality_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.HIGH_QUALITY, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.DETERMINISTIC, 2, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(gp.PresetType.DETERMINISTIC, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()

  def test_checks_if_deterministic_preset_produces_same_result_for_graph(self):
    partitioner = self.GraphPartitioner(gp.PresetType.DETERMINISTIC, 8, 0.03, gp.Objective.KM1, False)
    partitioner.partition()
    objective_1 = partitioner.partitioned_graph.cut()
    partitioner.partition()
    objective_2 = partitioner.partitioned_graph.cut()
    partitioner.partition()
    objective_3 = partitioner.partitioned_graph.cut()
    self.assertEqual(objective_1, objective_2)
    self.assertEqual(objective_1, objective_3)

  def test_improves_a_partition_with_one_vcycle(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_improves_a_partition_with_one_vcycle_and_different_preset_type(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()
    partitioner.context.loadPreset(gp.PresetType.HIGH_QUALITY)
    partitioner.improvePartition(1)

  def test_improves_a_partition_with_three_vcycle(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(3)

  def test_partitions_a_graph_with_individual_block_weights(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([11201,4384,14174,3989])
    partitioner.partition()

  def test_partitions_a_graph_with_individual_block_weights_and_one_vcycle(self):
    partitioner = self.GraphPartitioner(gp.PresetType.SPEED, 4, 0.03, gp.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([11201,4384,14174,3989])
    partitioner.partition()
    partitioner.improvePartition(1)

if __name__ == '__main__':
  unittest.main()
