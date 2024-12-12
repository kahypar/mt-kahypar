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

import mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))
logging = False

class MainTest(unittest.TestCase):

  def setUp(self):
    mtkahypar.initialize(multiprocessing.cpu_count())

  def test_set_partitioning_parameters_in_context(self):
    context = mtkahypar.Context(mtkahypar.PresetType.DEFAULT)
    context.set_partitioning_parameters(2, 0.03, mtkahypar.Objective.KM1)
    self.assertEqual(context.k, 2)
    self.assertEqual(context.epsilon, 0.03)
    self.assertEqual(context.objective, mtkahypar.Objective.KM1)

  def test_set_partitioning_parameters_over_properties(self):
    context = mtkahypar.Context(mtkahypar.PresetType.DEFAULT)
    context.k = 4
    context.epsilon = 0.05
    context.objective = mtkahypar.Objective.CUT
    context.num_vcycles = 5
    context.logging = True
    context.max_block_weights = [100, 200, 300, 400]

    self.assertEqual(context.k, 4)
    self.assertEqual(context.epsilon, 0.05)
    self.assertEqual(context.objective, mtkahypar.Objective.CUT)
    self.assertEqual(context.num_vcycles, 5)
    self.assertEqual(context.logging, True)
    self.assertEqual(context.max_block_weights[0], 100)
    self.assertEqual(context.max_block_weights[1], 200)
    self.assertEqual(context.max_block_weights[2], 300)
    self.assertEqual(context.max_block_weights[3], 400)

  def test_check_graph_stats(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.num_nodes(), 5)
    self.assertEqual(graph.num_edges(), 6)
    self.assertEqual(graph.num_directed_edges(), 12)
    self.assertEqual(graph.total_weight(), 5)

  def test_check_graph_node_degrees(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.node_degree(0), 2)
    self.assertEqual(graph.node_degree(1), 3)
    self.assertEqual(graph.node_degree(2), 3)
    self.assertEqual(graph.node_degree(3), 3)
    self.assertEqual(graph.node_degree(4), 1)

  def test_check_graph_node_weights(self):
    graph = mtkahypar.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,1,1,1,1,1])

    self.assertEqual(graph.total_weight(), 15)
    self.assertEqual(graph.node_weight(0), 1)
    self.assertEqual(graph.node_weight(1), 2)
    self.assertEqual(graph.node_weight(2), 3)
    self.assertEqual(graph.node_weight(3), 4)
    self.assertEqual(graph.node_weight(4), 5)

  def test_check_graph_edge_weights(self):
    graph = mtkahypar.Graph(
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.edge_weight(0),  1) # (0,1)
    self.assertEqual(graph.edge_weight(1),  2) # (0,2)
    self.assertEqual(graph.edge_weight(2),  1) # (1,0)
    self.assertEqual(graph.edge_weight(3),  3) # (1,2)
    self.assertEqual(graph.edge_weight(4),  4) # (1,3)
    self.assertEqual(graph.edge_weight(5),  2) # (2,0)
    self.assertEqual(graph.edge_weight(6),  3) # (2,1)
    self.assertEqual(graph.edge_weight(7),  5) # (2,3)
    self.assertEqual(graph.edge_weight(8),  4) # (3,1)
    self.assertEqual(graph.edge_weight(9),  5) # (3,2)
    self.assertEqual(graph.edge_weight(10), 6) # (3,4)
    self.assertEqual(graph.edge_weight(11), 6) # (4,3)

  def test_source_nodes_of_edges(self):
    graph = mtkahypar.Graph(
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
    graph = mtkahypar.Graph(
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
    graph = mtkahypar.Graph(
      mydir + "/test_instances/delaunay_n15.graph", mtkahypar.FileFormat.METIS)

    self.assertEqual(graph.num_nodes(), 32768)
    self.assertEqual(graph.num_edges(), 98274)
    self.assertEqual(graph.num_directed_edges(), 2 * 98274)
    self.assertEqual(graph.total_weight(), 32768)

  def test_check_hypergraph_stats(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.num_nodes(), 7)
    self.assertEqual(hypergraph.num_edges(), 4)
    self.assertEqual(hypergraph.num_pins(), 12)
    self.assertEqual(hypergraph.total_weight(), 7)

  def test_check_hypergraph_node_degrees(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.node_degree(0), 2)
    self.assertEqual(hypergraph.node_degree(1), 1)
    self.assertEqual(hypergraph.node_degree(2), 2)
    self.assertEqual(hypergraph.node_degree(3), 2)
    self.assertEqual(hypergraph.node_degree(4), 2)
    self.assertEqual(hypergraph.node_degree(5), 1)
    self.assertEqual(hypergraph.node_degree(6), 2)

  def test_check_hypergraph_edge_sizes(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.edge_size(0), 2)
    self.assertEqual(hypergraph.edge_size(1), 4)
    self.assertEqual(hypergraph.edge_size(2), 3)
    self.assertEqual(hypergraph.edge_size(3), 3)

  def test_check_hypergraph_node_weights(self):
    hypergraph = mtkahypar.Hypergraph(
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,1,1,1])

    self.assertEqual(hypergraph.total_weight(), 28)
    self.assertEqual(hypergraph.node_weight(0), 1)
    self.assertEqual(hypergraph.node_weight(1), 2)
    self.assertEqual(hypergraph.node_weight(2), 3)
    self.assertEqual(hypergraph.node_weight(3), 4)
    self.assertEqual(hypergraph.node_weight(4), 5)
    self.assertEqual(hypergraph.node_weight(5), 6)
    self.assertEqual(hypergraph.node_weight(6), 7)

  def test_check_hypergraph_edge_weights(self):
    hypergraph = mtkahypar.Hypergraph(
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,2,3,4])

    self.assertEqual(hypergraph.edge_weight(0), 1)
    self.assertEqual(hypergraph.edge_weight(1), 2)
    self.assertEqual(hypergraph.edge_weight(2), 3)
    self.assertEqual(hypergraph.edge_weight(3), 4)

  def test_load_hypergraph_in_hmetis_file_format(self):
    hypergraph = mtkahypar.Hypergraph(
      mydir + "/test_instances/ibm01.hgr", mtkahypar.FileFormat.HMETIS)

    self.assertEqual(hypergraph.num_nodes(), 12752)
    self.assertEqual(hypergraph.num_edges(), 14111)
    self.assertEqual(hypergraph.num_pins(), 50566)
    self.assertEqual(hypergraph.total_weight(), 12752)

  def test_load_hypergraph_in_metis_file_format(self):
    hypergraph = mtkahypar.Hypergraph(
      mydir + "/test_instances/delaunay_n15.graph", mtkahypar.FileFormat.METIS)

    self.assertEqual(hypergraph.num_nodes(), 32768)
    self.assertEqual(hypergraph.num_edges(), 98274)
    self.assertEqual(hypergraph.total_weight(), 32768)

  def test_for_graph_if_all_nodes_in_correct_block(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.block_id(0), 0)
    self.assertEqual(partitioned_graph.block_id(1), 1)
    self.assertEqual(partitioned_graph.block_id(2), 1)
    self.assertEqual(partitioned_graph.block_id(3), 2)
    self.assertEqual(partitioned_graph.block_id(4), 2)

  def test_for_graph_if_block_have_correct_weight(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.block_weight(0), 1)
    self.assertEqual(partitioned_graph.block_weight(1), 2)
    self.assertEqual(partitioned_graph.block_weight(2), 2)

  def test_cut_metric_for_graph(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.cut(), 4)

  def test_for_graph_if_all_nodes_contains_correct_number_of_incident_cut_edges(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.num_incident_cut_edges(0), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(1), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(2), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(4), 0)

  def test_for_graph_if_all_edges_have_correct_connectivity(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,1,1,2,2])

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

  def test_load_graph_partition_from_file(self):
    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3,
      mydir + "/test_instances/test_graph_partition.part3")

    self.assertEqual(partitioned_graph.block_id(0), 0)
    self.assertEqual(partitioned_graph.block_id(1), 0)
    self.assertEqual(partitioned_graph.block_id(2), 1)
    self.assertEqual(partitioned_graph.block_id(3), 1)
    self.assertEqual(partitioned_graph.block_id(4), 2)

  def test_write_graph_partition_to_file(self):
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    graph = mtkahypar.Graph(5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = mtkahypar.PartitionedGraph(graph, 3, [0,0,1,2,2])

    partitioned_graph.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_graph_2 = mtkahypar.PartitionedGraph(graph, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_graph_2.block_id(0), 0)
    self.assertEqual(partitioned_graph_2.block_id(1), 0)
    self.assertEqual(partitioned_graph_2.block_id(2), 1)
    self.assertEqual(partitioned_graph_2.block_id(3), 2)
    self.assertEqual(partitioned_graph_2.block_id(4), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  def test_for_hypergraph_if_all_nodes_are_in_correct_block(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_for_hypergraph_if_blocks_have_correct_weight(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_weight(0), 3)
    self.assertEqual(partitioned_hg.block_weight(1), 3)
    self.assertEqual(partitioned_hg.block_weight(2), 1)

  def test_all_metrics_for_hypergraph(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.cut(), 3)
    self.assertEqual(partitioned_hg.km1(), 4)

  def test_for_hypergraph_if_all_nodes_contains_correct_number_of_incident_cut_hyperedges(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.num_incident_cut_edges(0), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(1), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(2), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(4), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(5), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(6), 2)

  def test_for_hypergraph_if_all_edges_have_correct_connectivity(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.connectivity(0), 1)
    self.assertEqual(partitioned_hg.connectivity(1), 2)
    self.assertEqual(partitioned_hg.connectivity(2), 2)
    self.assertEqual(partitioned_hg.connectivity(3), 3)

  def test_for_hypergraph_if_all_edges_have_correct_number_of_pins_in_blocks(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.num_pins_in_block(0,0), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(0,1), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(0,2), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,0), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,1), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,2), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,0), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,1), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,2), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,0), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,1), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,2), 1)

  def test_load_hypergraph_partition_from_file(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3,
      mydir + "/test_instances/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_write_hypergraph_partition_to_file(self):
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.PartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,2,2])

    partitioned_hg.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_hg_2 = mtkahypar.PartitionedHypergraph(hypergraph, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 2)
    self.assertEqual(partitioned_hg.block_id(6), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  def test_for_sparse_hypergraph_if_all_nodes_are_in_correct_block(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_for_sparse_hypergraph_if_blocks_have_correct_weight(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_weight(0), 3)
    self.assertEqual(partitioned_hg.block_weight(1), 3)
    self.assertEqual(partitioned_hg.block_weight(2), 1)

  def test_all_metrics_for_sparse_hypergraph(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.cut(), 3)
    self.assertEqual(partitioned_hg.km1(), 4)

  def test_for_sparse_hypergraph_if_all_nodes_contains_correct_number_of_incident_cut_hyperedges(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.num_incident_cut_edges(0), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(1), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(2), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(4), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(5), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(6), 2)

  def test_for_sparse_hypergraph_if_all_edges_have_correct_connectivity(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.connectivity(0), 1)
    self.assertEqual(partitioned_hg.connectivity(1), 2)
    self.assertEqual(partitioned_hg.connectivity(2), 2)
    self.assertEqual(partitioned_hg.connectivity(3), 3)

  def test_for_sparse_hypergraph_if_all_edges_have_correct_number_of_pins_in_blocks(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.num_pins_in_block(0,0), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(0,1), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(0,2), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,0), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,1), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(1,2), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,0), 0)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,1), 2)
    self.assertEqual(partitioned_hg.num_pins_in_block(2,2), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,0), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,1), 1)
    self.assertEqual(partitioned_hg.num_pins_in_block(3,2), 1)

  def test_load_sparse_hypergraph_partition_from_file(self):
    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3,
      mydir + "/test_instances/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_write_sparse_hypergraph_partition_to_file(self):
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    hypergraph = mtkahypar.Hypergraph(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3, [0,0,0,1,1,2,2])

    partitioned_hg.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_hg_2 = mtkahypar.SparsePartitionedHypergraph(hypergraph, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 2)
    self.assertEqual(partitioned_hg.block_id(6), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  class GraphPartitioner(unittest.TestCase):

    def __init__(self, preset_type, num_blocks, epsilon, objective, force_logging):
      self.context = mtkahypar.Context(preset_type)
      self.context.set_partitioning_parameters(num_blocks, epsilon, objective)
      mtkahypar.set_seed(42)
      self.context.logging = logging or force_logging
      self.target_graph = mtkahypar.Graph(
        mydir + "/test_instances/target.graph",  mtkahypar.FileFormat.METIS)
      self.graph = mtkahypar.Graph(
        mydir + "/test_instances/delaunay_n15.graph", mtkahypar.FileFormat.METIS)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.graph.total_weight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.max_block_weights = individualBlockWeights

    def addFixedVertices(self):
      self.graph.add_fixed_vertices_from_file(mydir + "/test_instances/delaunay_n15.k4.p1.fix", self.k)

    def partition(self):
      self.partitioned_graph = self.graph.partition(self.context)
      self.__verifyPartition()

    def map_onto_graph(self):
      self.partitioned_graph = self.graph.map_onto_graph(self.target_graph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      objective_before = self.partitioned_graph.cut()
      self.partitioned_graph.improve_partition(self.context, num_vcycles)
      objective_after = self.partitioned_graph.cut()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def improveMapping(self, num_vcycles):
      objective_before = self.partitioned_graph.steiner_tree(self.target_graph)
      self.partitioned_graph.improve_mapping(self.target_graph, self.context, num_vcycles)
      objective_after = self.partitioned_graph.steiner_tree(self.target_graph)
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def __verifyPartition(self):
      if not self.useIndividualBlockWeights:
        # Check if imbalance is smaller than allowed imbalance
        self.assertLessEqual(self.partitioned_graph.imbalance(self.context), self.epsilon)
        # Check if block weights are smaller than maximum allowed block weight
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_graph.block_weight(block), self.maxAllowedBlockWeight)
      else:
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_graph.block_weight(block), self.individualBlockWeights[block])

      # Verify block IDs of nodes
      self.graph.do_for_all_nodes(lambda hn : (
        self.assertGreaterEqual(self.partitioned_graph.block_id(hn), 0),
        self.assertLess(self.partitioned_graph.block_id(hn), self.k)
      ))

      # Verify block IDs of fixed vertices
      for hn in range(0, self.graph.num_nodes()):
        if self.partitioned_graph.is_fixed(hn):
          if self.partitioned_graph.block_id(hn) != self.partitioned_graph.fixed_vertex_block(hn):
            print("Wrong fixed vertex assignment")
            self.assertTrue(False)

  def test_partitions_a_graph_with_default_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_default_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_quality_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.QUALITY, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_quality_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.QUALITY, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_into_a_large_number_of_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.LARGE_K, 1024, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_checks_if_deterministic_preset_produces_same_result_for_graph(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 8, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    objective_1 = partitioner.partitioned_graph.cut()
    partitioner.partition()
    objective_2 = partitioner.partitioned_graph.cut()
    partitioner.partition()
    objective_3 = partitioner.partitioned_graph.cut()
    self.assertEqual(objective_1, objective_2)
    self.assertEqual(objective_1, objective_3)

  def test_improves_a_graph_partition_with_one_vcycle(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_improves_a_graph_partition_with_one_vcycle_and_different_preset_type(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    partitioner.context = mtkahypar.Context(mtkahypar.PresetType.QUALITY)
    partitioner.context.set_partitioning_parameters(4, 0.03, mtkahypar.Objective.CUT)
    partitioner.improvePartition(1)

  def test_improves_a_graph_partition_with_three_vcycle(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    partitioner.improvePartition(3)

  def test_partitions_a_graph_with_individual_block_weights(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.setIndividualBlockWeights([11201,4384,14174,3989])
    partitioner.partition()

  def test_partitions_a_graph_with_individual_block_weights_and_one_vcycle(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.setIndividualBlockWeights([11201,4384,14174,3989])
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_maps_a_graph_with_default_preset_onto_target_graph(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 8, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.map_onto_graph()

  def test_maps_a_graph_with_quality_preset_onto_target_graph(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.QUALITY, 8, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.map_onto_graph()

  def test_improves_mapping_of_a_graph_with_default_preset(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 8, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.map_onto_graph()
    partitioner.improveMapping(1)

  def test_partitions_a_graph_with_fixed_vertices_and_default_preset(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.addFixedVertices()
    partitioner.partition()

  def test_partitions_a_graph_with_fixed_vertices_and_quality_preset(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.QUALITY, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.addFixedVertices()
    partitioner.partition()

  def test_improve_graph_partition_with_fixed_vertices(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.addFixedVertices()
    partitioner.partition()
    partitioner.improvePartition(1)

  class HypergraphPartitioner(unittest.TestCase):

    def __init__(self, preset_type, num_blocks, epsilon, objective, force_logging):
      self.context = mtkahypar.Context(preset_type)
      self.context.set_partitioning_parameters(num_blocks, epsilon, objective)
      mtkahypar.set_seed(42)
      self.context.logging = logging or force_logging
      self.target_graph = mtkahypar.Graph(
        mydir + "/test_instances/target.graph",  mtkahypar.FileFormat.METIS)
      self.hypergraph = mtkahypar.Hypergraph(
        mydir + "/test_instances/ibm01.hgr", mtkahypar.FileFormat.HMETIS)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.preset_type = preset_type
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.hypergraph.total_weight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.max_block_weights = individualBlockWeights

    def addFixedVertices(self):
      self.hypergraph.add_fixed_vertices_from_file(mydir + "/test_instances/ibm01.k4.p1.fix", self.k)

    def partition(self):
      if self.preset_type == mtkahypar.PresetType.LARGE_K:
        self.partitioned_hg = self.hypergraph.partition_into_large_k(self.context)
      else:
        self.partitioned_hg = self.hypergraph.partition(self.context)
      self.__verifyPartition()

    def map_onto_graph(self):
      self.partitioned_hg = self.hypergraph.map_onto_graph(self.target_graph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      objective_before = self.partitioned_hg.km1()
      self.partitioned_hg.improve_partition(self.context, num_vcycles)
      objective_after = self.partitioned_hg.km1()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def improveMapping(self, num_vcycles):
      objective_before = self.partitioned_hg.steiner_tree(self.target_graph)
      self.partitioned_hg.improve_mapping(self.target_graph, self.context, num_vcycles)
      objective_after = self.partitioned_hg.steiner_tree(self.target_graph)
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def __verifyPartition(self):
      if not self.useIndividualBlockWeights:
        # Check if imbalance is smaller than allowed imbalance
        self.assertLessEqual(self.partitioned_hg.imbalance(self.context), self.epsilon)
        # Check if block weights are smaller than maximum allowed block weight
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_hg.block_weight(block), self.maxAllowedBlockWeight)
      else:
        for block in range(self.k):
          self.assertLessEqual(self.partitioned_hg.block_weight(block), self.individualBlockWeights[block])

      # Verify block IDs of nodes
      self.hypergraph.do_for_all_nodes(lambda hn : (
        self.assertGreaterEqual(self.partitioned_hg.block_id(hn), 0),
        self.assertLess(self.partitioned_hg.block_id(hn), self.k)
      ))

      # Verify block IDs of fixed vertices
      for hn in range(0, self.hypergraph.num_nodes()):
        if self.partitioned_hg.is_fixed(hn):
          if self.partitioned_hg.block_id(hn) != self.partitioned_hg.fixed_vertex_block(hn):
            print("Wrong fixed vertex assignment")
            self.assertTrue(False)

  def test_partitions_a_hypergraph_with_default_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_default_preset_into_four_blocks_km1(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_default_preset_into_four_blocks_soed(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.SOED, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_quality_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.QUALITY, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_quality_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.QUALITY, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_into_a_large_number_of_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.LARGE_K, 512, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_checks_if_deterministic_preset_produces_same_result_for_hypergraphs(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 8, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()
    objective_1 = partitioner.partitioned_hg.km1()
    partitioner.partition()
    objective_2 = partitioner.partitioned_hg.km1()
    partitioner.partition()
    objective_3 = partitioner.partitioned_hg.km1()
    self.assertEqual(objective_1, objective_2)
    self.assertEqual(objective_1, objective_3)

  def test_improves_a_hypergraph_partition_with_one_vcycle(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_improves_a_hypergraph_partition_with_one_vcycle_and_different_preset_type(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()
    partitioner.context = mtkahypar.Context(mtkahypar.PresetType.QUALITY)
    partitioner.context.set_partitioning_parameters(4, 0.03, mtkahypar.Objective.KM1)
    partitioner.improvePartition(1)

  def test_improves_a_hypergraph_partition_with_three_vcycle(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()
    partitioner.improvePartition(3)

  def test_partitions_a_hypergraph_with_individual_block_weights(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([2131,1213,7287,2501])
    partitioner.partition()

  def test_partitions_a_hypergraph_with_individual_block_weights_and_one_vcycle(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.setIndividualBlockWeights([2131,1213,7287,2501])
    partitioner.partition()
    partitioner.improvePartition(1)

  def test_maps_a_hypergraph_with_default_preset_onto_target_graph(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 8, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.map_onto_graph()

  def test_maps_a_hypergraph_with_quality_preset_onto_target_graph(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.QUALITY, 8, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.map_onto_graph()

  def test_improves_mapping_of_a_hypergraph_with_default_preset(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 8, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.map_onto_graph()
    partitioner.improveMapping(1)

  def test_partitions_a_hypergraph_with_fixed_vertices_and_default_preset(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.addFixedVertices()
    partitioner.partition()

  def test_partitions_a_hypergraph_with_fixed_vertices_and_quality_preset(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.QUALITY, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.addFixedVertices()
    partitioner.partition()

  def test_improve_hypergraph_partition_with_fixed_vertices(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.addFixedVertices()
    partitioner.partition()
    partitioner.improvePartition(1)

if __name__ == '__main__':
  unittest.main()
