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

mtk = mtkahypar.initialize(multiprocessing.cpu_count())

class MainTest(unittest.TestCase):

  def test_set_partitioning_parameters_in_context(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    context.set_partitioning_parameters(2, 0.03, mtkahypar.Objective.KM1)
    self.assertEqual(context.k, 2)
    self.assertEqual(context.epsilon, 0.03)
    self.assertEqual(context.objective, mtkahypar.Objective.KM1)

  def test_set_partitioning_parameters_over_properties(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    context.k = 4
    context.epsilon = 0.05
    context.objective = mtkahypar.Objective.CUT
    context.num_vcycles = 5
    context.logging = True

    self.assertEqual(context.k, 4)
    self.assertEqual(context.epsilon, 0.05)
    self.assertEqual(context.objective, mtkahypar.Objective.CUT)
    self.assertEqual(context.num_vcycles, 5)
    self.assertEqual(context.logging, True)

  def test_get_and_set_max_block_weights(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    context.k = 4
    context.epsilon = 0.05
    perfect_balance_weights = context.compute_perfect_balance_block_weights(1000)
    self.assertEqual(perfect_balance_weights[0], 250)
    self.assertEqual(perfect_balance_weights[1], 250)
    self.assertEqual(perfect_balance_weights[2], 250)
    self.assertEqual(perfect_balance_weights[3], 250)

    max_block_weights = context.compute_max_block_weights(1000)
    self.assertEqual(max_block_weights[0], 262)
    self.assertEqual(max_block_weights[1], 262)
    self.assertEqual(max_block_weights[2], 262)
    self.assertEqual(max_block_weights[3], 262)

    context.set_individual_target_block_weights([100, 200, 300, 400])
    max_block_weights = context.compute_max_block_weights(1000)
    self.assertEqual(max_block_weights[0], 100)
    self.assertEqual(max_block_weights[1], 200)
    self.assertEqual(max_block_weights[2], 300)
    self.assertEqual(max_block_weights[3], 400)

  def test_check_graph_stats(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.num_nodes(), 5)
    self.assertEqual(graph.num_edges(), 12)
    self.assertEqual(graph.num_undirected_edges(), 6)
    self.assertEqual(graph.num_directed_edges(), 12)
    self.assertEqual(graph.total_weight(), 5)

  def test_check_graph_iterators(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual([hn for hn in graph.nodes()], [0,1,2,3,4])
    self.assertEqual([he for he in graph.edges()], [0,1,2,3,4,5,6,7,8,9,10,11])
    self.assertEqual([pin for pin in graph.pins(0)],  [0,1])
    self.assertEqual([pin for pin in graph.pins(1)],  [0,2])
    self.assertEqual([pin for pin in graph.pins(2)],  [1,0])
    self.assertEqual([pin for pin in graph.pins(3)],  [1,2])
    self.assertEqual([pin for pin in graph.pins(4)],  [1,3])
    self.assertEqual([pin for pin in graph.pins(5)],  [2,0])
    self.assertEqual([pin for pin in graph.pins(6)],  [2,1])
    self.assertEqual([pin for pin in graph.pins(7)],  [2,3])
    self.assertEqual([pin for pin in graph.pins(8)],  [3,1])
    self.assertEqual([pin for pin in graph.pins(9)],  [3,2])
    self.assertEqual([pin for pin in graph.pins(10)], [3,4])
    self.assertEqual([pin for pin in graph.pins(11)], [4,3])
    self.assertEqual([he for he in graph.incident_edges(0)], [0,1])
    self.assertEqual([he for he in graph.incident_edges(1)], [2,3,4])
    self.assertEqual([he for he in graph.incident_edges(2)], [5,6,7])
    self.assertEqual([he for he in graph.incident_edges(3)], [8,9,10])
    self.assertEqual([he for he in graph.incident_edges(4)], [11])

  def test_check_graph_node_degrees(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    self.assertEqual(graph.node_degree(0), 2)
    self.assertEqual(graph.node_degree(1), 3)
    self.assertEqual(graph.node_degree(2), 3)
    self.assertEqual(graph.node_degree(3), 3)
    self.assertEqual(graph.node_degree(4), 1)

  def test_check_graph_node_weights(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context,
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,1,1,1,1,1])

    self.assertEqual(graph.total_weight(), 15)
    self.assertEqual(graph.node_weight(0), 1)
    self.assertEqual(graph.node_weight(1), 2)
    self.assertEqual(graph.node_weight(2), 3)
    self.assertEqual(graph.node_weight(3), 4)
    self.assertEqual(graph.node_weight(4), 5)

  def test_check_graph_edge_weights(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context,
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
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context,
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.edge_source(0),  0) # (0,1)
    self.assertEqual(graph.edge_source(1),  0) # (0,2)
    self.assertEqual(graph.edge_source(2),  1) # (1,0)
    self.assertEqual(graph.edge_source(3),  1) # (1,2)
    self.assertEqual(graph.edge_source(4),  1) # (1,3)
    self.assertEqual(graph.edge_source(5),  2) # (2,0)
    self.assertEqual(graph.edge_source(6),  2) # (2,1)
    self.assertEqual(graph.edge_source(7),  2) # (2,3)
    self.assertEqual(graph.edge_source(8),  3) # (3,1)
    self.assertEqual(graph.edge_source(9),  3) # (3,2)
    self.assertEqual(graph.edge_source(10), 3) # (3,4)
    self.assertEqual(graph.edge_source(11), 4) # (4,3)

  def test_target_nodes_of_edges(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context,
      5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)],
      [1,2,3,4,5], [1,2,3,4,5,6])

    self.assertEqual(graph.edge_target(0),  1) # (0,1)
    self.assertEqual(graph.edge_target(1),  2) # (0,2)
    self.assertEqual(graph.edge_target(2),  0) # (1,0)
    self.assertEqual(graph.edge_target(3),  2) # (1,2)
    self.assertEqual(graph.edge_target(4),  3) # (1,3)
    self.assertEqual(graph.edge_target(5),  0) # (2,0)
    self.assertEqual(graph.edge_target(6),  1) # (2,1)
    self.assertEqual(graph.edge_target(7),  3) # (2,3)
    self.assertEqual(graph.edge_target(8),  1) # (3,1)
    self.assertEqual(graph.edge_target(9),  2) # (3,2)
    self.assertEqual(graph.edge_target(10), 4) # (3,4)
    self.assertEqual(graph.edge_target(11), 3) # (4,3)

  def test_graph_applies_bounds_checking(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])

    # call all the methods with invalid IDs
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.node_degree(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.node_weight(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.is_fixed(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.fixed_vertex_block(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.edge_size(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.edge_weight(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.edge_source(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.edge_target(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.incident_edges(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: graph.pins(12))

  def test_load_graph_in_metis_file_format(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.graph_from_file(
      mydir + "/test_instances/delaunay_n15.graph", context, mtkahypar.FileFormat.METIS)

    self.assertEqual(graph.num_nodes(), 32768)
    self.assertEqual(graph.num_edges(), 2 * 98274)
    self.assertEqual(graph.num_undirected_edges(), 98274)
    self.assertEqual(graph.num_directed_edges(), 2 * 98274)
    self.assertEqual(graph.total_weight(), 32768)

  def test_check_hypergraph_stats(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.num_nodes(), 7)
    self.assertEqual(hypergraph.num_edges(), 4)
    self.assertEqual(hypergraph.num_pins(), 12)
    self.assertEqual(hypergraph.total_weight(), 7)

  def test_check_hypergraph_iterators(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual([hn for hn in hypergraph.nodes()], [0,1,2,3,4,5,6])
    self.assertEqual([he for he in hypergraph.edges()], [0,1,2,3])
    self.assertEqual([pin for pin in hypergraph.pins(0)], [0,2])
    self.assertEqual([pin for pin in hypergraph.pins(1)], [0,1,3,4])
    self.assertEqual([pin for pin in hypergraph.pins(2)], [3,4,6])
    self.assertEqual([pin for pin in hypergraph.pins(3)], [2,5,6])
    self.assertEqual([he for he in hypergraph.incident_edges(0)], [0,1])
    self.assertEqual([he for he in hypergraph.incident_edges(1)], [1])
    self.assertEqual([he for he in hypergraph.incident_edges(2)], [0,3])
    self.assertEqual([he for he in hypergraph.incident_edges(3)], [1,2])
    self.assertEqual([he for he in hypergraph.incident_edges(4)], [1,2])
    self.assertEqual([he for he in hypergraph.incident_edges(5)], [3])
    self.assertEqual([he for he in hypergraph.incident_edges(6)], [2,3])

  def test_check_hypergraph_node_degrees(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.node_degree(0), 2)
    self.assertEqual(hypergraph.node_degree(1), 1)
    self.assertEqual(hypergraph.node_degree(2), 2)
    self.assertEqual(hypergraph.node_degree(3), 2)
    self.assertEqual(hypergraph.node_degree(4), 2)
    self.assertEqual(hypergraph.node_degree(5), 1)
    self.assertEqual(hypergraph.node_degree(6), 2)

  def test_check_hypergraph_edge_sizes(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.edge_size(0), 2)
    self.assertEqual(hypergraph.edge_size(1), 4)
    self.assertEqual(hypergraph.edge_size(2), 3)
    self.assertEqual(hypergraph.edge_size(3), 3)

  def test_check_hypergraph_node_weights(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context,
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
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context,
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,2,3,4])

    self.assertEqual(hypergraph.edge_weight(0), 1)
    self.assertEqual(hypergraph.edge_weight(1), 2)
    self.assertEqual(hypergraph.edge_weight(2), 3)
    self.assertEqual(hypergraph.edge_weight(3), 4)

  def test_hypergraph_applies_bounds_checking(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context,
      7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,2,3,4])

    # call all the methods with invalid IDs
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.node_degree(8))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.node_weight(8))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.is_fixed(8))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.fixed_vertex_block(8))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.edge_size(4))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.edge_weight(4))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.incident_edges(8))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: hypergraph.pins(4))

  def test_load_hypergraph_in_hmetis_file_format(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    # default file format is HMETIS
    hypergraph = mtk.hypergraph_from_file(mydir + "/test_instances/ibm01.hgr", context)

    self.assertEqual(hypergraph.num_nodes(), 12752)
    self.assertEqual(hypergraph.num_edges(), 14111)
    self.assertEqual(hypergraph.num_pins(), 50566)
    self.assertEqual(hypergraph.total_weight(), 12752)

  def test_load_hypergraph_in_metis_file_format(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.hypergraph_from_file(
      mydir + "/test_instances/delaunay_n15.graph", context, mtkahypar.FileFormat.METIS)

    self.assertEqual(hypergraph.num_nodes(), 32768)
    self.assertEqual(hypergraph.num_edges(), 98274)
    self.assertEqual(hypergraph.total_weight(), 32768)

  def test_for_graph_if_all_nodes_in_correct_block(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.block_id(0), 0)
    self.assertEqual(partitioned_graph.block_id(1), 1)
    self.assertEqual(partitioned_graph.block_id(2), 1)
    self.assertEqual(partitioned_graph.block_id(3), 2)
    self.assertEqual(partitioned_graph.block_id(4), 2)

  def test_for_graph_if_block_have_correct_weight(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.block_weight(0), 1)
    self.assertEqual(partitioned_graph.block_weight(1), 2)
    self.assertEqual(partitioned_graph.block_weight(2), 2)

  def test_cut_metric_for_graph(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.cut(), 4)

  def test_phg_keeps_graph_alive(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graphs = [mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])]
    partitioned_graph = graphs[0].create_partitioned_hypergraph(context, 3, [0,1,1,2,2])
    graphs = []  # graph is freed if it is not kept alive by the phg

    partitioned_graph.cut()  # segfaults if graph is freed

  def test_for_graph_if_all_nodes_contains_correct_number_of_incident_cut_edges(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

    self.assertEqual(partitioned_graph.num_incident_cut_edges(0), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(1), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(2), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_graph.num_incident_cut_edges(4), 0)

  def test_for_graph_if_all_edges_have_correct_connectivity(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

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

  def test_partitioned_graph_applies_bounds_checking(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,1,1,2,2])

    # call all the methods with invalid IDs
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.block_weight(3))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.block_id(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.is_fixed(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.fixed_vertex_block(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.is_incident_to_cut_edge(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.num_incident_cut_edges(5))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.connectivity(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.connectivity_set(12))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.num_pins_in_block(12, 0))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_graph.num_pins_in_block(0, 3))

  def test_load_graph_partition_from_file(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.partitioned_hypergraph_from_file(context, 3,
      mydir + "/test_instances/test_graph_partition.part3")

    self.assertEqual(partitioned_graph.block_id(0), 0)
    self.assertEqual(partitioned_graph.block_id(1), 0)
    self.assertEqual(partitioned_graph.block_id(2), 1)
    self.assertEqual(partitioned_graph.block_id(3), 1)
    self.assertEqual(partitioned_graph.block_id(4), 2)

  def test_get_and_write_graph_partition_to_file(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    graph = mtk.create_graph(context, 5, 6, [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4)])
    partitioned_graph = graph.create_partitioned_hypergraph(context, 3, [0,0,1,2,2])

    partition = partitioned_graph.get_partition()
    self.assertEqual(partition[0], 0)
    self.assertEqual(partition[1], 0)
    self.assertEqual(partition[2], 1)
    self.assertEqual(partition[3], 2)
    self.assertEqual(partition[4], 2)

    partitioned_graph.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_graph_2 = graph.partitioned_hypergraph_from_file(context, 3,
      mydir + "/test_partition.part3")

    self.assertEqual(partitioned_graph_2.block_id(0), 0)
    self.assertEqual(partitioned_graph_2.block_id(1), 0)
    self.assertEqual(partitioned_graph_2.block_id(2), 1)
    self.assertEqual(partitioned_graph_2.block_id(3), 2)
    self.assertEqual(partitioned_graph_2.block_id(4), 2)

    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

  def test_for_hypergraph_if_all_nodes_are_in_correct_block(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_for_hypergraph_if_blocks_have_correct_weight(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.block_weight(0), 3)
    self.assertEqual(partitioned_hg.block_weight(1), 3)
    self.assertEqual(partitioned_hg.block_weight(2), 1)

  def test_all_metrics_for_hypergraph(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.cut(), 3)
    self.assertEqual(partitioned_hg.km1(), 4)

  def test_phg_keeps_hypergraph_alive(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraphs = [mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])]
    partitioned_hg = hypergraphs[0].create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])
    hypergraphs = []  # hypergraph is freed if it is not kept alive by the phg

    partitioned_hg.cut()  # segfaults if hypergraph is freed

  def test_for_hypergraph_if_all_nodes_contains_correct_number_of_incident_cut_hyperedges(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.num_incident_cut_edges(0), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(1), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(2), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(4), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(5), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(6), 2)

  def test_for_hypergraph_if_all_edges_have_correct_connectivity(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual(partitioned_hg.connectivity(0), 1)
    self.assertEqual(partitioned_hg.connectivity(1), 2)
    self.assertEqual(partitioned_hg.connectivity(2), 2)
    self.assertEqual(partitioned_hg.connectivity(3), 3)

  def test_for_hypergraph_if_all_edges_have_correct_number_of_pins_in_blocks(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

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

  def test_partitioned_hypergraph_iterators(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    self.assertEqual([block for block in partitioned_hg.blocks()], [0,1,2])
    self.assertEqual([block for block in partitioned_hg.connectivity_set(0)], [0])
    self.assertEqual([block for block in partitioned_hg.connectivity_set(1)], [0,1])
    self.assertEqual([block for block in partitioned_hg.connectivity_set(2)], [1,2])
    self.assertEqual([block for block in partitioned_hg.connectivity_set(3)], [0,1,2])

  def test_partitioned_hypergraph_applies_bounds_checking(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])

    # call all the methods with invalid IDs
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.block_weight(3))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.block_id(7))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.is_fixed(7))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.fixed_vertex_block(7))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.is_incident_to_cut_edge(7))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.num_incident_cut_edges(7))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.connectivity(4))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.connectivity_set(4))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.num_pins_in_block(4, 0))
    self.assertRaises(mtkahypar.InvalidInputError, lambda: partitioned_hg.num_pins_in_block(0, 3))

  def test_load_hypergraph_partition_from_file(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.partitioned_hypergraph_from_file(context, 3,
      mydir + "/test_instances/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_get_and_write_hypergraph_partition_to_file(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,2,2])

    partition = partitioned_hg.get_partition()
    self.assertEqual(partition[0], 0)
    self.assertEqual(partition[1], 0)
    self.assertEqual(partition[2], 0)
    self.assertEqual(partition[3], 1)
    self.assertEqual(partition[4], 1)
    self.assertEqual(partition[5], 2)
    self.assertEqual(partition[6], 2)

    partitioned_hg.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_hg_2 = hypergraph.partitioned_hypergraph_from_file(context, 3,
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


  def create_sparse_phg(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.LARGE_K)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,1,2])
    self.assertTrue(partitioned_hg.is_compatible(mtkahypar.PresetType.LARGE_K))
    self.assertFalse(partitioned_hg.is_compatible(mtkahypar.PresetType.DEFAULT))
    return partitioned_hg

  def test_for_sparse_hypergraph_if_all_nodes_are_in_correct_block(self):
    partitioned_hg = self.create_sparse_phg()
    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_for_sparse_hypergraph_if_blocks_have_correct_weight(self):
    partitioned_hg = self.create_sparse_phg()
    self.assertEqual(partitioned_hg.block_weight(0), 3)
    self.assertEqual(partitioned_hg.block_weight(1), 3)
    self.assertEqual(partitioned_hg.block_weight(2), 1)

  def test_all_metrics_for_sparse_hypergraph(self):
    partitioned_hg = self.create_sparse_phg()
    self.assertEqual(partitioned_hg.cut(), 3)
    self.assertEqual(partitioned_hg.km1(), 4)

  def test_for_sparse_hypergraph_if_all_nodes_contains_correct_number_of_incident_cut_hyperedges(self):
    partitioned_hg = self.create_sparse_phg()
    self.assertEqual(partitioned_hg.num_incident_cut_edges(0), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(1), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(2), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(3), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(4), 2)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(5), 1)
    self.assertEqual(partitioned_hg.num_incident_cut_edges(6), 2)

  def test_for_sparse_hypergraph_if_all_edges_have_correct_connectivity(self):
    partitioned_hg = self.create_sparse_phg()
    self.assertEqual(partitioned_hg.connectivity(0), 1)
    self.assertEqual(partitioned_hg.connectivity(1), 2)
    self.assertEqual(partitioned_hg.connectivity(2), 2)
    self.assertEqual(partitioned_hg.connectivity(3), 3)

  def test_for_sparse_hypergraph_if_all_edges_have_correct_number_of_pins_in_blocks(self):
    partitioned_hg = self.create_sparse_phg()
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
    context = mtk.context_from_preset(mtkahypar.PresetType.LARGE_K)
    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.partitioned_hypergraph_from_file(context, 3,
      mydir + "/test_instances/test_partition.part3")

    self.assertEqual(partitioned_hg.block_id(0), 0)
    self.assertEqual(partitioned_hg.block_id(1), 0)
    self.assertEqual(partitioned_hg.block_id(2), 0)
    self.assertEqual(partitioned_hg.block_id(3), 1)
    self.assertEqual(partitioned_hg.block_id(4), 1)
    self.assertEqual(partitioned_hg.block_id(5), 1)
    self.assertEqual(partitioned_hg.block_id(6), 2)

  def test_write_sparse_hypergraph_partition_to_file(self):
    context = mtk.context_from_preset(mtkahypar.PresetType.LARGE_K)
    if os.path.isfile(mydir + "/test_partition.part3"):
      os.remove(mydir + "/test_partition.part3")

    hypergraph = mtk.create_hypergraph(context, 7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])
    partitioned_hg = hypergraph.create_partitioned_hypergraph(context, 3, [0,0,0,1,1,2,2])

    partitioned_hg.write_partition_to_file(mydir + "/test_partition.part3")
    partitioned_hg_2 = hypergraph.partitioned_hypergraph_from_file(context, 3,
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
      self.context = mtk.context_from_preset(preset_type)
      self.context.set_partitioning_parameters(num_blocks, epsilon, objective)
      mtkahypar.set_seed(42)
      self.context.logging = logging or force_logging
      self.target_graph = mtk.target_graph_from_file(mydir + "/test_instances/target.graph", self.context)
      self.graph = mtk.graph_from_file(mydir + "/test_instances/delaunay_n15.graph", self.context)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.graph.total_weight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.set_individual_target_block_weights(individualBlockWeights)

    def addFixedVertices(self):
      self.graph.add_fixed_vertices_from_file(mydir + "/test_instances/delaunay_n15.k4.p1.fix", self.k)

    def partition(self):
      self.assertTrue(self.graph.is_compatible(self.context.preset))
      self.partitioned_graph = self.graph.partition(self.context)
      self.__verifyPartition()

    def map_onto_graph(self):
      self.assertTrue(self.graph.is_compatible(self.context.preset))
      self.partitioned_graph = self.graph.map_onto_graph(self.target_graph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      self.assertTrue(self.partitioned_graph.is_compatible(self.context.preset))
      objective_before = self.partitioned_graph.cut()
      self.partitioned_graph.improve_partition(self.context, num_vcycles)
      objective_after = self.partitioned_graph.cut()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def improveMapping(self, num_vcycles):
      self.assertTrue(self.partitioned_graph.is_compatible(self.context.preset))
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
      for hn in self.graph.nodes():
        self.assertGreaterEqual(self.partitioned_graph.block_id(hn), 0),
        self.assertLess(self.partitioned_graph.block_id(hn), self.k)

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

  def test_partitions_a_graph_with_highest_quality_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_highest_quality_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_quality_preset_into_two_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 2, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()

  def test_partitions_a_graph_with_deterministic_quality_preset_into_four_blocks(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 4, 0.03, mtkahypar.Objective.CUT, False)
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

  def test_checks_if_deterministic_quality_preset_produces_same_result_for_graph(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 8, 0.03, mtkahypar.Objective.CUT, False)
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
    partitioner.context = mtk.context_from_preset(mtkahypar.PresetType.QUALITY)
    partitioner.context.set_partitioning_parameters(4, 0.03, mtkahypar.Objective.CUT)
    partitioner.improvePartition(1)

  def test_improves_a_graph_partition_with_three_vcycle(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.DEFAULT, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    partitioner.improvePartition(3)

  def test_improves_a_graph_partition_with_one_vcycle_and_highest_quality_preset(self):
    partitioner = self.GraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 4, 0.03, mtkahypar.Objective.CUT, False)
    partitioner.partition()
    partitioner.improvePartition(1)

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
      self.context = mtk.context_from_preset(preset_type)
      self.context.set_partitioning_parameters(num_blocks, epsilon, objective)
      mtkahypar.set_seed(42)
      self.context.logging = logging or force_logging
      self.target_graph = mtk.target_graph_from_file(mydir + "/test_instances/target.graph", self.context)
      self.hypergraph = mtk.hypergraph_from_file(mydir + "/test_instances/ibm01.hgr", self.context)
      self.useIndividualBlockWeights = False
      self.k = num_blocks
      self.epsilon = epsilon
      self.preset_type = preset_type
      self.maxAllowedBlockWeight = math.floor((1.0 + epsilon) *
        math.ceil(self.hypergraph.total_weight() / float(num_blocks)))

    def setIndividualBlockWeights(self, individualBlockWeights):
      self.useIndividualBlockWeights = True
      self.individualBlockWeights = individualBlockWeights
      self.context.set_individual_target_block_weights(individualBlockWeights)

    def addFixedVertices(self):
      self.hypergraph.add_fixed_vertices_from_file(mydir + "/test_instances/ibm01.k4.p1.fix", self.k)

    def partition(self):
      self.assertTrue(self.hypergraph.is_compatible(self.context.preset))
      self.partitioned_hg = self.hypergraph.partition(self.context)
      self.__verifyPartition()

    def map_onto_graph(self):
      self.assertTrue(self.hypergraph.is_compatible(self.context.preset))
      self.partitioned_hg = self.hypergraph.map_onto_graph(self.target_graph, self.context)
      self.__verifyPartition()

    def improvePartition(self, num_vcycles):
      self.assertTrue(self.partitioned_hg.is_compatible(self.context.preset))
      objective_before = self.partitioned_hg.km1()
      self.partitioned_hg.improve_partition(self.context, num_vcycles)
      objective_after = self.partitioned_hg.km1()
      self.assertLessEqual(objective_after, objective_before)
      self.__verifyPartition()

    def improveMapping(self, num_vcycles):
      self.assertTrue(self.partitioned_hg.is_compatible(self.context.preset))
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
      for hn in self.hypergraph.nodes():
        self.assertGreaterEqual(self.partitioned_hg.block_id(hn), 0),
        self.assertLess(self.partitioned_hg.block_id(hn), self.k)

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

  def test_partitions_a_hypergraph_with_highest_quality_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_highest_quality_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_quality_preset_into_two_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 2, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()

  def test_partitions_a_hypergraph_with_deterministic_quality_preset_into_four_blocks(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 4, 0.03, mtkahypar.Objective.KM1, False)
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

  def test_checks_if_deterministic_quality_preset_produces_same_result_for_hypergraphs(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.DETERMINISTIC_QUALITY, 8, 0.03, mtkahypar.Objective.KM1, False)
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
    partitioner.context = mtk.context_from_preset(mtkahypar.PresetType.QUALITY)
    partitioner.context.set_partitioning_parameters(4, 0.03, mtkahypar.Objective.KM1)
    partitioner.improvePartition(1)

  def test_improves_a_hypergraph_partition_with_one_vcycle_and_highest_quality_preset(self):
    partitioner = self.HypergraphPartitioner(mtkahypar.PresetType.HIGHEST_QUALITY, 4, 0.03, mtkahypar.Objective.KM1, False)
    partitioner.partition()
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
