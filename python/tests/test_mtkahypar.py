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

import mtkahypar as mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

class MainTest(unittest.TestCase):

  def test_check_hypergraph_stats(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.numNodes(), 7)
    self.assertEqual(hypergraph.numEdges(), 4)
    self.assertEqual(hypergraph.numPins(), 12)
    self.assertEqual(hypergraph.totalWeight(), 7)

  def test_check_node_degrees(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.nodeDegree(0), 2)
    self.assertEqual(hypergraph.nodeDegree(1), 1)
    self.assertEqual(hypergraph.nodeDegree(2), 2)
    self.assertEqual(hypergraph.nodeDegree(3), 2)
    self.assertEqual(hypergraph.nodeDegree(4), 2)
    self.assertEqual(hypergraph.nodeDegree(5), 1)
    self.assertEqual(hypergraph.nodeDegree(6), 2)

  def test_check_edge_sizes(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]])

    self.assertEqual(hypergraph.edgeSize(0), 2)
    self.assertEqual(hypergraph.edgeSize(1), 4)
    self.assertEqual(hypergraph.edgeSize(2), 3)
    self.assertEqual(hypergraph.edgeSize(3), 3)

  def test_check_node_weights(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
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
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(7, 4, [[0,2],[0,1,3,4],[3,4,6],[2,5,6]],
      [1,2,3,4,5,6,7], [1,2,3,4])

    self.assertEqual(hypergraph.edgeWeight(0), 1)
    self.assertEqual(hypergraph.edgeWeight(1), 2)
    self.assertEqual(hypergraph.edgeWeight(2), 3)
    self.assertEqual(hypergraph.edgeWeight(3), 4)

  def test_load_hypergraph_in_hmetis_file_format(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(mydir + "/test_instances/ibm01.hgr", mtkahypar.FileFormat.HMETIS)

    self.assertEqual(hypergraph.numNodes(), 12752)
    self.assertEqual(hypergraph.numEdges(), 14111)
    self.assertEqual(hypergraph.numPins(), 50566)
    self.assertEqual(hypergraph.totalWeight(), 12752)

  def test_load_hypergraph_in_metis_file_format(self):
    hypergraph = mtkahypar.Hypergraph()
    hypergraph.construct(mydir + "/test_instances/delaunay_n15.graph", mtkahypar.FileFormat.METIS)

    self.assertEqual(hypergraph.numNodes(), 32768)
    self.assertEqual(hypergraph.numEdges(), 98274)
    self.assertEqual(hypergraph.totalWeight(), 32768)

if __name__ == '__main__':
  unittest.main()
