/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/separated_nodes.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using Edge = SeparatedNodes::Edge;

class ASeparatedNodes: public Test { };

void verifyIncidentEgdes(IteratorRange<SeparatedNodes::IncidenceIterator> it,
                         const std::set<std::pair<HypernodeID, HyperedgeWeight>>& edges) {
  size_t count = 0;
  for (const Edge& edge : it) {
    ASSERT_TRUE(edges.find({edge.target, edge.weight}) != edges.end()) << V(edge.target);
    count++;
  }
  ASSERT_EQ(count, edges.size());
}

TEST_F(ASeparatedNodes, HasCorrectStats) {
  SeparatedNodes nodes(10);
  ASSERT_EQ(0,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(0,  nodes.numEdges());
}

TEST_F(ASeparatedNodes, AddsNodes1) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {2, 2}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(2,  nodes.numEdges());

  ASSERT_EQ(1,  nodes.nodeWeight(0));
  ASSERT_EQ(1,  nodes.nodeWeight(1));
  ASSERT_EQ(2,  nodes.nodeWeight(2));

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(3));
}

TEST_F(ASeparatedNodes, AddsNodes2) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {1, 1}};
  vec<Edge> new_edges { Edge(8, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  new_nodes = { {3, 1}, {3, 1}, {4, 1}};
  new_edges = { Edge(1, 1), Edge(8, 1), Edge(2, 1), Edge(8, 1) };
  nodes.addNodes(new_nodes, new_edges);

  ASSERT_EQ(5,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(6,  nodes.numEdges());

  ASSERT_EQ(1,  nodes.nodeWeight(0));
  ASSERT_EQ(1,  nodes.nodeWeight(1));
  ASSERT_EQ(1,  nodes.nodeWeight(2));
  ASSERT_EQ(1,  nodes.nodeWeight(3));
  ASSERT_EQ(1,  nodes.nodeWeight(4));

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(2));
  ASSERT_EQ(3,  nodes.outwardIncidentWeight(8));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(9));
}

TEST_F(ASeparatedNodes, HasCorrectNodeIterator) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  HypernodeID expected_node = 0;
  for ( const HypernodeID& node : nodes.nodes() ) {
    ASSERT_EQ(expected_node++, node);
  }
  ASSERT_EQ(3, expected_node);
}

TEST_F(ASeparatedNodes, HasCorrectInwardEdgeIterator) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  verifyIncidentEgdes(nodes.inwardEdges(0), { {0, 1}, {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, MapsNodes1) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 1, 1, 0, 2, kInvalidHypernode, 2, 0, 0, kInvalidHypernode };
  nodes.contract(communities, 3);

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(3,  nodes.numGraphNodes());
  ASSERT_EQ(1,  nodes.numEdges());

  ASSERT_EQ(2,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(2));

  verifyIncidentEgdes(nodes.inwardEdges(0), { {0, 2} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { });
  verifyIncidentEgdes(nodes.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, MapsNodes2) {
  SeparatedNodes nodes(10);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {3, 1}, {4, 1}, {6, 1}};
  vec<Edge> new_edges { Edge(3, 2), Edge(1, 1), Edge(1, 2), Edge(1, 1), Edge(3, 1), Edge(4, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 3, 1, 3, kInvalidHypernode, kInvalidHypernode, 0, 1, 4, kInvalidHypernode };
  nodes.contract(communities, 5);

  ASSERT_EQ(5,  nodes.numNodes());
  ASSERT_EQ(5,  nodes.numGraphNodes());
  ASSERT_EQ(4,  nodes.numEdges());

  ASSERT_EQ(0,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(2));
  ASSERT_EQ(7,  nodes.outwardIncidentWeight(3));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(4));

  verifyIncidentEgdes(nodes.inwardEdges(0), { {3, 3} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { {3, 2} });
  verifyIncidentEgdes(nodes.inwardEdges(2), { {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(3), { {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(4), { });
}

TEST_F(ASeparatedNodes, InitializesEdges) {
  SeparatedNodes nodes(5);
  vec<std::pair<HyperedgeID, HypernodeWeight>> new_nodes { {0, 1}, {2, 1}, {4, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  nodes.initializeOutwardEdges();

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(5,  nodes.numGraphNodes());
  ASSERT_EQ(4,  nodes.numEdges());

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(2));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(3));
  ASSERT_EQ(2,  nodes.outwardIncidentWeight(4));

  verifyIncidentEgdes(nodes.outwardEdges(0), { {0, 1} });
  verifyIncidentEgdes(nodes.outwardEdges(1), { });
  verifyIncidentEgdes(nodes.outwardEdges(2), { });
  verifyIncidentEgdes(nodes.outwardEdges(3), { {1, 1} });
  verifyIncidentEgdes(nodes.outwardEdges(4), { {0, 1}, {1, 1} });
}

}
} // namespace mt_kahypar
