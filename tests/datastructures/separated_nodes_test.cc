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

#include <algorithm>

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

void verifyGraphNeighbors(const Hypergraph& graph, const HypernodeID& node, vec<HypernodeID> expected) {
  std::sort(expected.begin(), expected.end());
  vec<HypernodeID> neighbors;
  for (HyperedgeID e: graph.incidentEdges(node)) {
    neighbors.push_back(graph.edgeTarget(e));
  }
  std::sort(neighbors.begin(), neighbors.end());
  ASSERT_EQ(neighbors.size(), expected.size());
  for (size_t i = 0; i < neighbors.size(); ++i) {
    ASSERT_EQ(neighbors[i], expected[i]);
  }
}

TEST_F(ASeparatedNodes, HasCorrectStats) {
  SeparatedNodes nodes(10);
  ASSERT_EQ(0,  nodes.numNodes());
  ASSERT_EQ(0,  nodes.numVisibleNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(0,  nodes.numEdges());
  ASSERT_EQ(0,  nodes.currentBatchIndex());
  ASSERT_EQ(0,  nodes.totalWeight());
}

TEST_F(ASeparatedNodes, AddsNodes1) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 2, 2}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  nodes.revealAll();

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(2,  nodes.numEdges());
  ASSERT_EQ(4,  nodes.totalWeight());

  ASSERT_EQ(1,  nodes.nodeWeight(0));
  ASSERT_EQ(1,  nodes.nodeWeight(1));
  ASSERT_EQ(2,  nodes.nodeWeight(2));

  ASSERT_EQ(0,  nodes.originalHypernodeID(0));
  ASSERT_EQ(1,  nodes.originalHypernodeID(1));
  ASSERT_EQ(2,  nodes.originalHypernodeID(2));

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(3));

  ASSERT_EQ(0,  nodes.currentBatchIndex());
  ASSERT_EQ(3,  nodes.numVisibleNodes());
}

TEST_F(ASeparatedNodes, AddsNodes2) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 1, 1}};
  vec<Edge> new_edges { Edge(8, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  new_nodes = { {2, 0, 1}, {3, 3, 1}, {4, 4, 1}};
  new_edges = { Edge(1, 1), Edge(8, 1), Edge(2, 1), Edge(8, 1) };
  nodes.addNodes(new_nodes, new_edges);
  nodes.revealNextBatch();

  ASSERT_EQ(5,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(6,  nodes.numEdges());
  ASSERT_EQ(5,  nodes.totalWeight());

  ASSERT_EQ(1,  nodes.nodeWeight(0));
  ASSERT_EQ(1,  nodes.nodeWeight(1));
  ASSERT_EQ(1,  nodes.nodeWeight(2));
  ASSERT_EQ(1,  nodes.nodeWeight(3));
  ASSERT_EQ(1,  nodes.nodeWeight(4));

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(2));
  ASSERT_EQ(3,  nodes.outwardIncidentWeight(8));
  ASSERT_EQ(1,  nodes.outwardIncidentWeight(9));

  ASSERT_EQ(2,  nodes.currentBatchIndex());
  ASSERT_EQ(2,  nodes.numVisibleNodes());
}

TEST_F(ASeparatedNodes, PopsBatches) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 1, 1}};
  vec<Edge> new_edges { Edge(8, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  new_nodes = { {2, 0, 1}, {3, 3, 1}, {4, 4, 1}};
  new_edges = { Edge(1, 1), Edge(8, 1), Edge(2, 1), Edge(8, 1) };
  nodes.addNodes(new_nodes, new_edges);

  ASSERT_EQ(2,  nodes.popBatch());

  ASSERT_EQ(2,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(2,  nodes.numEdges());
  ASSERT_EQ(2,  nodes.totalWeight());

  ASSERT_EQ(1,  nodes.nodeWeight(0));
  ASSERT_EQ(1,  nodes.nodeWeight(1));
  
  // TODO: incident weight?

  ASSERT_EQ(0,  nodes.currentBatchIndex());
}

TEST_F(ASeparatedNodes, HasCorrectNodeIterator) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}};
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
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  verifyIncidentEgdes(nodes.inwardEdges(0), { {0, 1}, {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, MapsNodes1) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 1, 1, 0, 2, kInvalidHypernode, 2, 0, 0, kInvalidHypernode };
  nodes.contract(communities, 3);

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(3,  nodes.numGraphNodes());
  ASSERT_EQ(1,  nodes.numEdges());
  ASSERT_EQ(3,  nodes.totalWeight());

  ASSERT_EQ(2,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(1));
  ASSERT_EQ(0,  nodes.outwardIncidentWeight(2));

  verifyIncidentEgdes(nodes.inwardEdges(0), { {0, 2} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { });
  verifyIncidentEgdes(nodes.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, MapsNodes2) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes
          { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 6, 1}};
  vec<Edge> new_edges { Edge(3, 2), Edge(1, 1), Edge(1, 2), Edge(1, 1), Edge(3, 1), Edge(4, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 3, 1, 3, kInvalidHypernode, kInvalidHypernode, 0, 1, 4, kInvalidHypernode };
  nodes.contract(communities, 5);

  ASSERT_EQ(5,  nodes.numNodes());
  ASSERT_EQ(5,  nodes.numGraphNodes());
  ASSERT_EQ(4,  nodes.numEdges());
  ASSERT_EQ(5,  nodes.totalWeight());

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

TEST_F(ASeparatedNodes, Coarsens1) {
  SeparatedNodes nodes(4);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}};
  vec<Edge> new_edges { Edge(2, 1), Edge(1, 1), Edge(0, 1), Edge(2, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 1, 2, 2 };
  nodes.revealAll();
  SeparatedNodes other = nodes.coarsen(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(1, communities[1]);
  ASSERT_EQ(1, communities[2]);

  ASSERT_EQ(2,  other.numNodes());
  ASSERT_EQ(4,  other.numGraphNodes());
  ASSERT_EQ(4,  other.numEdges());
  ASSERT_EQ(3,  other.totalWeight());
  ASSERT_EQ(2,  other.numVisibleNodes());

  ASSERT_EQ(1,  other.outwardIncidentWeight(0));
  ASSERT_EQ(1,  other.outwardIncidentWeight(1));
  ASSERT_EQ(2,  other.outwardIncidentWeight(2));
  ASSERT_EQ(0,  other.outwardIncidentWeight(3));

  verifyIncidentEgdes(other.inwardEdges(0), { {2, 1}, {1, 1} });
  verifyIncidentEgdes(other.inwardEdges(1), { {0, 1}, {2, 1} });

  communities = { 0, 0 };
  SeparatedNodes next = other.coarsen(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(0, communities[1]);

  ASSERT_EQ(1,  next.numNodes());
  ASSERT_EQ(4,  next.numGraphNodes());
  ASSERT_EQ(3,  next.numEdges());
  ASSERT_EQ(3,  next.totalWeight());
  ASSERT_EQ(1,  next.numVisibleNodes());

  ASSERT_EQ(1,  next.outwardIncidentWeight(0));
  ASSERT_EQ(1,  next.outwardIncidentWeight(1));
  ASSERT_EQ(2,  next.outwardIncidentWeight(2));
  ASSERT_EQ(0,  next.outwardIncidentWeight(3));

  verifyIncidentEgdes(next.inwardEdges(0), { {0, 1}, {1, 1}, {2, 2} });
}

TEST_F(ASeparatedNodes, Coarsens2) {
  SeparatedNodes nodes(4);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes
           { {0, 0, 1}, {1, 3, 1}, {2, 5, 1}, {3, 5, 1}, {4, 7, 1} };
  vec<Edge> new_edges { Edge(3, 1), Edge(1, 1), Edge(0, 1), Edge(1, 1),
                        Edge(3, 1), Edge(0, 1), Edge(1, 1), Edge(0, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 0, 1, 2, 2 };
  nodes.revealAll();
  SeparatedNodes other = nodes.coarsen(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(0, communities[1]);
  ASSERT_EQ(1, communities[2]);
  ASSERT_EQ(2, communities[3]);
  ASSERT_EQ(2, communities[4]);

  ASSERT_EQ(3,  other.numNodes());
  ASSERT_EQ(4,  other.numGraphNodes());
  ASSERT_EQ(5,  other.numEdges());
  ASSERT_EQ(5,  other.totalWeight());
  ASSERT_EQ(3,  other.numVisibleNodes());

  ASSERT_EQ(3,  other.outwardIncidentWeight(0));
  ASSERT_EQ(3,  other.outwardIncidentWeight(1));
  ASSERT_EQ(0,  other.outwardIncidentWeight(2));
  ASSERT_EQ(2,  other.outwardIncidentWeight(3));

  verifyIncidentEgdes(other.inwardEdges(0), { {0, 1}, {1, 2}, {3, 2} });
  verifyIncidentEgdes(other.inwardEdges(1), { });
  verifyIncidentEgdes(other.inwardEdges(2), { {0, 2}, {1, 1} });
}

TEST_F(ASeparatedNodes, CoarsensWithHiddenNodes) {
  SeparatedNodes nodes(4);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes
           { {0, 0, 1}, {1, 3, 1}, {2, 5, 1} };
  vec<Edge> new_edges { Edge(3, 1), Edge(1, 1), Edge(0, 1), Edge(1, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  new_nodes = { {3, 0, 1}, {4, 2, 1} };
  new_edges = { Edge(0, 1), Edge(1, 1), Edge(0, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 0, 1 };
  nodes.revealNextBatch();
  SeparatedNodes other = nodes.coarsen(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(0, communities[1]);
  ASSERT_EQ(1, communities[2]);

  ASSERT_EQ(4,  other.numNodes());
  ASSERT_EQ(4,  other.numGraphNodes());
  ASSERT_EQ(6,  other.numEdges());
  ASSERT_EQ(5,  other.totalWeight());
  ASSERT_EQ(2,  other.numVisibleNodes());

  ASSERT_EQ(3,  other.outwardIncidentWeight(0));
  ASSERT_EQ(3,  other.outwardIncidentWeight(1));
  ASSERT_EQ(0,  other.outwardIncidentWeight(2));
  ASSERT_EQ(2,  other.outwardIncidentWeight(3));

  verifyIncidentEgdes(other.inwardEdges(0), { {0, 1}, {1, 2}, {3, 2} });
  verifyIncidentEgdes(other.inwardEdges(1), { });
  verifyIncidentEgdes(other.inwardEdges(2), { {0, 1}, {1, 1} });
  verifyIncidentEgdes(other.inwardEdges(3), { {0, 1} });
}

TEST_F(ASeparatedNodes, ContractsStack) {
  SeparatedNodes nodes(4);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
  vec<Edge> new_edges { Edge(2, 1), Edge(1, 1), Edge(0, 1), Edge(2, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  nodes.revealAll();

  SepNodesStack stack(std::move(nodes));
  stack.coarsen({ 1, 2, 2, 0 });
  stack.coarsen({ 0, 0, 1 });

  ASSERT_EQ(3, stack.numLevels());
  stack.contractToNLevels(2);
  ASSERT_EQ(2, stack.numLevels());

  const SeparatedNodes& coarse = stack.coarsest();
  ASSERT_EQ(2,  coarse.numNodes());
  ASSERT_EQ(4,  coarse.numGraphNodes());
  ASSERT_EQ(5,  coarse.numEdges());
  ASSERT_EQ(4,  coarse.totalWeight());
  ASSERT_EQ(2,  coarse.numVisibleNodes());

  ASSERT_EQ(1,  coarse.outwardIncidentWeight(0));
  ASSERT_EQ(1,  coarse.outwardIncidentWeight(1));
  ASSERT_EQ(2,  coarse.outwardIncidentWeight(2));
  ASSERT_EQ(1,  coarse.outwardIncidentWeight(3));

  verifyIncidentEgdes(coarse.inwardEdges(0), { {2, 1}, {1, 1}, {3, 1} });
  verifyIncidentEgdes(coarse.inwardEdges(1), { {0, 1}, {2, 1} });

  const vec<HypernodeID>& mapping = stack.mapping(0);
  ASSERT_EQ(4, mapping.size());
  ASSERT_EQ(0, mapping[0]);
  ASSERT_EQ(1, mapping[1]);
  ASSERT_EQ(1, mapping[2]);
  ASSERT_EQ(0, mapping[3]);
}

TEST_F(ASeparatedNodes, InitializesEdges) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 4, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  nodes.initializeOutwardEdges();

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(5,  nodes.numGraphNodes());
  ASSERT_EQ(4,  nodes.numEdges());
  ASSERT_EQ(3,  nodes.totalWeight());

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

TEST_F(ASeparatedNodes, CopiesSequential) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 4, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  nodes.initializeOutwardEdges();

  SeparatedNodes other = nodes.copy();

  ASSERT_EQ(3,  other.numNodes());
  ASSERT_EQ(5,  other.numGraphNodes());
  ASSERT_EQ(4,  other.numEdges());
  ASSERT_EQ(3,  other.totalWeight());

  verifyIncidentEgdes(other.inwardEdges(0), { {0, 1}, {4, 1} });
  verifyIncidentEgdes(other.inwardEdges(1), { {4, 1}, {3, 1} });
  verifyIncidentEgdes(other.inwardEdges(2), { });

  ASSERT_EQ(1,  other.outwardIncidentWeight(0));
  ASSERT_EQ(0,  other.outwardIncidentWeight(1));
  ASSERT_EQ(0,  other.outwardIncidentWeight(2));
  ASSERT_EQ(1,  other.outwardIncidentWeight(3));
  ASSERT_EQ(2,  other.outwardIncidentWeight(4));

  verifyIncidentEgdes(other.outwardEdges(0), { {0, 1} });
  verifyIncidentEgdes(other.outwardEdges(1), { });
  verifyIncidentEgdes(other.outwardEdges(2), { });
  verifyIncidentEgdes(other.outwardEdges(3), { {1, 1} });
  verifyIncidentEgdes(other.outwardEdges(4), { {0, 1}, {1, 1} });
}

TEST_F(ASeparatedNodes, CopiesParallel) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 4, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  nodes.initializeOutwardEdges();

  SeparatedNodes other = nodes.copy(parallel_tag_t());

  ASSERT_EQ(3,  other.numNodes());
  ASSERT_EQ(5,  other.numGraphNodes());
  ASSERT_EQ(4,  other.numEdges());
  ASSERT_EQ(3,  other.totalWeight());

  verifyIncidentEgdes(other.inwardEdges(0), { {0, 1}, {4, 1} });
  verifyIncidentEgdes(other.inwardEdges(1), { {4, 1}, {3, 1} });
  verifyIncidentEgdes(other.inwardEdges(2), { });

  ASSERT_EQ(1,  other.outwardIncidentWeight(0));
  ASSERT_EQ(0,  other.outwardIncidentWeight(1));
  ASSERT_EQ(0,  other.outwardIncidentWeight(2));
  ASSERT_EQ(1,  other.outwardIncidentWeight(3));
  ASSERT_EQ(2,  other.outwardIncidentWeight(4));

  verifyIncidentEgdes(other.outwardEdges(0), { {0, 1} });
  verifyIncidentEgdes(other.outwardEdges(1), { });
  verifyIncidentEgdes(other.outwardEdges(2), { });
  verifyIncidentEgdes(other.outwardEdges(3), { {1, 1} });
  verifyIncidentEgdes(other.outwardEdges(4), { {0, 1}, {1, 1} });
}

TEST_F(ASeparatedNodes, PerformsAlternateAddingAndMapping) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 4, 1}, {3, 6, 1}};
  vec<Edge> new_edges { Edge(3, 2), Edge(1, 1), Edge(2, 1), Edge(1, 2), Edge(3, 1), Edge(4, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeID> communities { 0, 3, 2, 3, kInvalidHypernode, 4, 5, 0, 0, 0 };
  nodes.contract(communities, 6);

  new_nodes = { {4, 0, 1}, {5, 4, 1}, {6, 5, 1} };
  new_edges = { Edge(0, 1), Edge(1, 1), Edge(2, 1), Edge(3, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  communities = { 0, kInvalidHypernode, 1, 1, kInvalidHypernode, kInvalidHypernode };
  nodes.contract(communities, 2);

  nodes.initializeOutwardEdges();

  ASSERT_EQ(7,  nodes.numNodes());
  ASSERT_EQ(2,  nodes.numGraphNodes());
  ASSERT_EQ(6,  nodes.numEdges());
  ASSERT_EQ(7,  nodes.totalWeight());

  verifyIncidentEgdes(nodes.inwardEdges(0), { {1, 3} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { {1, 3} });
  verifyIncidentEgdes(nodes.inwardEdges(2), { {1, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(3), { });
  verifyIncidentEgdes(nodes.inwardEdges(4), { {0, 1}, {1, 2} });
  verifyIncidentEgdes(nodes.inwardEdges(5), { {1, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(6), { });

  ASSERT_EQ(1,  nodes.outwardIncidentWeight(0));
  ASSERT_EQ(10,  nodes.outwardIncidentWeight(1));

  verifyIncidentEgdes(nodes.outwardEdges(0), { {4, 1} });
  verifyIncidentEgdes(nodes.outwardEdges(1), { {0, 3}, {1, 3}, {2, 1}, {4, 2}, {5, 1} });
}

TEST_F(ASeparatedNodes, ExtractsBlockWhileGraphNodesAreKept) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 1, 1}, {2, 2, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  vec<CAtomic<PartitionID>> part_ids;
  part_ids.emplace_back(0);
  part_ids.emplace_back(0);
  part_ids.emplace_back(1);

  SeparatedNodes extracted = nodes.extract(0, {0, 1, 2, 3, 4}, part_ids);

  ASSERT_EQ(2,  extracted.numNodes());
  ASSERT_EQ(5,  extracted.numGraphNodes());
  ASSERT_EQ(2,  extracted.numEdges());
  ASSERT_EQ(2,  extracted.totalWeight());

  verifyIncidentEgdes(extracted.inwardEdges(0), { {0, 1} });
  verifyIncidentEgdes(extracted.inwardEdges(1), { {4, 1} });

  ASSERT_EQ(1,  extracted.outwardIncidentWeight(0));
  ASSERT_EQ(0,  extracted.outwardIncidentWeight(1));
  ASSERT_EQ(0,  extracted.outwardIncidentWeight(2));
  ASSERT_EQ(0,  extracted.outwardIncidentWeight(3));
  ASSERT_EQ(1,  extracted.outwardIncidentWeight(4));
}

TEST_F(ASeparatedNodes, ExtractsBlockWhileRemovingGraphNodes) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 1, 1}, {2, 2, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);
  vec<CAtomic<PartitionID>> part_ids;
  part_ids.emplace_back(0);
  part_ids.emplace_back(0);
  part_ids.emplace_back(1);

  SeparatedNodes extracted = nodes.extract(0, {0, 1, kInvalidHypernode, 2, kInvalidHypernode}, part_ids);

  ASSERT_EQ(2,  extracted.numNodes());
  ASSERT_EQ(3,  extracted.numGraphNodes());
  ASSERT_EQ(1,  extracted.numEdges());
  ASSERT_EQ(2,  extracted.totalWeight());

  verifyIncidentEgdes(extracted.inwardEdges(0), { {0, 1} });
  verifyIncidentEgdes(extracted.inwardEdges(1), { });

  ASSERT_EQ(1,  extracted.outwardIncidentWeight(0));
  ASSERT_EQ(0,  extracted.outwardIncidentWeight(1));
  ASSERT_EQ(0,  extracted.outwardIncidentWeight(2));
}

TEST_F(ASeparatedNodes, RestoresSavepoint) {
  SeparatedNodes nodes(10);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 3, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(3, 1), Edge(9, 1) };
  nodes.addNodes(new_nodes, new_edges);

  nodes.setSavepoint();

  vec<HypernodeID> communities { 0, 1, 1, 0, 2, kInvalidHypernode, 2, 0, 0, kInvalidHypernode };
  nodes.contract(communities, 3);

  nodes.restoreSavepoint();

  ASSERT_EQ(3,  nodes.numNodes());
  ASSERT_EQ(10,  nodes.numGraphNodes());
  ASSERT_EQ(3,  nodes.numEdges());
  ASSERT_EQ(3,  nodes.totalWeight());

  verifyIncidentEgdes(nodes.inwardEdges(0), { {0, 1}, {3, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(1), { {9, 1} });
  verifyIncidentEgdes(nodes.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, CreatesCopyFromSavepoint) {
  SeparatedNodes nodes(5);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 2, 1}, {2, 4, 1}};
  vec<Edge> new_edges { Edge(0, 1), Edge(4, 1), Edge(4, 1), Edge(3, 1) };
  nodes.addNodes(new_nodes, new_edges);

  nodes.setSavepoint();

  vec<HypernodeID> communities { 0, 1, 1, kInvalidHypernode, 2 };
  nodes.contract(communities, 4);

  SeparatedNodes other = nodes.createCopyFromSavepoint();

  ASSERT_EQ(3,  other.numNodes());
  ASSERT_EQ(5,  other.numGraphNodes());
  ASSERT_EQ(4,  other.numEdges());
  ASSERT_EQ(3,  other.totalWeight());

  verifyIncidentEgdes(other.inwardEdges(0), { {0, 1}, {4, 1} });
  verifyIncidentEgdes(other.inwardEdges(1), { {4, 1}, {3, 1} });
  verifyIncidentEgdes(other.inwardEdges(2), { });
}

TEST_F(ASeparatedNodes, ReinsertsSeparated) {
  SeparatedNodes nodes(3);
  vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> new_nodes { {0, 0, 1}, {1, 1, 1}, {2, 2, 2}};
  vec<Edge> new_edges { Edge(0, 1), Edge(2, 2), Edge(2, 1), Edge(1, 1) };
  nodes.addNodes(new_nodes, new_edges);

  vec<HypernodeWeight> weight = {1, 2, 3};
  Hypergraph original = HypergraphFactory::construct(3, 2, {{0, 1}, {1, 2}}, nullptr, weight.data());
  Hypergraph graph = HypergraphFactory::reinsertSeparatedNodes(original, nodes);

  ASSERT_EQ(6, graph.initialNumNodes());
  ASSERT_EQ(12, graph.initialNumEdges());
  ASSERT_EQ(10, graph.totalWeight());
  ASSERT_EQ(3, graph.numIncludedSeparated());

  verifyGraphNeighbors(graph, 0, {1, 3});
  verifyGraphNeighbors(graph, 1, {0, 2, 5});
  verifyGraphNeighbors(graph, 2, {1, 4, 5});
  verifyGraphNeighbors(graph, 3, {0});
  verifyGraphNeighbors(graph, 4, {2});
  verifyGraphNeighbors(graph, 5, {1, 2});
}

}
} // namespace mt_kahypar
