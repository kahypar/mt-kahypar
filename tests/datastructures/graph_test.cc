/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/tmp_graph.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using AGraph = HypergraphFixture<StaticHypergraph, StaticHypergraphFactory>;
using Graph = GraphT<StaticHypergraph>;
using Arc = typename Graph::Arc;
using ArcWeight = typename Graph::ArcWeight;

void verifyArcIterator(const Graph& graph,
                       const NodeID u,
                       const std::vector<NodeID>& arcs,
                       const std::vector<ArcWeight>& weights) {
  size_t size = 0;
  std::vector<bool> vis(arcs.size(), false);
  for ( const Arc& arc : graph.arcsOf(u) ) {
    for ( size_t pos = 0; pos < arcs.size(); ++pos ) {
      if ( arcs[pos] == arc.head ) {
        ASSERT_FALSE(vis[pos]);
        ASSERT_EQ(arcs[pos], arc.head);
        ASSERT_EQ(weights[pos], arc.weight);
        vis[pos] = true;
        ++size;
      }
    }
  }
  ASSERT_EQ(arcs.size(), size);
}

TEST_F(AGraph, HasCorrectNumNodesAndArcs) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(11, graph.numNodes());
  ASSERT_EQ(24, graph.numArcs());
}

TEST_F(AGraph, IteratesOverAllNodes) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<bool> vis(graph.numNodes(), false);
  for ( const NodeID& u : graph.nodes() ) {
    ASSERT_LE(u, graph.numNodes() - 1);
    vis[u] = true;
  }

  for ( size_t i = 0; i < vis.size(); ++i ) {
    ASSERT_TRUE(vis[i]);
  }
}

TEST_F(AGraph, VerifyTotalVolumeForUniformEdgeWeight) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(24, graph.totalVolume());
}

TEST_F(AGraph, VerifyNodeVolumeForUniformEdgeWeight) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(2, graph.nodeVolume(0));
  ASSERT_EQ(1, graph.nodeVolume(1));
  ASSERT_EQ(2, graph.nodeVolume(3));
  ASSERT_EQ(1, graph.nodeVolume(5));
  ASSERT_EQ(4, graph.nodeVolume(8));
  ASSERT_EQ(3, graph.nodeVolume(10));
}

TEST_F(AGraph, VerifyNodeVolumeForNonUniformEdgeWeight) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(0.75, graph.nodeVolume(0));
  ASSERT_EQ(0.25, graph.nodeVolume(1));
  ASSERT_EQ(0.25 + ( 1.0 / 3.0 ), graph.nodeVolume(3));
  ASSERT_EQ(1.0 / 3.0, graph.nodeVolume(5));
  ASSERT_EQ(1.0, graph.nodeVolume(8));
  ASSERT_EQ(1.0, graph.nodeVolume(10));
}

TEST_F(AGraph, WithCorrectVertexDegrees) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(2, graph.degree(0));
  ASSERT_EQ(1, graph.degree(1));
  ASSERT_EQ(2, graph.degree(2));
  ASSERT_EQ(2, graph.degree(3));
  ASSERT_EQ(2, graph.degree(4));
  ASSERT_EQ(1, graph.degree(5));
  ASSERT_EQ(2, graph.degree(6));
  ASSERT_EQ(2, graph.degree(7));
  ASSERT_EQ(4, graph.degree(8));
  ASSERT_EQ(3, graph.degree(9));
  ASSERT_EQ(3, graph.degree(10));
}

TEST_F(AGraph, HasCorrectAdjacentVertices1a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 0, {7, 8}, {1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices1b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 0, {7, 8}, {0.5, 0.25});
}

TEST_F(AGraph, HasCorrectAdjacentVertices1c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 0, {7, 8}, {1.0, 0.5});
}

TEST_F(AGraph, HasCorrectAdjacentVertices2a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 2, {7, 10}, {1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices2b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 2, {7, 10}, {0.5, 1.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices2c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 2, {7, 10}, {1.0, 2.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices3a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 5, {10}, {1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices3b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 5, {10}, {1.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices3c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 5, {10}, {1.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices4a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 6, {9, 10}, {1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices4b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 6, {9, 10}, {1.0 / 3.0, 1.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices4c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 6, {9, 10}, {2.0 / 3.0, 2.0 / 3.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices5a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 7, {0, 2}, {1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices5b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform,
    TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 7, {0, 2}, {0.5, 0.5});
}

TEST_F(AGraph, HasCorrectAdjacentVertices5c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 7, {0, 2}, {1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices6a) {
  Graph graph(hypergraph, LouvainEdgeWeight::uniform, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 8, {0, 1, 3, 4}, {1.0, 1.0, 1.0, 1.0});
}

TEST_F(AGraph, HasCorrectAdjacentVertices6b) {
  Graph graph(hypergraph, LouvainEdgeWeight::non_uniform, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 8, {0, 1, 3, 4}, {0.25, 0.25, 0.25, 0.25});
}

TEST_F(AGraph, HasCorrectAdjacentVertices6c) {
  Graph graph(hypergraph, LouvainEdgeWeight::degree, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyArcIterator(graph, 8, {0, 1, 3, 4}, {0.5, 0.25, 0.5, 0.5});
}

}
} // namespace mt_kahypar