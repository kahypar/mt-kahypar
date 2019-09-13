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

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using AHypergraphWithTwoStreamingHypergraphs = AHypergraph<2>;
using TestHypergraph = typename AHypergraphWithTwoStreamingHypergraphs::TestHypergraph;

TestHypergraph construct_test_hypergraph(const AHypergraphWithTwoStreamingHypergraphs& test) {
  return test.construct_hypergraph(7, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
                                      { 0, 0, 0, 1, 1, 1, 1 },
                                      { 0, 0, 1, 1 } );
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContainsCorrectNumberofNodesEdgesAndPins) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  ASSERT_EQ(7, hypergraph.initialNumNodes());
  ASSERT_EQ(4, hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIterators) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  std::set<HypernodeID> node_0 = {0, 1, 2};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(0) ) {
    ASSERT_TRUE(node_0.find(hn) != node_0.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_0.size());

  std::set<HypernodeID> node_1 = {281474976710656, 281474976710657, 281474976710658, 281474976710659};
  num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(1) ) {
    ASSERT_TRUE(node_1.find(hn) != node_1.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_1.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(1);

  std::set<HypernodeID> node_0 = {0, 2};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(0) ) {
    ASSERT_TRUE(node_0.find(hn) != node_0.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_0.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710658);
  
  std::set<HypernodeID> node_1 = {281474976710656, 281474976710659};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(1) ) {
    ASSERT_TRUE(node_1.find(hn) != node_1.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_1.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksOriginalNodeIds) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  std::set<HypernodeID> node_0 = {0, 1, 2};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(0) ) {
    ASSERT_TRUE(node_0.find(hypergraph.originalNodeID(hn)) != node_0.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_0.size());

  std::set<HypernodeID> node_1 = {3, 4, 5, 6};
  num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(1) ) {
    ASSERT_TRUE(node_1.find(hypergraph.originalNodeID(hn)) != node_1.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, node_1.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIterators) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  std::set<HyperedgeID> edge_0 = {0, 1};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges(0) ) {
    ASSERT_TRUE(edge_0.find(he) != edge_0.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edge_0.size());

  std::set<HyperedgeID> edge_1 = {281474976710656, 281474976710657};
  num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges(1) ) {
    ASSERT_TRUE(edge_1.find(he) != edge_1.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edge_1.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);
  
  std::set<HyperedgeID> edge_0 = {0};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges(0) ) {
    ASSERT_TRUE(edge_0.find(he) != edge_0.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edge_0.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(281474976710656);
  
  std::set<HyperedgeID> edge_1 = {281474976710657};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges(1) ) {
    ASSERT_TRUE(edge_1.find(he) != edge_1.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edge_1.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIterator) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  std::set<HypernodeID> nodes = {0, 1, 2, 281474976710656, 281474976710657, 281474976710658, 281474976710659};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_TRUE(nodes.find(hn) != nodes.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, nodes.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710659);
  
  std::set<HypernodeID> nodes = {0, 1, 281474976710656, 281474976710658};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_TRUE(nodes.find(hn) != nodes.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, nodes.size());
}
TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(0);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710656);
  hypergraph.disableHypernode(281474976710659);
  
  std::set<HypernodeID> nodes = {1, 281474976710657, 281474976710658};
  size_t num_vertices = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_TRUE(nodes.find(hn) != nodes.end()) << V(hn);
    num_vertices++;
  }
  ASSERT_EQ(num_vertices, nodes.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIterators) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  
  std::set<HyperedgeID> edges = {0, 1, 281474976710656, 281474976710657};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_TRUE(edges.find(he) != edges.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edges.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);
  hypergraph.disableHyperedge(281474976710656);
  
  std::set<HyperedgeID> edges = {0, 281474976710657};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_TRUE(edges.find(he) != edges.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edges.size());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(0);
  hypergraph.disableHyperedge(281474976710657);
  
  std::set<HyperedgeID> edges = {1, 281474976710656};
  size_t num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_TRUE(edges.find(he) != edges.end()) << V(he);
    num_edges++;
  }
  ASSERT_EQ(num_edges, edges.size());
}


} // namespace ds
} // namespace mt_kahypar
