/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>

#include "mt-kahypar/definitions.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph_factory.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace ds {

using ADynamicGraph = HypergraphFixture<DynamicGraph, DynamicGraphFactory, true>;

template<typename F, typename K>
void executeParallel(const F& f1, const K& f2) {
  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while ( cnt < 2 ) { }
    f1();
  }, [&] {
    ++cnt;
    while ( cnt < 2 ) { }
    f2();
  });
}

std::vector<HyperedgeID> expected_edges() {
  const HyperedgeID offset = DynamicAdjacencyArray::index_offset_per_node;
  std::vector<HyperedgeID> expected_iter;
  HyperedgeID current_offset = 0;
  // node 0
  current_offset += offset;
  // node 1
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 2
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 3
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 4
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 5
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  // node 6
  current_offset += offset;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  expected_iter.push_back(current_offset);
  current_offset += 1;
  return expected_iter;
}

TEST_F(ADynamicGraph, HasCorrectStats) {
  ASSERT_EQ(7,  hypergraph.initialNumNodes());
  ASSERT_EQ(12,  hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree());
  ASSERT_EQ(7,  hypergraph.totalWeight());
  ASSERT_EQ(2,  hypergraph.maxEdgeSize());
}

TEST_F(ADynamicGraph, HasCorrectInitialNodeIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_hn++, hn);
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(ADynamicGraph, HasCorrectNodeIteratorIfVerticesAreDisabled) {
  hypergraph.removeDegreeZeroHypernode(0);
  hypergraph.disableHypernode(5);
  const std::vector<HypernodeID> expected_iter =
    { 1, 2, 3, 4, 6 };
  HypernodeID pos = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_iter[pos++], hn);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicGraph, HasCorrectInitialEdgeIterator) {
  auto expected_iter = expected_edges();
  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicGraph, VerifiesIncidentEdges) {
  auto edges = expected_edges();
  verifyIncidentNets(0, { });
  verifyIncidentNets(1, { edges[0], edges[1] });
  verifyIncidentNets(2, { edges[2], edges[3] });
  verifyIncidentNets(3, { edges[4] });
  verifyIncidentNets(4, { edges[5], edges[6], edges[7] });
  verifyIncidentNets(5, { edges[8], edges[9] });
  verifyIncidentNets(6, { edges[10], edges[11] });
}

TEST_F(ADynamicGraph, VerifiesPinsOfEdges) {
  auto edges = expected_edges();
  verifyPins({ edges[0], edges[1], edges[3], edges[6], edges[7], edges[9] },
    { {1, 2}, {1, 4}, {2, 3}, {4, 5}, {4, 6}, {5, 6} });
}

TEST_F(ADynamicGraph, VerifiesVertexWeights) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(1, hypergraph.nodeWeight(hn));
  }
}

TEST_F(ADynamicGraph, HasCorrectEdgeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHypernode(5);
  hypergraph.disableHypernode(6);
  std::vector<HyperedgeID> expected_iter = expected_edges();
  expected_iter.resize(8);
  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicGraph, IteratesParallelOverAllNodes) {
  std::vector<uint8_t> visited(7, false);
  hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      visited[hn] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(ADynamicGraph, ModifiesNodeWeight) {
  hypergraph.setNodeWeight(0, 2);
  hypergraph.setNodeWeight(6, 2);
  ASSERT_EQ(2, hypergraph.nodeWeight(0));
  ASSERT_EQ(2, hypergraph.nodeWeight(6));
  hypergraph.updateTotalWeight(parallel_tag_t());
  ASSERT_EQ(9, hypergraph.totalWeight());
}


TEST_F(ADynamicGraph, VerifiesVertexDegrees) {
  ASSERT_EQ(0, hypergraph.nodeDegree(0));
  ASSERT_EQ(2, hypergraph.nodeDegree(1));
  ASSERT_EQ(2, hypergraph.nodeDegree(2));
  ASSERT_EQ(1, hypergraph.nodeDegree(3));
  ASSERT_EQ(3, hypergraph.nodeDegree(4));
  ASSERT_EQ(2, hypergraph.nodeDegree(5));
  ASSERT_EQ(2, hypergraph.nodeDegree(6));
}

TEST_F(ADynamicGraph, RemovesVertices) {
  hypergraph.removeHypernode(0);
  hypergraph.removeHypernode(5);
  ASSERT_EQ(2, hypergraph.numRemovedHypernodes());
}

TEST_F(ADynamicGraph, VerifiesEdgeWeights) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(1, hypergraph.edgeWeight(he));
  }
}

TEST_F(ADynamicGraph, ModifiesEdgeWeight) {
  hypergraph.setEdgeWeight(0, 2);
  hypergraph.setEdgeWeight(2, 2);
  ASSERT_EQ(2, hypergraph.edgeWeight(0));
  ASSERT_EQ(2, hypergraph.edgeWeight(2));
}

TEST_F(ADynamicGraph, VerifiesEdgeSizes) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(2, hypergraph.edgeSize(he));
  }
}

TEST_F(ADynamicGraph, SetsCommunityIDsForEachVertex) {
  hypergraph.setCommunityID(0, 1);
  hypergraph.setCommunityID(1, 1);
  hypergraph.setCommunityID(2, 1);
  hypergraph.setCommunityID(3, 2);
  hypergraph.setCommunityID(4, 2);
  hypergraph.setCommunityID(5, 3);
  hypergraph.setCommunityID(6, 3);

  ASSERT_EQ(1, hypergraph.communityID(0));
  ASSERT_EQ(1, hypergraph.communityID(1));
  ASSERT_EQ(1, hypergraph.communityID(2));
  ASSERT_EQ(2, hypergraph.communityID(3));
  ASSERT_EQ(2, hypergraph.communityID(4));
  ASSERT_EQ(3, hypergraph.communityID(5));
  ASSERT_EQ(3, hypergraph.communityID(6));
}

TEST_F(ADynamicGraph, ComparesStatsIfCopiedParallel) {
  DynamicGraph copy_hg = hypergraph.copy(parallel_tag_t());
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(ADynamicGraph, ComparesStatsIfCopiedSequential) {
  DynamicGraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(ADynamicGraph, ComparesIncidentEdgesIfCopiedParallel) {
  DynamicGraph copy_hg = hypergraph.copy(parallel_tag_t());
  auto edges = expected_edges();
  verifyIncidentNets(copy_hg, 0, { });
  verifyIncidentNets(copy_hg, 1, { edges[0], edges[1] });
  verifyIncidentNets(copy_hg, 2, { edges[2], edges[3] });
  verifyIncidentNets(copy_hg, 3, { edges[4] });
  verifyIncidentNets(copy_hg, 4, { edges[5], edges[6], edges[7] });
  verifyIncidentNets(copy_hg, 5, { edges[8], edges[9] });
  verifyIncidentNets(copy_hg, 6, { edges[10], edges[11] });
}

TEST_F(ADynamicGraph, ComparesIncidentEdgesIfCopiedSequential) {
  DynamicGraph copy_hg = hypergraph.copy();
  auto edges = expected_edges();
  verifyIncidentNets(copy_hg, 0, { });
  verifyIncidentNets(copy_hg, 1, { edges[0], edges[1] });
  verifyIncidentNets(copy_hg, 2, { edges[2], edges[3] });
  verifyIncidentNets(copy_hg, 3, { edges[4] });
  verifyIncidentNets(copy_hg, 4, { edges[5], edges[6], edges[7] });
  verifyIncidentNets(copy_hg, 5, { edges[8], edges[9] });
  verifyIncidentNets(copy_hg, 6, { edges[10], edges[11] });
}

} // namespace ds
} // namespace mt_kahypar
