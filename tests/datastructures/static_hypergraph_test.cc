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
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using AStaticHypergraph = HypergraphFixture<StaticHypergraph, StaticHypergraphFactory>;

TEST_F(AStaticHypergraph, HasCorrectStats) {
  ASSERT_EQ(7,  hypergraph.initialNumNodes());
  ASSERT_EQ(4,  hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree());
  ASSERT_EQ(7,  hypergraph.initialNumNodes(0));
  ASSERT_EQ(4,  hypergraph.initialNumEdges(0));
  ASSERT_EQ(12, hypergraph.initialNumPins(0));
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree(0));
  ASSERT_EQ(7,  hypergraph.totalWeight());
  ASSERT_EQ(4,  hypergraph.maxEdgeSize());
}

TEST_F(AStaticHypergraph, HasCorrectInitialNodeIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_hn++, hn);
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(AStaticHypergraph, HasCorrectNodeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHypernode(0);
  hypergraph.disableHypernode(5);
  const std::vector<HypernodeID> expected_iter =
    { 1, 2, 3, 4, 6 };
  HypernodeID pos = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_iter[pos++], hn);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, HasCorrectInitialEdgeIterator) {
  HyperedgeID expected_he = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_he++, he);
  }
  ASSERT_EQ(4, expected_he);
}

TEST_F(AStaticHypergraph, HasCorrectEdgeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHyperedge(0);
  hypergraph.disableHyperedge(2);
  const std::vector<HyperedgeID> expected_iter = { 1, 3 };
  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, IteratesParallelOverAllNodes) {
  std::vector<uint8_t> visited(7, false);
  hypergraph.doParallelForAllNodes(TBBNumaArena::GLOBAL_TASK_GROUP,
    [&](const HypernodeID hn) {
      visited[hn] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(AStaticHypergraph, IteratesParallelOverAllEdges) {
  std::vector<uint8_t> visited(4, false);
  hypergraph.doParallelForAllEdges(TBBNumaArena::GLOBAL_TASK_GROUP,
    [&](const HyperedgeID he) {
      visited[he] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(AStaticHypergraph, VerifiesIncidentNets1) {
  verifyIncidentNets(0, { 0, 1 });
}

TEST_F(AStaticHypergraph, VerifiesIncidentNets2) {
  verifyIncidentNets(2, { 0, 3 });
}

TEST_F(AStaticHypergraph, VerifiesIncidentNets3) {
  verifyIncidentNets(3, { 1, 2 });
}

TEST_F(AStaticHypergraph, VerifiesIncidentNets4) {
  verifyIncidentNets(6, { 2, 3 });
}

TEST_F(AStaticHypergraph, VerifiesPinsOfHyperedges) {
  verifyPins({ 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(AStaticHypergraph, VerifiesOriginalNodeIDs) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(hn, hypergraph.originalNodeID(hn));
  }
}

TEST_F(AStaticHypergraph, VerifiesVertexWeights) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(1, hypergraph.nodeWeight(hn));
  }
}

TEST_F(AStaticHypergraph, ModifiesNodeWeight) {
  hypergraph.setNodeWeight(0, 2);
  hypergraph.setNodeWeight(6, 2);
  ASSERT_EQ(2, hypergraph.nodeWeight(0));
  ASSERT_EQ(2, hypergraph.nodeWeight(6));
  hypergraph.updateTotalWeight(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(9, hypergraph.totalWeight());
}


TEST_F(AStaticHypergraph, VerifiesVertexDegrees) {
  ASSERT_EQ(2, hypergraph.nodeDegree(0));
  ASSERT_EQ(1, hypergraph.nodeDegree(1));
  ASSERT_EQ(2, hypergraph.nodeDegree(2));
  ASSERT_EQ(2, hypergraph.nodeDegree(3));
  ASSERT_EQ(2, hypergraph.nodeDegree(4));
  ASSERT_EQ(1, hypergraph.nodeDegree(5));
  ASSERT_EQ(2, hypergraph.nodeDegree(6));
}

TEST_F(AStaticHypergraph, RemovesVertices) {
  hypergraph.removeHypernode(0);
  hypergraph.removeHypernode(5);
  ASSERT_EQ(2, hypergraph.numRemovedHypernodes());
}

TEST_F(AStaticHypergraph, VerifiesOriginalEdgeIDs) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(he, hypergraph.originalEdgeID(he));
  }
}

TEST_F(AStaticHypergraph, VerifiesEdgeWeights) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(1, hypergraph.edgeWeight(he));
  }
}

TEST_F(AStaticHypergraph, ModifiesEdgeWeight) {
  hypergraph.setEdgeWeight(0, 2);
  hypergraph.setEdgeWeight(2, 2);
  ASSERT_EQ(2, hypergraph.edgeWeight(0));
  ASSERT_EQ(2, hypergraph.edgeWeight(2));
}

TEST_F(AStaticHypergraph, VerifiesEdgeSizes) {
  ASSERT_EQ(2, hypergraph.edgeSize(0));
  ASSERT_EQ(4, hypergraph.edgeSize(1));
  ASSERT_EQ(3, hypergraph.edgeSize(2));
  ASSERT_EQ(3, hypergraph.edgeSize(3));
}

TEST_F(AStaticHypergraph, SetsCommunityIDsForEachVertex) {
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

TEST_F(AStaticHypergraph, ComputesCorrectNumberOfCommunities) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunities());
}

TEST_F(AStaticHypergraph, ComputesCorrectNumberOfHypernodesInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunityHypernodes(0));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(1));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(2));
}

TEST_F(AStaticHypergraph, ComputesCorrectNumberOfPinsInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.numCommunityPins(0));
  ASSERT_EQ(4, hypergraph.numCommunityPins(1));
  ASSERT_EQ(3, hypergraph.numCommunityPins(2));
}

TEST_F(AStaticHypergraph, ComputesCorrectCommunityDegreeInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.communityDegree(0));
  ASSERT_EQ(4, hypergraph.communityDegree(1));
  ASSERT_EQ(3, hypergraph.communityDegree(2));
}

TEST_F(AStaticHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityZero) {
  assignCommunityIds();
  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(0), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(0)]);
  flag[hypergraph.communityNodeId(0)] = true;
  ASSERT_LE(hypergraph.communityNodeId(1), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(1)]);
  flag[hypergraph.communityNodeId(1)] = true;
  ASSERT_LE(hypergraph.communityNodeId(2), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(2)]);
  flag[hypergraph.communityNodeId(2)] = true;
}

TEST_F(AStaticHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityOne) {
  assignCommunityIds();
  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(3), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(3)]);
  flag[hypergraph.communityNodeId(3)] = true;
  ASSERT_LE(hypergraph.communityNodeId(4), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(4)]);
  flag[hypergraph.communityNodeId(4)] = true;
}

TEST_F(AStaticHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityTwo) {
  assignCommunityIds();
  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(5), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(5)]);
  flag[hypergraph.communityNodeId(5)] = true;
  ASSERT_LE(hypergraph.communityNodeId(6), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(6)]);
  flag[hypergraph.communityNodeId(6)] = true;
}

TEST_F(AStaticHypergraph, VerifiesNumberOfCommunitiesInHyperedges) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.numCommunitiesInHyperedge(0));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(1));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(2));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(3));
}

TEST_F(AStaticHypergraph, ChecksHyperedgeCommunityIterator1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(0) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, ChecksHyperedgeCommunityIterator2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 1 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(1) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, ChecksHyperedgeCommunityIterator3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 1, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(2) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, ChecksHyperedgeCommunityIterator4) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(3) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticHypergraph, VerifiesCommunityHyperedgeInternals1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(0, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(0, 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(1, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(3, 0));
  ASSERT_EQ(1, hypergraph.edgeSize(3, 0));
}

TEST_F(AStaticHypergraph, VerifiesCommunityHyperedgeInternals2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(1, 1));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 1));
  ASSERT_EQ(1, hypergraph.edgeWeight(2, 1));
  ASSERT_EQ(2, hypergraph.edgeSize(2, 1));
}

TEST_F(AStaticHypergraph, VerifiesCommunityHyperedgeInternals3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(2, 2));
  ASSERT_EQ(1, hypergraph.edgeSize(2, 2));
  ASSERT_EQ(1, hypergraph.edgeWeight(3, 2));
  ASSERT_EQ(2, hypergraph.edgeSize(3, 2));
}

TEST_F(AStaticHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
}

TEST_F(AStaticHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
}

TEST_F(AStaticHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

TEST_F(AStaticHypergraph, RemovesAHyperedgeFromTheHypergraph1) {
  hypergraph.removeEdge(0);
  verifyIncidentNets(0, { 1 });
  verifyIncidentNets(2, { 3 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(0, he);
  }
}

TEST_F(AStaticHypergraph, RemovesAHyperedgeFromTheHypergraph2) {
  hypergraph.removeEdge(1);
  verifyIncidentNets(0, { 0 });
  verifyIncidentNets(1, { });
  verifyIncidentNets(3, { 2 });
  verifyIncidentNets(4, { 2 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(1, he);
  }
}

TEST_F(AStaticHypergraph, RemovesAHyperedgeFromTheHypergraph3) {
  hypergraph.removeEdge(2);
  verifyIncidentNets(3, { 1 });
  verifyIncidentNets(4, { 1 });
  verifyIncidentNets(6, { 3 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(2, he);
  }
}

TEST_F(AStaticHypergraph, RemovesAHyperedgeFromTheHypergraph4) {
  hypergraph.removeEdge(3);
  verifyIncidentNets(2, { 0 });
  verifyIncidentNets(5, { });
  verifyIncidentNets(6, { 2 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(3, he);
  }
}

TEST_F(AStaticHypergraph, RestoresARemovedHyperedge1) {
  hypergraph.removeEdge(0);
  hypergraph.restoreEdge(0, 2);
  verifyIncidentNets(0, { 0, 1 });
  verifyIncidentNets(2, { 0, 3 });
  verifyPins({ 0 }, { {0, 2} });
}

TEST_F(AStaticHypergraph, RestoresARemovedHyperedge2) {
  hypergraph.removeEdge(1);
  hypergraph.restoreEdge(1, 4);
  verifyIncidentNets(0, { 0, 1 });
  verifyIncidentNets(1, { 1 });
  verifyIncidentNets(3, { 1, 2 });
  verifyIncidentNets(4, { 1, 2 });
  verifyPins({ 1 }, { {0, 1, 3, 4} });
}

TEST_F(AStaticHypergraph, RestoresARemovedHyperedge3) {
  hypergraph.removeEdge(2);
  hypergraph.restoreEdge(2, 3);
  verifyIncidentNets(3, { 1, 2 });
  verifyIncidentNets(4, { 1, 2 });
  verifyIncidentNets(6, { 2, 3 });
  verifyPins({ 2 }, { {3, 4, 6} });
}

TEST_F(AStaticHypergraph, RestoresARemovedHyperedge4) {
  hypergraph.removeEdge(3);
  hypergraph.restoreEdge(3, 3);
  verifyIncidentNets(2, { 0, 3 });
  verifyIncidentNets(5, { 3 });
  verifyIncidentNets(6, { 2, 3 });
  verifyPins({ 3 }, { {2, 5, 6} });
}

TEST_F(AStaticHypergraph, ComparesStatsIfCopiedParallel) {
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(AStaticHypergraph, ComparesStatsIfCopiedSequential) {
  StaticHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(AStaticHypergraph, ComparesIncidentNetsIfCopiedParallel) {
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyIncidentNets(copy_hg, 0, { 0, 1 });
  verifyIncidentNets(copy_hg, 1, { 1 });
  verifyIncidentNets(copy_hg, 2, { 0, 3 });
  verifyIncidentNets(copy_hg, 3, { 1, 2 });
  verifyIncidentNets(copy_hg, 4, { 1, 2 });
  verifyIncidentNets(copy_hg, 5, { 3 });
  verifyIncidentNets(copy_hg, 6, { 2, 3 });
}

TEST_F(AStaticHypergraph, ComparesIncidentNetsIfCopiedSequential) {
  StaticHypergraph copy_hg = hypergraph.copy();
  verifyIncidentNets(copy_hg, 0, { 0, 1 });
  verifyIncidentNets(copy_hg, 1, { 1 });
  verifyIncidentNets(copy_hg, 2, { 0, 3 });
  verifyIncidentNets(copy_hg, 3, { 1, 2 });
  verifyIncidentNets(copy_hg, 4, { 1, 2 });
  verifyIncidentNets(copy_hg, 5, { 3 });
  verifyIncidentNets(copy_hg, 6, { 2, 3 });
}

TEST_F(AStaticHypergraph, ComparesPinsOfHyperedgesIfCopiedParallel) {
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPins(copy_hg, { 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(AStaticHypergraph, ComparesPinsOfHyperedgesIfCopiedSequential) {
  StaticHypergraph copy_hg = hypergraph.copy();
  verifyPins(copy_hg, { 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(AStaticHypergraph, ComparesCommunityStatsIfCopiedParallel) {
  assignCommunityIds();
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.numCommunityHypernodes(0), copy_hg.numCommunityHypernodes(0));
  ASSERT_EQ(hypergraph.numCommunityHypernodes(1), copy_hg.numCommunityHypernodes(1));
  ASSERT_EQ(hypergraph.numCommunityHypernodes(2), copy_hg.numCommunityHypernodes(2));
  ASSERT_EQ(hypergraph.numCommunityPins(0), copy_hg.numCommunityPins(0));
  ASSERT_EQ(hypergraph.numCommunityPins(1), copy_hg.numCommunityPins(1));
  ASSERT_EQ(hypergraph.numCommunityPins(2), copy_hg.numCommunityPins(2));
  ASSERT_EQ(hypergraph.communityDegree(0), copy_hg.communityDegree(0));
  ASSERT_EQ(hypergraph.communityDegree(1), copy_hg.communityDegree(1));
  ASSERT_EQ(hypergraph.communityDegree(2), copy_hg.communityDegree(2));
}

TEST_F(AStaticHypergraph, ComparesCommunityStatsIfCopiedSequential) {
  assignCommunityIds();
  StaticHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.numCommunityHypernodes(0), copy_hg.numCommunityHypernodes(0));
  ASSERT_EQ(hypergraph.numCommunityHypernodes(1), copy_hg.numCommunityHypernodes(1));
  ASSERT_EQ(hypergraph.numCommunityHypernodes(2), copy_hg.numCommunityHypernodes(2));
  ASSERT_EQ(hypergraph.numCommunityPins(0), copy_hg.numCommunityPins(0));
  ASSERT_EQ(hypergraph.numCommunityPins(1), copy_hg.numCommunityPins(1));
  ASSERT_EQ(hypergraph.numCommunityPins(2), copy_hg.numCommunityPins(2));
  ASSERT_EQ(hypergraph.communityDegree(0), copy_hg.communityDegree(0));
  ASSERT_EQ(hypergraph.communityDegree(1), copy_hg.communityDegree(1));
  ASSERT_EQ(hypergraph.communityDegree(2), copy_hg.communityDegree(2));
}

TEST_F(AStaticHypergraph, ComparesCommunityIdsIfCopiedParallel) {
  assignCommunityIds();
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.communityID(0), copy_hg.communityID(0));
  ASSERT_EQ(hypergraph.communityID(1), copy_hg.communityID(1));
  ASSERT_EQ(hypergraph.communityID(2), copy_hg.communityID(2));
  ASSERT_EQ(hypergraph.communityID(3), copy_hg.communityID(3));
  ASSERT_EQ(hypergraph.communityID(4), copy_hg.communityID(4));
  ASSERT_EQ(hypergraph.communityID(5), copy_hg.communityID(5));
  ASSERT_EQ(hypergraph.communityID(6), copy_hg.communityID(6));
}

TEST_F(AStaticHypergraph, ComparesCommunityIdsIfCopiedSequential) {
  assignCommunityIds();
  StaticHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.communityID(0), copy_hg.communityID(0));
  ASSERT_EQ(hypergraph.communityID(1), copy_hg.communityID(1));
  ASSERT_EQ(hypergraph.communityID(2), copy_hg.communityID(2));
  ASSERT_EQ(hypergraph.communityID(3), copy_hg.communityID(3));
  ASSERT_EQ(hypergraph.communityID(4), copy_hg.communityID(4));
  ASSERT_EQ(hypergraph.communityID(5), copy_hg.communityID(5));
  ASSERT_EQ(hypergraph.communityID(6), copy_hg.communityID(6));
}

TEST_F(AStaticHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(0), copy_hg.numCommunitiesInHyperedge(0));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(1), copy_hg.numCommunitiesInHyperedge(1));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(2), copy_hg.numCommunitiesInHyperedge(2));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(3), copy_hg.numCommunitiesInHyperedge(3));
}

TEST_F(AStaticHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  StaticHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(0), copy_hg.numCommunitiesInHyperedge(0));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(1), copy_hg.numCommunitiesInHyperedge(1));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(2), copy_hg.numCommunitiesInHyperedge(2));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(3), copy_hg.numCommunitiesInHyperedge(3));
}

TEST_F(AStaticHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  StaticHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

TEST_F(AStaticHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  StaticHypergraph copy_hg = hypergraph.copy();
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

TEST_F(AStaticHypergraph, ContractsCommunities1) {
  parallel::scalable_vector<HypernodeID> c_mapping = {1, 4, 1, 5, 5, 4, 5};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Mapping
  ASSERT_EQ(0, c_mapping[0]);
  ASSERT_EQ(1, c_mapping[1]);
  ASSERT_EQ(0, c_mapping[2]);
  ASSERT_EQ(2, c_mapping[3]);
  ASSERT_EQ(2, c_mapping[4]);
  ASSERT_EQ(1, c_mapping[5]);
  ASSERT_EQ(2, c_mapping[6]);

  // Verify Stats
  ASSERT_EQ(3, c_hypergraph.initialNumNodes());
  ASSERT_EQ(1, c_hypergraph.initialNumEdges());
  ASSERT_EQ(3, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(3, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(0));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(1));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(2));

  // Verify Hyperedge Weights
  ASSERT_EQ(2, c_hypergraph.edgeWeight(0));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, 0, { 0 });
  verifyIncidentNets(c_hypergraph, 1, { 0 });
  verifyIncidentNets(c_hypergraph, 2, { 0 });
  verifyPins(c_hypergraph, { 0 }, { {0, 1, 2} });
}

TEST_F(AStaticHypergraph, ContractsCommunities2) {
  parallel::scalable_vector<HypernodeID> c_mapping = {1, 4, 1, 5, 5, 6, 5};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Mapping
  ASSERT_EQ(0, c_mapping[0]);
  ASSERT_EQ(1, c_mapping[1]);
  ASSERT_EQ(0, c_mapping[2]);
  ASSERT_EQ(2, c_mapping[3]);
  ASSERT_EQ(2, c_mapping[4]);
  ASSERT_EQ(3, c_mapping[5]);
  ASSERT_EQ(2, c_mapping[6]);

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(2, c_hypergraph.initialNumEdges());
  ASSERT_EQ(6, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(3, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(0));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(1));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(2));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(3));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(0));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(1));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, 0, { 0, 1 });
  verifyIncidentNets(c_hypergraph, 1, { 0 });
  verifyIncidentNets(c_hypergraph, 2, { 0, 1 });
  verifyIncidentNets(c_hypergraph, 3, { 1 });
  verifyPins(c_hypergraph, { 0, 1 }, { {0, 1, 2}, {0, 2, 3} });
}

TEST_F(AStaticHypergraph, ContractsCommunities3) {
  parallel::scalable_vector<HypernodeID> c_mapping = {2, 2, 0, 5, 5, 1, 1};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Mapping
  ASSERT_EQ(2, c_mapping[0]);
  ASSERT_EQ(2, c_mapping[1]);
  ASSERT_EQ(0, c_mapping[2]);
  ASSERT_EQ(3, c_mapping[3]);
  ASSERT_EQ(3, c_mapping[4]);
  ASSERT_EQ(1, c_mapping[5]);
  ASSERT_EQ(1, c_mapping[6]);

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(4, c_hypergraph.initialNumEdges());
  ASSERT_EQ(8, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(1, c_hypergraph.nodeWeight(0));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(1));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(2));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(3));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(0));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(1));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(2));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(3));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, 0, { 0, 3 });
  verifyIncidentNets(c_hypergraph, 1, { 2, 3 });
  verifyIncidentNets(c_hypergraph, 2, { 0, 1 });
  verifyIncidentNets(c_hypergraph, 3, { 1, 2 });
  verifyPins(c_hypergraph, { 0, 1, 2, 3 },
    { {0, 2}, {2, 3}, {1, 3}, {0, 1} });
}

TEST_F(AStaticHypergraph, ContractsCommunitiesWithDisabledHypernodes) {
  hypergraph.disableHypernode(0);
  hypergraph.disableHypernode(6);

  parallel::scalable_vector<HypernodeID> c_mapping = {0, 1, 1, 2, 2, 2, 6};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Mapping
  ASSERT_EQ(0, c_mapping[1]);
  ASSERT_EQ(0, c_mapping[2]);
  ASSERT_EQ(1, c_mapping[3]);
  ASSERT_EQ(1, c_mapping[4]);
  ASSERT_EQ(1, c_mapping[5]);

  // Verify Stats
  ASSERT_EQ(2, c_hypergraph.initialNumNodes());
  ASSERT_EQ(1, c_hypergraph.initialNumEdges());
  ASSERT_EQ(2, c_hypergraph.initialNumPins());
  ASSERT_EQ(5, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(0));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(1));

  // Verify Hyperedge Weights
  ASSERT_EQ(2, c_hypergraph.edgeWeight(0));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, 0, { 0 });
  verifyIncidentNets(c_hypergraph, 1, { 0 });
  verifyPins(c_hypergraph, { 0 }, { {0, 1} });
}

TEST_F(AStaticHypergraph, ContractsCommunitiesWithDisabledHyperedges) {
  hypergraph.disableHyperedge(3);

  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 0, 1, 1, 2, 3};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Mapping
  ASSERT_EQ(0, c_mapping[0]);
  ASSERT_EQ(0, c_mapping[1]);
  ASSERT_EQ(0, c_mapping[2]);
  ASSERT_EQ(1, c_mapping[3]);
  ASSERT_EQ(1, c_mapping[4]);
  ASSERT_EQ(2, c_mapping[5]);
  ASSERT_EQ(3, c_mapping[6]);

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(2, c_hypergraph.initialNumEdges());
  ASSERT_EQ(4, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(3, c_hypergraph.nodeWeight(0));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(1));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(2));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(3));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(0));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(1));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, 0, { 0 });
  verifyIncidentNets(c_hypergraph, 1, { 0, 1 });
  verifyIncidentNets(c_hypergraph, 2, { });
  verifyIncidentNets(c_hypergraph, 3, { 1 });
  verifyPins(c_hypergraph, { 0, 1 },
    { {0, 1}, {1, 3} });
}

TEST_F(AStaticHypergraph, ContractCommunitiesIfCommunityInformationAreAvailable) {
  assignCommunityIds();
  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 1, 2, 2, 3, 3};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Community Ids
  ASSERT_EQ(0, c_hypergraph.communityID(0));
  ASSERT_EQ(0, c_hypergraph.communityID(1));
  ASSERT_EQ(1, c_hypergraph.communityID(2));
  ASSERT_EQ(2, c_hypergraph.communityID(3));

  // Verify Community Stats
  ASSERT_EQ(2, c_hypergraph.numCommunityHypernodes(0));
  ASSERT_EQ(1, c_hypergraph.numCommunityHypernodes(1));
  ASSERT_EQ(1, c_hypergraph.numCommunityHypernodes(2));
  ASSERT_EQ(4, c_hypergraph.numCommunityPins(0));
  ASSERT_EQ(2, c_hypergraph.numCommunityPins(1));
  ASSERT_EQ(2, c_hypergraph.numCommunityPins(2));
  ASSERT_EQ(4, c_hypergraph.communityDegree(0));
  ASSERT_EQ(2, c_hypergraph.communityDegree(1));
  ASSERT_EQ(2, c_hypergraph.communityDegree(2));
}

TEST_F(AStaticHypergraph, ContractCommunitiesIfCommunityHyperedgesAreAvailable) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 1, 2, 2, 3, 3};
  StaticHypergraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);

  // Verify Hypergraph Structure
  verifyCommunityPins(c_hypergraph, 0, { 0, 1, 3 },
    { {0, 1}, {0}, {1} });
  verifyCommunityPins(c_hypergraph, 1, { 1, 2 },
    { {2}, {2} });
  verifyCommunityPins(c_hypergraph, 2, { 2, 3 },
    { {3}, {3} });
}

}
} // namespace mt_kahypar