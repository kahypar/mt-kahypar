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

#include "mt-kahypar/definitions.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"

namespace mt_kahypar {
namespace ds {

using ADynamicHypergraph = HypergraphFixture<DynamicHypergraph, DynamicHypergraphFactory>;

TEST_F(ADynamicHypergraph, HasCorrectStats) {
  ASSERT_EQ(7,  hypergraph.initialNumNodes());
  ASSERT_EQ(4,  hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree());
  ASSERT_EQ(7,  hypergraph.totalWeight());
  ASSERT_EQ(4,  hypergraph.maxEdgeSize());
}

TEST_F(ADynamicHypergraph, HasCorrectInitialNodeIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_hn++, hn);
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(ADynamicHypergraph, HasCorrectNodeIteratorIfVerticesAreDisabled) {
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

TEST_F(ADynamicHypergraph, HasCorrectInitialEdgeIterator) {
  HyperedgeID expected_he = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_he++, he);
  }
  ASSERT_EQ(4, expected_he);
}

TEST_F(ADynamicHypergraph, HasCorrectEdgeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHyperedge(0);
  hypergraph.disableHyperedge(2);
  const std::vector<HyperedgeID> expected_iter = { 1, 3 };
  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicHypergraph, IteratesParallelOverAllNodes) {
  std::vector<uint8_t> visited(7, false);
  hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      visited[hn] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(ADynamicHypergraph, IteratesParallelOverAllEdges) {
  std::vector<uint8_t> visited(4, false);
  hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
      visited[he] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(ADynamicHypergraph, VerifiesIncidentNets1) {
  verifyIncidentNets(0, { 0, 1 });
}

TEST_F(ADynamicHypergraph, VerifiesIncidentNets2) {
  verifyIncidentNets(2, { 0, 3 });
}

TEST_F(ADynamicHypergraph, VerifiesIncidentNets3) {
  verifyIncidentNets(3, { 1, 2 });
}

TEST_F(ADynamicHypergraph, VerifiesIncidentNets4) {
  verifyIncidentNets(6, { 2, 3 });
}

TEST_F(ADynamicHypergraph, VerifiesPinsOfHyperedges) {
  verifyPins({ 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, VerifiesVertexWeights) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(1, hypergraph.nodeWeight(hn));
  }
}

TEST_F(ADynamicHypergraph, ModifiesNodeWeight) {
  hypergraph.setNodeWeight(0, 2);
  hypergraph.setNodeWeight(6, 2);
  ASSERT_EQ(2, hypergraph.nodeWeight(0));
  ASSERT_EQ(2, hypergraph.nodeWeight(6));
  hypergraph.updateTotalWeight(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(9, hypergraph.totalWeight());
}


TEST_F(ADynamicHypergraph, VerifiesVertexDegrees) {
  ASSERT_EQ(2, hypergraph.nodeDegree(0));
  ASSERT_EQ(1, hypergraph.nodeDegree(1));
  ASSERT_EQ(2, hypergraph.nodeDegree(2));
  ASSERT_EQ(2, hypergraph.nodeDegree(3));
  ASSERT_EQ(2, hypergraph.nodeDegree(4));
  ASSERT_EQ(1, hypergraph.nodeDegree(5));
  ASSERT_EQ(2, hypergraph.nodeDegree(6));
}

TEST_F(ADynamicHypergraph, RemovesVertices) {
  hypergraph.removeHypernode(0);
  hypergraph.removeHypernode(5);
  ASSERT_EQ(2, hypergraph.numRemovedHypernodes());
}

TEST_F(ADynamicHypergraph, VerifiesEdgeWeights) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(1, hypergraph.edgeWeight(he));
  }
}

TEST_F(ADynamicHypergraph, ModifiesEdgeWeight) {
  hypergraph.setEdgeWeight(0, 2);
  hypergraph.setEdgeWeight(2, 2);
  ASSERT_EQ(2, hypergraph.edgeWeight(0));
  ASSERT_EQ(2, hypergraph.edgeWeight(2));
}

TEST_F(ADynamicHypergraph, VerifiesEdgeSizes) {
  ASSERT_EQ(2, hypergraph.edgeSize(0));
  ASSERT_EQ(4, hypergraph.edgeSize(1));
  ASSERT_EQ(3, hypergraph.edgeSize(2));
  ASSERT_EQ(3, hypergraph.edgeSize(3));
}

TEST_F(ADynamicHypergraph, SetsCommunityIDsForEachVertex) {
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


TEST_F(ADynamicHypergraph, ComputesCorrectNumberOfCommunities) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunities());
}

TEST_F(ADynamicHypergraph, ComputesCorrectNumberOfHypernodesInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunityHypernodes(0));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(1));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(2));
}

TEST_F(ADynamicHypergraph, ComputesCorrectNumberOfPinsInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.numCommunityPins(0));
  ASSERT_EQ(4, hypergraph.numCommunityPins(1));
  ASSERT_EQ(3, hypergraph.numCommunityPins(2));
}

TEST_F(ADynamicHypergraph, ComputesCorrectCommunityDegreeInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.communityDegree(0));
  ASSERT_EQ(4, hypergraph.communityDegree(1));
  ASSERT_EQ(3, hypergraph.communityDegree(2));
}

TEST_F(ADynamicHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityZero) {
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

TEST_F(ADynamicHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityOne) {
  assignCommunityIds();
  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(3), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(3)]);
  flag[hypergraph.communityNodeId(3)] = true;
  ASSERT_LE(hypergraph.communityNodeId(4), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(4)]);
  flag[hypergraph.communityNodeId(4)] = true;
}

TEST_F(ADynamicHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityTwo) {
  assignCommunityIds();
  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(5), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(5)]);
  flag[hypergraph.communityNodeId(5)] = true;
  ASSERT_LE(hypergraph.communityNodeId(6), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(6)]);
  flag[hypergraph.communityNodeId(6)] = true;
}

TEST_F(ADynamicHypergraph, VerifiesNumberOfCommunitiesInHyperedges) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.numCommunitiesInHyperedge(0));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(1));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(2));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(3));
}

TEST_F(ADynamicHypergraph, ChecksHyperedgeCommunityIterator1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(0) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicHypergraph, ChecksHyperedgeCommunityIterator2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 1 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(1) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicHypergraph, ChecksHyperedgeCommunityIterator3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 1, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(2) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicHypergraph, ChecksHyperedgeCommunityIterator4) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(3) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(ADynamicHypergraph, VerifiesCommunityHyperedgeInternals1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(0, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(0, 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(1, 0));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(3, 0));
  ASSERT_EQ(1, hypergraph.edgeSize(3, 0));
}

TEST_F(ADynamicHypergraph, VerifiesCommunityHyperedgeInternals2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(1, 1));
  ASSERT_EQ(2, hypergraph.edgeSize(1, 1));
  ASSERT_EQ(1, hypergraph.edgeWeight(2, 1));
  ASSERT_EQ(2, hypergraph.edgeSize(2, 1));
}

TEST_F(ADynamicHypergraph, VerifiesCommunityHyperedgeInternals3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(2, 2));
  ASSERT_EQ(1, hypergraph.edgeSize(2, 2));
  ASSERT_EQ(1, hypergraph.edgeWeight(3, 2));
  ASSERT_EQ(2, hypergraph.edgeSize(3, 2));
}

TEST_F(ADynamicHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
}

TEST_F(ADynamicHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
}

TEST_F(ADynamicHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

TEST_F(ADynamicHypergraph, RemovesAHyperedgeFromTheHypergraph1) {
  hypergraph.removeEdge(0);
  verifyIncidentNets(0, { 1 });
  verifyIncidentNets(2, { 3 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(0, he);
  }
}

TEST_F(ADynamicHypergraph, RemovesAHyperedgeFromTheHypergraph2) {
  hypergraph.removeEdge(1);
  verifyIncidentNets(0, { 0 });
  verifyIncidentNets(1, { });
  verifyIncidentNets(3, { 2 });
  verifyIncidentNets(4, { 2 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(1, he);
  }
}

TEST_F(ADynamicHypergraph, RemovesAHyperedgeFromTheHypergraph3) {
  hypergraph.removeEdge(2);
  verifyIncidentNets(3, { 1 });
  verifyIncidentNets(4, { 1 });
  verifyIncidentNets(6, { 3 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(2, he);
  }
}

TEST_F(ADynamicHypergraph, RemovesAHyperedgeFromTheHypergraph4) {
  hypergraph.removeEdge(3);
  verifyIncidentNets(2, { 0 });
  verifyIncidentNets(5, { });
  verifyIncidentNets(6, { 2 });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(3, he);
  }
}

TEST_F(ADynamicHypergraph, ComparesStatsIfCopiedParallel) {
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(ADynamicHypergraph, ComparesStatsIfCopiedSequential) {
  DynamicHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(ADynamicHypergraph, ComparesIncidentNetsIfCopiedParallel) {
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyIncidentNets(copy_hg, 0, { 0, 1 });
  verifyIncidentNets(copy_hg, 1, { 1 });
  verifyIncidentNets(copy_hg, 2, { 0, 3 });
  verifyIncidentNets(copy_hg, 3, { 1, 2 });
  verifyIncidentNets(copy_hg, 4, { 1, 2 });
  verifyIncidentNets(copy_hg, 5, { 3 });
  verifyIncidentNets(copy_hg, 6, { 2, 3 });
}

TEST_F(ADynamicHypergraph, ComparesIncidentNetsIfCopiedSequential) {
  DynamicHypergraph copy_hg = hypergraph.copy();
  verifyIncidentNets(copy_hg, 0, { 0, 1 });
  verifyIncidentNets(copy_hg, 1, { 1 });
  verifyIncidentNets(copy_hg, 2, { 0, 3 });
  verifyIncidentNets(copy_hg, 3, { 1, 2 });
  verifyIncidentNets(copy_hg, 4, { 1, 2 });
  verifyIncidentNets(copy_hg, 5, { 3 });
  verifyIncidentNets(copy_hg, 6, { 2, 3 });
}

TEST_F(ADynamicHypergraph, ComparesPinsOfHyperedgesIfCopiedParallel) {
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyPins(copy_hg, { 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, ComparesPinsOfHyperedgesIfCopiedSequential) {
  DynamicHypergraph copy_hg = hypergraph.copy();
  verifyPins(copy_hg, { 0, 1, 2, 3 },
    { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, ComparesCommunityStatsIfCopiedParallel) {
  assignCommunityIds();
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
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

TEST_F(ADynamicHypergraph, ComparesCommunityStatsIfCopiedSequential) {
  assignCommunityIds();
  DynamicHypergraph copy_hg = hypergraph.copy();
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

TEST_F(ADynamicHypergraph, ComparesCommunityIdsIfCopiedParallel) {
  assignCommunityIds();
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.communityID(0), copy_hg.communityID(0));
  ASSERT_EQ(hypergraph.communityID(1), copy_hg.communityID(1));
  ASSERT_EQ(hypergraph.communityID(2), copy_hg.communityID(2));
  ASSERT_EQ(hypergraph.communityID(3), copy_hg.communityID(3));
  ASSERT_EQ(hypergraph.communityID(4), copy_hg.communityID(4));
  ASSERT_EQ(hypergraph.communityID(5), copy_hg.communityID(5));
  ASSERT_EQ(hypergraph.communityID(6), copy_hg.communityID(6));
}

TEST_F(ADynamicHypergraph, ComparesCommunityIdsIfCopiedSequential) {
  assignCommunityIds();
  DynamicHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.communityID(0), copy_hg.communityID(0));
  ASSERT_EQ(hypergraph.communityID(1), copy_hg.communityID(1));
  ASSERT_EQ(hypergraph.communityID(2), copy_hg.communityID(2));
  ASSERT_EQ(hypergraph.communityID(3), copy_hg.communityID(3));
  ASSERT_EQ(hypergraph.communityID(4), copy_hg.communityID(4));
  ASSERT_EQ(hypergraph.communityID(5), copy_hg.communityID(5));
  ASSERT_EQ(hypergraph.communityID(6), copy_hg.communityID(6));
}

TEST_F(ADynamicHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(0), copy_hg.numCommunitiesInHyperedge(0));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(1), copy_hg.numCommunitiesInHyperedge(1));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(2), copy_hg.numCommunitiesInHyperedge(2));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(3), copy_hg.numCommunitiesInHyperedge(3));
}

TEST_F(ADynamicHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  DynamicHypergraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(0), copy_hg.numCommunitiesInHyperedge(0));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(1), copy_hg.numCommunitiesInHyperedge(1));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(2), copy_hg.numCommunitiesInHyperedge(2));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(3), copy_hg.numCommunitiesInHyperedge(3));
}

TEST_F(ADynamicHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  DynamicHypergraph copy_hg = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

TEST_F(ADynamicHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  DynamicHypergraph copy_hg = hypergraph.copy();
  verifyCommunityPins(0, { 0, 1, 3 },
    { {0, 2}, {0, 1}, {2} });
  verifyCommunityPins(1, { 1, 2 },
    { {3, 4}, {3, 4} });
  verifyCommunityPins(2, { 2, 3 },
    { {6}, {5, 6} });
}

} // namespace ds
} // namespace mt_kahypar