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

#include <atomic>

#include "mt-kahypar/definitions.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"

namespace mt_kahypar {
namespace ds {

using ADynamicHypergraph = HypergraphFixture<DynamicHypergraph, DynamicHypergraphFactory>;

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

TEST_F(ADynamicHypergraph, RegistersAContraction1) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_EQ(1, hypergraph.contractionTree(0));
  ASSERT_EQ(1, hypergraph.contractionTree(1));
  ASSERT_EQ(0, hypergraph.referenceCount(0));
  ASSERT_EQ(1, hypergraph.referenceCount(1));
}

TEST_F(ADynamicHypergraph, RegistersAContraction2) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_EQ(4, hypergraph.contractionTree(3));
  ASSERT_EQ(4, hypergraph.contractionTree(4));
  ASSERT_EQ(0, hypergraph.referenceCount(3));
  ASSERT_EQ(1, hypergraph.referenceCount(4));
}

TEST_F(ADynamicHypergraph, RegistersAContraction3) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_EQ(4, hypergraph.contractionTree(2));
  ASSERT_EQ(4, hypergraph.contractionTree(3));
  ASSERT_EQ(4, hypergraph.contractionTree(4));
  ASSERT_EQ(0, hypergraph.referenceCount(2));
  ASSERT_EQ(0, hypergraph.referenceCount(3));
  ASSERT_EQ(2, hypergraph.referenceCount(4));
}

TEST_F(ADynamicHypergraph, RegistersAContraction4) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_TRUE(hypergraph.registerContraction(6, 4));
  ASSERT_EQ(4, hypergraph.contractionTree(2));
  ASSERT_EQ(4, hypergraph.contractionTree(3));
  ASSERT_EQ(6, hypergraph.contractionTree(4));
  ASSERT_EQ(6, hypergraph.contractionTree(6));
  ASSERT_EQ(0, hypergraph.referenceCount(2));
  ASSERT_EQ(0, hypergraph.referenceCount(3));
  ASSERT_EQ(2, hypergraph.referenceCount(4));
  ASSERT_EQ(1, hypergraph.referenceCount(6));
}

TEST_F(ADynamicHypergraph, RegistersAContraction5) {
  ASSERT_TRUE(hypergraph.registerContraction(6, 4));
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_EQ(6, hypergraph.contractionTree(2));
  ASSERT_EQ(6, hypergraph.contractionTree(3));
  ASSERT_EQ(6, hypergraph.contractionTree(4));
  ASSERT_EQ(6, hypergraph.contractionTree(6));
  ASSERT_EQ(0, hypergraph.referenceCount(2));
  ASSERT_EQ(0, hypergraph.referenceCount(3));
  ASSERT_EQ(0, hypergraph.referenceCount(4));
  ASSERT_EQ(3, hypergraph.referenceCount(6));
}

TEST_F(ADynamicHypergraph, RegistersAContraction6) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(6, 4));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_EQ(4, hypergraph.contractionTree(2));
  ASSERT_EQ(4, hypergraph.contractionTree(3));
  ASSERT_EQ(6, hypergraph.contractionTree(4));
  ASSERT_EQ(6, hypergraph.contractionTree(6));
  ASSERT_EQ(0, hypergraph.referenceCount(2));
  ASSERT_EQ(0, hypergraph.referenceCount(3));
  ASSERT_EQ(2, hypergraph.referenceCount(4));
  ASSERT_EQ(1, hypergraph.referenceCount(6));
}

TEST_F(ADynamicHypergraph, RegistersAContraction7) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  ASSERT_EQ(1, hypergraph.contractionTree(0));
  ASSERT_EQ(2, hypergraph.contractionTree(1));
  ASSERT_EQ(3, hypergraph.contractionTree(2));
  ASSERT_EQ(3, hypergraph.contractionTree(3));
  ASSERT_EQ(0, hypergraph.referenceCount(0));
  ASSERT_EQ(1, hypergraph.referenceCount(1));
  ASSERT_EQ(1, hypergraph.referenceCount(2));
  ASSERT_EQ(1, hypergraph.referenceCount(3));
}

TEST_F(ADynamicHypergraph, RegistersAContractionThatInducesACycle1) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_FALSE(hypergraph.registerContraction(0, 1));
}

TEST_F(ADynamicHypergraph, RegistersAContractionThatInducesACycle2) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  hypergraph.decrementReferenceCount(1);
  ASSERT_FALSE(hypergraph.registerContraction(0, 2));
}

TEST_F(ADynamicHypergraph, RegistersAContractionThatInducesACycle3) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  hypergraph.decrementReferenceCount(1);
  hypergraph.decrementReferenceCount(2);
  ASSERT_FALSE(hypergraph.registerContraction(0, 3));
}

TEST_F(ADynamicHypergraph, RegistersAContractionThatInducesACycle4) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_TRUE(hypergraph.registerContraction(6, 4));
  ASSERT_FALSE(hypergraph.registerContraction(2, 6));
}

TEST_F(ADynamicHypergraph, RegistersAContractionThatInducesACycle5) {
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));
  ASSERT_TRUE(hypergraph.registerContraction(4, 2));
  ASSERT_TRUE(hypergraph.registerContraction(6, 4));
  ASSERT_TRUE(hypergraph.registerContraction(5, 6));
  ASSERT_FALSE(hypergraph.registerContraction(4, 5));
}

TEST_F(ADynamicHypergraph, RegisterContractionsInParallel1) {
  executeParallel([&] {
    ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  }, [&] {
    ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  });

  ASSERT_TRUE(
    // In case (0,1) is executed before (1,2)
    ( hypergraph.contractionTree(0) == 1 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 2 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 1 &&
      hypergraph.referenceCount(2) == 1 ) ||
    // In case (1,2) is executed before (0,1)
    ( hypergraph.contractionTree(0) == 2 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 2 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 2 )
  ) << V(hypergraph.contractionTree(0)) << " "
    << V(hypergraph.contractionTree(1)) << " "
    << V(hypergraph.contractionTree(2)) << " "
    << V(hypergraph.contractionTree(3));
}

TEST_F(ADynamicHypergraph, RegisterContractionsInParallel2) {
  executeParallel([&] {
    ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  }, [&] {
    ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  });

  ASSERT_EQ(2, hypergraph.contractionTree(0));
  ASSERT_EQ(2, hypergraph.contractionTree(1));
  ASSERT_EQ(2, hypergraph.contractionTree(2));
  ASSERT_EQ(0, hypergraph.referenceCount(0));
  ASSERT_EQ(0, hypergraph.referenceCount(1));
  ASSERT_EQ(2, hypergraph.referenceCount(2));
}

TEST_F(ADynamicHypergraph, RegisterContractionsInParallel3) {
  executeParallel([&] {
    ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  }, [&] {
    ASSERT_TRUE(hypergraph.registerContraction(3, 2));
    ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  });

  ASSERT_TRUE(
    // In case (0,2) is executed before (2,3)
    ( hypergraph.contractionTree(0) == 2 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 3 &&
      hypergraph.contractionTree(3) == 3 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 2 &&
      hypergraph.referenceCount(3) == 1 ) ||
    // In case (2,3) is executed before (0,2)
    ( hypergraph.contractionTree(0) == 3 &&
      hypergraph.contractionTree(1) == 3 &&
      hypergraph.contractionTree(2) == 3 &&
      hypergraph.contractionTree(3) == 3 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 0 &&
      hypergraph.referenceCount(3) == 3  )
  ) << V(hypergraph.contractionTree(0)) << " "
    << V(hypergraph.contractionTree(1)) << " "
    << V(hypergraph.contractionTree(2)) << " "
    << V(hypergraph.contractionTree(3));
}

TEST_F(ADynamicHypergraph, RegisterContractionsInParallel4) {
  executeParallel([&] {
    ASSERT_TRUE(hypergraph.registerContraction(2, 0)); // (0)
    ASSERT_TRUE(hypergraph.registerContraction(4, 3)); // (1)
  }, [&] {
    ASSERT_TRUE(hypergraph.registerContraction(3, 2)); // (2)
    ASSERT_TRUE(hypergraph.registerContraction(2, 1)); // (3)
  });

  ASSERT_TRUE(
    // Execution order 0, 1, 2, 3
    ( hypergraph.contractionTree(0) == 2 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 4 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 2 &&
      hypergraph.referenceCount(3) == 0 &&
      hypergraph.referenceCount(4) == 2) ||
    // Execution order 0, 2, 1, 3
    ( hypergraph.contractionTree(0) == 2 &&
      hypergraph.contractionTree(1) == 4 &&
      hypergraph.contractionTree(2) == 3 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 1 &&
      hypergraph.referenceCount(3) == 1 &&
      hypergraph.referenceCount(4) == 2) ||
    // Execution order 0, 2, 3, 1
    ( hypergraph.contractionTree(0) == 2 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 3 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 2 &&
      hypergraph.referenceCount(3) == 1 &&
      hypergraph.referenceCount(4) == 1) ||
    // Execution order 2, 0, 1, 3 or 2, 0, 3, 1 or 2, 3, 0, 1
    ( hypergraph.contractionTree(0) == 3 &&
      hypergraph.contractionTree(1) == 3 &&
      hypergraph.contractionTree(2) == 3 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 0 &&
      hypergraph.referenceCount(2) == 0 &&
      hypergraph.referenceCount(3) == 3 &&
      hypergraph.referenceCount(4) == 1)
  ) << V(hypergraph.contractionTree(0)) << " "
    << V(hypergraph.contractionTree(1)) << " "
    << V(hypergraph.contractionTree(2)) << " "
    << V(hypergraph.contractionTree(3)) << " "
    << V(hypergraph.contractionTree(4));
}

TEST_F(ADynamicHypergraph, RegisterContractionsThatInducesACycleInParallel1) {
  bool succeded_1 = false;
  bool succeded_2 = false;
  executeParallel([&] {
    succeded_1 = hypergraph.registerContraction(0, 1);
  }, [&] {
    succeded_2 = hypergraph.registerContraction(1, 0);
  });

  ASSERT_TRUE((succeded_1 && !succeded_2) || (!succeded_1 && succeded_2));
  if ( succeded_1 ) {
    ASSERT_EQ(0, hypergraph.contractionTree(0));
    ASSERT_EQ(0, hypergraph.contractionTree(1));
    ASSERT_EQ(1, hypergraph.referenceCount(0));
    ASSERT_EQ(0, hypergraph.referenceCount(1));
  } else {
    ASSERT_EQ(1, hypergraph.contractionTree(0));
    ASSERT_EQ(1, hypergraph.contractionTree(1));
    ASSERT_EQ(0, hypergraph.referenceCount(0));
    ASSERT_EQ(1, hypergraph.referenceCount(1));
  }
}

TEST_F(ADynamicHypergraph, RegisterContractionsThatInducesACycleInParallel2) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));

  bool succeded_1 = false;
  bool succeded_2 = false;
  executeParallel([&] {
    succeded_1 = hypergraph.registerContraction(0, 3);
  }, [&] {
    succeded_2 = hypergraph.registerContraction(2, 1);
  });

  ASSERT_TRUE((succeded_1 && !succeded_2) || (!succeded_1 && succeded_2));
  if ( succeded_1 ) {
    ASSERT_EQ(1, hypergraph.contractionTree(0));
    ASSERT_EQ(1, hypergraph.contractionTree(1));
    ASSERT_EQ(3, hypergraph.contractionTree(2));
    ASSERT_EQ(1, hypergraph.contractionTree(3));
    ASSERT_EQ(0, hypergraph.referenceCount(0));
    ASSERT_EQ(2, hypergraph.referenceCount(1));
    ASSERT_EQ(0, hypergraph.referenceCount(2));
    ASSERT_EQ(1, hypergraph.referenceCount(3));
  } else {
    ASSERT_EQ(1, hypergraph.contractionTree(0));
    ASSERT_EQ(3, hypergraph.contractionTree(1));
    ASSERT_EQ(3, hypergraph.contractionTree(2));
    ASSERT_EQ(3, hypergraph.contractionTree(3));
    ASSERT_EQ(0, hypergraph.referenceCount(0));
    ASSERT_EQ(1, hypergraph.referenceCount(1));
    ASSERT_EQ(0, hypergraph.referenceCount(2));
    ASSERT_EQ(2, hypergraph.referenceCount(3));
  }
}

TEST_F(ADynamicHypergraph, RegisterContractionsThatInducesACycleInParallel3) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(4, 3));

  bool succeded_1 = false;
  bool succeded_2 = false;
  bool succeded_3 = false;
  executeParallel([&] {
    succeded_1 = hypergraph.registerContraction(2, 1);
  }, [&] {
    succeded_2 = hypergraph.registerContraction(3, 2);
    succeded_3 = hypergraph.registerContraction(0, 4);
  });

  const size_t num_succeded = succeded_1 + succeded_2 + succeded_3;
  ASSERT_EQ(2, num_succeded);
  if ( succeded_1 ) {
  ASSERT_TRUE(
    // In case (1,2) is executed before (2,3)
    ( hypergraph.contractionTree(0) == 1 &&
      hypergraph.contractionTree(1) == 2 &&
      hypergraph.contractionTree(2) == 4 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 1 &&
      hypergraph.referenceCount(2) == 1 &&
      hypergraph.referenceCount(3) == 0 &&
      hypergraph.referenceCount(4) == 2 ) ||
    // In case (2,3) is executed before (1,2)
    ( hypergraph.contractionTree(0) == 1 &&
      hypergraph.contractionTree(1) == 4 &&
      hypergraph.contractionTree(2) == 4 &&
      hypergraph.contractionTree(3) == 4 &&
      hypergraph.contractionTree(4) == 4 &&
      hypergraph.referenceCount(0) == 0 &&
      hypergraph.referenceCount(1) == 1 &&
      hypergraph.referenceCount(2) == 0 &&
      hypergraph.referenceCount(3) == 0 &&
      hypergraph.referenceCount(4) == 3 )
  ) << V(hypergraph.contractionTree(0)) << " "
    << V(hypergraph.contractionTree(1)) << " "
    << V(hypergraph.contractionTree(2)) << " "
    << V(hypergraph.contractionTree(3)) << " "
    << V(hypergraph.contractionTree(4));
  } else {
    ASSERT_EQ(1, hypergraph.contractionTree(0));
    ASSERT_EQ(1, hypergraph.contractionTree(1));
    ASSERT_EQ(4, hypergraph.contractionTree(2));
    ASSERT_EQ(4, hypergraph.contractionTree(3));
    ASSERT_EQ(1, hypergraph.contractionTree(4));
    ASSERT_EQ(0, hypergraph.referenceCount(0));
    ASSERT_EQ(2, hypergraph.referenceCount(1));
    ASSERT_EQ(0, hypergraph.referenceCount(2));
    ASSERT_EQ(0, hypergraph.referenceCount(3));
    ASSERT_EQ(2, hypergraph.referenceCount(4));
  }
}


using MementoVector = parallel::scalable_vector<Memento>;

void assertEqual(MementoVector actual, MementoVector expected) {
  auto compare = [&](const Memento& lhs, const Memento& rhs) {
    return lhs.u < rhs.u || (lhs.u == rhs.u && lhs.v < rhs.v);
  };
  std::sort(actual.begin(), actual.end(), compare);
  std::sort(expected.begin(), expected.end(), compare);

  ASSERT_EQ(actual.size(), expected.size());
  for ( size_t i = 0; i < actual.size(); ++i ) {
    ASSERT_EQ(actual[i].u, expected[i].u);
    ASSERT_EQ(actual[i].v, expected[i].v);
  }
}

bool assertEqualToOneAlternative(MementoVector actual,
                                 MementoVector alternative_1,
                                 MementoVector alternative_2) {
  auto compare = [&](const Memento& lhs, const Memento& rhs) {
    return lhs.u < rhs.u || (lhs.u == rhs.u && lhs.v < rhs.v);
  };
  std::sort(actual.begin(), actual.end(), compare);
  std::sort(alternative_1.begin(), alternative_1.end(), compare);
  std::sort(alternative_2.begin(), alternative_2.end(), compare);

  bool equal_to_alternative_1 = actual.size() == alternative_1.size();
  bool equal_to_alternative_2 = actual.size() == alternative_2.size();
  if ( equal_to_alternative_1 ) {
    for ( size_t i = 0; i < actual.size(); ++i ) {
      equal_to_alternative_1 = equal_to_alternative_1 &&
        actual[i].u == alternative_1[i].u && actual[i].v == alternative_1[i].v;
    }
  }

  if ( equal_to_alternative_2 ) {
    for ( size_t i = 0; i < actual.size(); ++i ) {
      equal_to_alternative_2 = equal_to_alternative_2 &&
        actual[i].u == alternative_2[i].u && actual[i].v == alternative_2[i].v;
    }
  }

  const size_t num_equal = UI64(equal_to_alternative_1) + UI64(equal_to_alternative_2);
  if ( num_equal == 1 ) {
    return true;
  } else {
    return false;
  }
}

TEST_F(ADynamicHypergraph, PerformsAContraction1) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  assertEqual(hypergraph.contract(0), { Memento { 1, 0 } });

  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_EQ(2, hypergraph.nodeWeight(1));
  ASSERT_EQ(0, hypergraph.referenceCount(1));

  verifyIncidentNets(1, {0, 1});
  verifyPins({ 0, 1, 2, 3 },
    { {1, 2}, {1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformsAContraction2) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  assertEqual(hypergraph.contract(1), { });
  assertEqual(hypergraph.contract(0), { Memento { 1, 0 },  Memento { 2, 1 } });

  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_EQ(3, hypergraph.nodeWeight(2));
  ASSERT_EQ(0, hypergraph.referenceCount(2));

  verifyIncidentNets(2, {0, 1, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {2}, {2, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformsAContraction3) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));

  assertEqual(hypergraph.contract(1), { Memento { 2, 1 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_EQ(2, hypergraph.nodeWeight(2));
  ASSERT_EQ(1, hypergraph.referenceCount(2));

  assertEqual(hypergraph.contract(0), { Memento { 2, 0 },  Memento { 3, 2 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(2));
  ASSERT_EQ(4, hypergraph.nodeWeight(3));
  ASSERT_EQ(0, hypergraph.referenceCount(3));

  verifyIncidentNets(3, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {3}, {3, 4}, {3, 4, 6}, {3, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformsAContraction4) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  ASSERT_TRUE(hypergraph.registerContraction(3, 4));

  assertEqual(hypergraph.contract(1), { Memento { 2, 1 } });
  assertEqual(hypergraph.contract(4), { Memento { 3, 4 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(4));
  ASSERT_EQ(2, hypergraph.nodeWeight(2));
  ASSERT_EQ(2, hypergraph.nodeWeight(3));
  ASSERT_EQ(1, hypergraph.referenceCount(2));
  ASSERT_EQ(1, hypergraph.referenceCount(3));

  assertEqual(hypergraph.contract(0), { Memento { 2, 0 },  Memento { 3, 2 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(2));
  ASSERT_EQ(5, hypergraph.nodeWeight(3));
  ASSERT_EQ(0, hypergraph.referenceCount(3));

  verifyIncidentNets(3, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {3}, {3}, {3, 6}, {3, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformsAContraction5) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  ASSERT_TRUE(hypergraph.registerContraction(3, 4));
  ASSERT_TRUE(hypergraph.registerContraction(6, 3));

  assertEqual(hypergraph.contract(1), { Memento { 2, 1 } });
  assertEqual(hypergraph.contract(4), { Memento { 3, 4 } });
  assertEqual(hypergraph.contract(0), { Memento { 2, 0 },  Memento { 3, 2 }, Memento { 6, 3 } });
  ASSERT_EQ(6, hypergraph.nodeWeight(6));

  verifyIncidentNets(6, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {6}, {6}, {6}, {5, 6} });
}

TEST_F(ADynamicHypergraph, PerformsAContractionWithWeightGreaterThanMaxNodeWeight1) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  ASSERT_EQ(1, hypergraph.contractionTree(0));
  ASSERT_EQ(1, hypergraph.referenceCount(1));
  assertEqual(hypergraph.contract(0, 1), { });
  ASSERT_TRUE(hypergraph.nodeIsEnabled(0));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(1));
  ASSERT_EQ(0, hypergraph.contractionTree(0));
  ASSERT_EQ(0, hypergraph.referenceCount(1));
}

TEST_F(ADynamicHypergraph, PerformsAContractionWithWeightGreaterThanMaxNodeWeight2) {
  ASSERT_TRUE(hypergraph.registerContraction(1, 0));
  assertEqual(hypergraph.contract(0, 2), { Memento { 1, 0 } });

  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  ASSERT_EQ(2, hypergraph.contractionTree(1));
  ASSERT_EQ(3, hypergraph.contractionTree(2));
  ASSERT_EQ(1, hypergraph.referenceCount(2));
  ASSERT_EQ(1, hypergraph.referenceCount(3));
  assertEqual(hypergraph.contract(1, 2), { Memento { 3, 2 } });
  ASSERT_EQ(1, hypergraph.contractionTree(1));
  ASSERT_EQ(0, hypergraph.referenceCount(3));

  verifyIncidentNets(1, {0, 1});
  verifyIncidentNets(3, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {1, 3}, {1, 3, 4}, {3, 4, 6}, {3, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformAContractionsInParallel1) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  executeParallel([&] {
    assertEqual(hypergraph.contract(0), { Memento { 2, 0 } });
  }, [&] {
    assertEqual(hypergraph.contract(1), { Memento { 2, 1 } });
  });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_EQ(3, hypergraph.nodeWeight(2));
  ASSERT_EQ(0, hypergraph.referenceCount(2));

  verifyIncidentNets(2, {0, 1, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {2}, {2, 3, 4}, {3, 4, 6}, {2, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformAContractionsInParallel2) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  MementoVector mementos_1;
  MementoVector mementos_2;
  executeParallel([&] {
    mementos_1 = hypergraph.contract(0);
    ASSERT_TRUE(
      assertEqualToOneAlternative(
        mementos_1,
        { Memento { 2, 0 } },
        { Memento { 2, 0 }, Memento { 3, 2 } }));
  }, [&] {
    mementos_2 = hypergraph.contract(1);
    ASSERT_TRUE(
      assertEqualToOneAlternative(
        mementos_2,
        { Memento { 2, 1 } },
        { Memento { 2, 1 }, Memento { 3, 2 } }));
  });
  MementoVector mementos;
  mementos.insert(mementos.end(), mementos_1.begin(), mementos_1.end());
  mementos.insert(mementos.end(), mementos_2.begin(), mementos_2.end());
  assertEqual(mementos, { Memento { 2, 0 }, Memento { 2, 1 }, Memento { 3, 2 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(2));
  ASSERT_EQ(4, hypergraph.nodeWeight(3));
  ASSERT_EQ(0, hypergraph.referenceCount(3));

  verifyIncidentNets(3, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {3}, {3, 4}, {3, 4, 6}, {3, 5, 6} });
}

TEST_F(ADynamicHypergraph, PerformAContractionsInParallel3) {
  ASSERT_TRUE(hypergraph.registerContraction(2, 0));
  ASSERT_TRUE(hypergraph.registerContraction(2, 1));
  ASSERT_TRUE(hypergraph.registerContraction(3, 2));
  ASSERT_TRUE(hypergraph.registerContraction(3, 4));
  ASSERT_TRUE(hypergraph.registerContraction(6, 3));
  MementoVector mementos_1;
  MementoVector mementos_2;
  std::atomic<size_t> cnt(2);
  executeParallel([&] {
    mementos_1 = hypergraph.contract(0);
    ASSERT_TRUE(
      assertEqualToOneAlternative(
        mementos_1,
        { Memento { 2, 0 } },
        { Memento { 2, 0 }, Memento { 3, 2 } }));
    ++cnt;
  }, [&] {
    mementos_2 = hypergraph.contract(1);
    ASSERT_TRUE(
      assertEqualToOneAlternative(
        mementos_2,
        { Memento { 2, 1 } },
        { Memento { 2, 1 }, Memento { 3, 2 } }));
    ++cnt;
    while ( cnt < 2 ) { }
    assertEqual(hypergraph.contract(4), { Memento { 3, 4 }, Memento { 6, 3 } });
  });
  MementoVector mementos;
  mementos.insert(mementos.end(), mementos_1.begin(), mementos_1.end());
  mementos.insert(mementos.end(), mementos_2.begin(), mementos_2.end());
  assertEqual(mementos, { Memento { 2, 0 }, Memento { 2, 1 }, Memento { 3, 2 } });
  ASSERT_FALSE(hypergraph.nodeIsEnabled(0));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(1));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(2));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(3));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(4));
  ASSERT_EQ(6, hypergraph.nodeWeight(6));
  ASSERT_EQ(0, hypergraph.referenceCount(6));

  verifyIncidentNets(6, {0, 1, 2, 3});
  verifyPins({ 0, 1, 2, 3 },
    { {6}, {6}, {6}, {5, 6} });
}
} // namespace ds
} // namespace mt_kahypar