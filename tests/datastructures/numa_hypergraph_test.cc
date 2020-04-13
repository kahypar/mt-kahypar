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
#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/datastructures/numa_hypergraph_factory.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

// Mocking Numa Architecture (=> 2 NUMA Nodes)
using TypeTraits = TestTypeTraits<2>;
using HwTopology = typename TypeTraits::HwTopology;
using TBB = typename TypeTraits::TBB;

// Define NUMA Hypergraph and Factory
using NumaHyperGraph = NumaHypergraph<StaticHypergraph, HwTopology, TBB>;
using NumaHyperGraphFactory = NumaHypergraphFactory<
  StaticHypergraph, StaticHypergraphFactory, HwTopology, TBB>;

// Define Test Fixture
using AStaticNumaHypergraph = HypergraphFixture<NumaHyperGraph, NumaHyperGraphFactory, TBB>;

TEST_F(AStaticNumaHypergraph, HasCorrectStats) {
  ASSERT_EQ(7,  hypergraph.initialNumNodes());
  ASSERT_EQ(4,  hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(12, hypergraph.initialTotalVertexDegree());
  ASSERT_EQ(5,  hypergraph.initialNumNodes(0));
  ASSERT_EQ(2,  hypergraph.initialNumEdges(0));
  ASSERT_EQ(6,  hypergraph.initialNumPins(0));
  ASSERT_EQ(9,  hypergraph.initialTotalVertexDegree(0));
  ASSERT_EQ(2,  hypergraph.initialNumNodes(1));
  ASSERT_EQ(2,  hypergraph.initialNumEdges(1));
  ASSERT_EQ(6,  hypergraph.initialNumPins(1));
  ASSERT_EQ(3,  hypergraph.initialTotalVertexDegree(1));
  ASSERT_EQ(7,  hypergraph.totalWeight());
  ASSERT_EQ(4,  hypergraph.maxEdgeSize());
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialNodeIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(hypergraph.globalNodeID(expected_hn), hn);
    ++expected_hn;
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(AStaticNumaHypergraph, HasCorrectNodeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHypernode(id[0]);
  hypergraph.disableHypernode(id[5]);
  const std::vector<HypernodeID> expected_iter =
    { id[1], id[2], id[3], id[4], id[6] };
  HypernodeID pos = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(expected_iter[pos++], hn);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialVertexNumaNodeZeroIterator) {
  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : hypergraph.nodes(0) ) {
    ASSERT_EQ(hypergraph.globalNodeID(expected_hn), hn);
    ++expected_hn;
  }
  ASSERT_EQ(5, expected_hn);
}

TEST_F(AStaticNumaHypergraph, HasCorrectEdgeIteratorIfVerticesAreDisabled) {
  hypergraph.disableHyperedge(GLOBAL_EDGE_ID(hypergraph, 0));
  hypergraph.disableHyperedge(GLOBAL_EDGE_ID(hypergraph, 2));
  const std::vector<HyperedgeID> expected_iter =
    { GLOBAL_EDGE_ID(hypergraph, 1), GLOBAL_EDGE_ID(hypergraph, 3) };
  HypernodeID pos = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(expected_iter[pos++], he);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialVertexNumaNodeOneIterator) {
  HypernodeID expected_hn = 5;
  for ( const HypernodeID& hn : hypergraph.nodes(1) ) {
    ASSERT_EQ(hypergraph.globalNodeID(expected_hn), hn);
    ++expected_hn;
  }
  ASSERT_EQ(7, expected_hn);
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialEdgeIterator) {
  HyperedgeID expected_he = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(hypergraph.globalEdgeID(expected_he), he);
    ++expected_he;
  }
  ASSERT_EQ(4, expected_he);
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialEdgeNumaNodeZeroIterator) {
  HyperedgeID expected_he = 0;
  for ( const HyperedgeID& he : hypergraph.edges(0) ) {
    ASSERT_EQ(hypergraph.globalEdgeID(expected_he), he);
    ++expected_he;
  }
  ASSERT_EQ(2, expected_he);
}

TEST_F(AStaticNumaHypergraph, HasCorrectInitialEdgeNumaNodeOneIterator) {
  HyperedgeID expected_he = 2;
  for ( const HyperedgeID& he : hypergraph.edges(1) ) {
    ASSERT_EQ(hypergraph.globalEdgeID(expected_he), he);
    ++expected_he;
  }
  ASSERT_EQ(4, expected_he);
}

TEST_F(AStaticNumaHypergraph, IteratesParallelOverAllNodes) {
  std::vector<uint8_t> visited(7, false);
  hypergraph.doParallelForAllNodes(TBB::GLOBAL_TASK_GROUP,
    [&](const HypernodeID hn) {
      visited[hypergraph.originalNodeID(hn)] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(AStaticNumaHypergraph, IteratesParallelOverAllEdges) {
  std::vector<uint8_t> visited(4, false);
  hypergraph.doParallelForAllEdges(TBB::GLOBAL_TASK_GROUP,
    [&](const HyperedgeID he) {
      visited[hypergraph.originalEdgeID(he)] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(AStaticNumaHypergraph, VerifiesIncidentNets1) {
  verifyIncidentNets(GLOBAL_NODE_ID(hypergraph, 0),
    { GLOBAL_EDGE_ID(hypergraph, 0), GLOBAL_EDGE_ID(hypergraph, 1) });
}

TEST_F(AStaticNumaHypergraph, VerifiesIncidentNets2) {
  verifyIncidentNets(GLOBAL_NODE_ID(hypergraph, 2),
    { GLOBAL_EDGE_ID(hypergraph, 0), GLOBAL_EDGE_ID(hypergraph, 3) });
}

TEST_F(AStaticNumaHypergraph, VerifiesIncidentNets3) {
  verifyIncidentNets(GLOBAL_NODE_ID(hypergraph, 3),
    { GLOBAL_EDGE_ID(hypergraph, 1), GLOBAL_EDGE_ID(hypergraph, 2) });
}

TEST_F(AStaticNumaHypergraph, VerifiesIncidentNets4) {
  verifyIncidentNets(GLOBAL_NODE_ID(hypergraph, 6),
    { GLOBAL_EDGE_ID(hypergraph, 2), GLOBAL_EDGE_ID(hypergraph, 3) });
}

TEST_F(AStaticNumaHypergraph, VerifiesOriginalNodeIDs) {
  for ( HypernodeID id = 0; id < hypergraph.initialNumNodes(); ++id ) {
    ASSERT_EQ(id, hypergraph.originalNodeID(hypergraph.globalNodeID(id)));
  }
}

TEST_F(AStaticNumaHypergraph, VerifiesPinsOfHyperedges) {
  verifyPins({ GLOBAL_EDGE_ID(hypergraph, 0),
               GLOBAL_EDGE_ID(hypergraph, 1),
               GLOBAL_EDGE_ID(hypergraph, 2),
               GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1], id[3], id[4]},
      {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, VerifiesVertexWeights) {
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(1, hypergraph.nodeWeight(hn));
  }
}

TEST_F(AStaticNumaHypergraph, ModifiesNodeWeight) {
  hypergraph.setNodeWeight(id[0], 2);
  hypergraph.setNodeWeight(id[6], 2);
  ASSERT_EQ(2, hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, hypergraph.nodeWeight(id[6]));
  hypergraph.updateTotalWeight(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(9, hypergraph.totalWeight());
}

TEST_F(AStaticNumaHypergraph, VerifiesVertexDegrees) {
  ASSERT_EQ(2, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 0)));
  ASSERT_EQ(1, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 1)));
  ASSERT_EQ(2, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 2)));
  ASSERT_EQ(2, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 3)));
  ASSERT_EQ(2, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 4)));
  ASSERT_EQ(1, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 5)));
  ASSERT_EQ(2, hypergraph.nodeDegree(GLOBAL_NODE_ID(hypergraph, 6)));
}

TEST_F(AStaticNumaHypergraph, RemovesVertices) {
  hypergraph.removeHypernode(id[0]);
  hypergraph.removeHypernode(id[5]);
  ASSERT_EQ(2, hypergraph.numRemovedHypernodes());
}

TEST_F(AStaticNumaHypergraph, VerifiesOriginalEdgeIDs) {
  for ( HyperedgeID id = 0; id < hypergraph.initialNumEdges(); ++id ) {
    ASSERT_EQ(id, hypergraph.originalEdgeID(hypergraph.globalEdgeID(id)));
  }
}

TEST_F(AStaticNumaHypergraph, VerifiesEdgeWeights) {
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(1, hypergraph.edgeWeight(he));
  }
}

TEST_F(AStaticNumaHypergraph, ModifiesEdgeWeight) {
  hypergraph.setEdgeWeight(GLOBAL_EDGE_ID(hypergraph, 0), 2);
  hypergraph.setEdgeWeight(GLOBAL_EDGE_ID(hypergraph, 2), 2);
  ASSERT_EQ(2, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(2, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 2)));
}

TEST_F(AStaticNumaHypergraph, VerifiesEdgeSizes) {
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(4, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 1)));
  ASSERT_EQ(3, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 2)));
  ASSERT_EQ(3, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 3)));
}

TEST_F(AStaticNumaHypergraph, SetsCommunityIDsForEachVertex) {
  hypergraph.setCommunityID(id[0], 1);
  hypergraph.setCommunityID(id[1], 1);
  hypergraph.setCommunityID(id[2], 1);
  hypergraph.setCommunityID(id[3], 2);
  hypergraph.setCommunityID(id[4], 2);
  hypergraph.setCommunityID(id[5], 3);
  hypergraph.setCommunityID(id[6], 3);

  ASSERT_EQ(1, hypergraph.communityID(id[0]));
  ASSERT_EQ(1, hypergraph.communityID(id[1]));
  ASSERT_EQ(1, hypergraph.communityID(id[2]));
  ASSERT_EQ(2, hypergraph.communityID(id[3]));
  ASSERT_EQ(2, hypergraph.communityID(id[4]));
  ASSERT_EQ(3, hypergraph.communityID(id[5]));
  ASSERT_EQ(3, hypergraph.communityID(id[6]));
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectNumberOfCommunities) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunities());
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectNumberOfHypernodesInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(3, hypergraph.numCommunityHypernodes(0));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(1));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(2));
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectNumberOfPinsInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.numCommunityPins(0));
  ASSERT_EQ(4, hypergraph.numCommunityPins(1));
  ASSERT_EQ(3, hypergraph.numCommunityPins(2));
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectCommunityDegreeInEachCommunity) {
  assignCommunityIds();
  ASSERT_EQ(5, hypergraph.communityDegree(0));
  ASSERT_EQ(4, hypergraph.communityDegree(1));
  ASSERT_EQ(3, hypergraph.communityDegree(2));
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityZero) {
  hypergraph.setCommunityID(id[0], 0);
  hypergraph.setCommunityID(id[1], 0);
  hypergraph.setCommunityID(id[2], 0);
  hypergraph.setCommunityID(id[3], 1);
  hypergraph.setCommunityID(id[4], 2);
  hypergraph.setCommunityID(id[5], 1);
  hypergraph.setCommunityID(id[6], 2);
  hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(id[0]), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[0])]);
  flag[hypergraph.communityNodeId(id[0])] = true;
  ASSERT_LE(hypergraph.communityNodeId(id[1]), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[1])]);
  flag[hypergraph.communityNodeId(id[1])] = true;
  ASSERT_LE(hypergraph.communityNodeId(id[2]), 2);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[2])]);
  flag[hypergraph.communityNodeId(id[2])] = true;
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityOne) {
  hypergraph.setCommunityID(id[0], 0);
  hypergraph.setCommunityID(id[1], 0);
  hypergraph.setCommunityID(id[2], 0);
  hypergraph.setCommunityID(id[3], 1);
  hypergraph.setCommunityID(id[4], 2);
  hypergraph.setCommunityID(id[5], 1);
  hypergraph.setCommunityID(id[6], 2);
  hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(id[3]), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[3])]);
  flag[hypergraph.communityNodeId(id[3])] = true;
  ASSERT_LE(hypergraph.communityNodeId(id[5]), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[5])]);
  flag[hypergraph.communityNodeId(id[5])] = true;
}

TEST_F(AStaticNumaHypergraph, ComputesCorrectCommunityNodeIdsOfCommunityTwo) {
  hypergraph.setCommunityID(id[0], 0);
  hypergraph.setCommunityID(id[1], 0);
  hypergraph.setCommunityID(id[2], 0);
  hypergraph.setCommunityID(id[3], 1);
  hypergraph.setCommunityID(id[4], 2);
  hypergraph.setCommunityID(id[5], 1);
  hypergraph.setCommunityID(id[6], 2);
  hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  std::vector<bool> flag(3, false);
  ASSERT_LE(hypergraph.communityNodeId(id[4]), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[4])]);
  flag[hypergraph.communityNodeId(id[4])] = true;
  ASSERT_LE(hypergraph.communityNodeId(id[6]), 1);
  ASSERT_FALSE(flag[hypergraph.communityNodeId(id[6])]);
  flag[hypergraph.communityNodeId(id[6])] = true;
}

TEST_F(AStaticNumaHypergraph, VerifiesNumberOfCommunitiesInHyperedges) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)));
}

TEST_F(AStaticNumaHypergraph, ChecksHyperedgeCommunityIterator1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(GLOBAL_EDGE_ID(hypergraph, 0)) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, ChecksHyperedgeCommunityIterator2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 1 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(GLOBAL_EDGE_ID(hypergraph, 1)) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, ChecksHyperedgeCommunityIterator3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 1, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(GLOBAL_EDGE_ID(hypergraph, 2)) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, ChecksHyperedgeCommunityIterator4) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  const std::vector<PartitionID> expected_iter = { 0, 2 };
  size_t pos = 0;
  for ( const PartitionID& community_id : hypergraph.communities(GLOBAL_EDGE_ID(hypergraph, 3)) ) {
    ASSERT_EQ(expected_iter[pos++], community_id);
  }
  ASSERT_EQ(expected_iter.size(), pos);
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 0), 0));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 0), 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 1), 0));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 1), 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 3), 0));
  ASSERT_EQ(1, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 3), 0));
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 1), 1));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 1), 1));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 2), 1));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 2), 1));
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 2), 2));
  ASSERT_EQ(1, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 2), 2));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 3), 2));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 3), 2));
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { GLOBAL_EDGE_ID(hypergraph, 0),
                           GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1]}, {id[2]} });
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  verifyCommunityPins(1, { GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 2) },
    { {id[3], id[4]}, {id[3], id[4]} });
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  verifyCommunityPins(2, { GLOBAL_EDGE_ID(hypergraph, 2),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[6]}, {id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, RemovesAHyperedgeFromTheHypergraph1) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 0));
  verifyIncidentNets(id[0], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[2], { GLOBAL_EDGE_ID(hypergraph, 3) });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(GLOBAL_EDGE_ID(hypergraph, 0), he);
  }
}

TEST_F(AStaticNumaHypergraph, RemovesAHyperedgeFromTheHypergraph2) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 1));
  verifyIncidentNets(id[0], { GLOBAL_EDGE_ID(hypergraph, 0) });
  verifyIncidentNets(id[1], { });
  verifyIncidentNets(id[3], { GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(id[4], { GLOBAL_EDGE_ID(hypergraph, 2) });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(GLOBAL_EDGE_ID(hypergraph, 1), he);
  }
}

TEST_F(AStaticNumaHypergraph, RemovesAHyperedgeFromTheHypergraph3) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 2));
  verifyIncidentNets(id[3], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[4], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[6], { GLOBAL_EDGE_ID(hypergraph, 3) });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(GLOBAL_EDGE_ID(hypergraph, 2), he);
  }
}

TEST_F(AStaticNumaHypergraph, RemovesAHyperedgeFromTheHypergraph4) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 3));
  verifyIncidentNets(id[2], { GLOBAL_EDGE_ID(hypergraph, 0) });
  verifyIncidentNets(id[5], { });
  verifyIncidentNets(id[6], { GLOBAL_EDGE_ID(hypergraph, 2) });
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_NE(GLOBAL_EDGE_ID(hypergraph, 3), he);
  }
}

TEST_F(AStaticNumaHypergraph, RestoresARemovedHyperedge1) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 0));
  hypergraph.restoreEdge(GLOBAL_EDGE_ID(hypergraph, 0), 2);
  verifyIncidentNets(id[0], { GLOBAL_EDGE_ID(hypergraph, 0),
                              GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[2], { GLOBAL_EDGE_ID(hypergraph, 0),
                              GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyPins({ GLOBAL_EDGE_ID(hypergraph, 0) },
    { {id[0], id[2]} });
}

TEST_F(AStaticNumaHypergraph, RestoresARemovedHyperedge2) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 1));
  hypergraph.restoreEdge(GLOBAL_EDGE_ID(hypergraph, 1), 4);
  verifyIncidentNets(id[0], { GLOBAL_EDGE_ID(hypergraph, 0),
                              GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[1], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(id[3], { GLOBAL_EDGE_ID(hypergraph, 1),
                              GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(id[4], { GLOBAL_EDGE_ID(hypergraph, 1),
                              GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyPins({ GLOBAL_EDGE_ID(hypergraph, 1) },
    { {id[0], id[1], id[3], id[4]} });
}

TEST_F(AStaticNumaHypergraph, RestoresARemovedHyperedge3) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 2));
  hypergraph.restoreEdge(GLOBAL_EDGE_ID(hypergraph, 2), 3);
  verifyIncidentNets(id[3], { GLOBAL_EDGE_ID(hypergraph, 1),
                              GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(id[4], { GLOBAL_EDGE_ID(hypergraph, 1),
                              GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(id[6], { GLOBAL_EDGE_ID(hypergraph, 2),
                              GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyPins({ GLOBAL_EDGE_ID(hypergraph, 2) },
    { {id[3], id[4], id[6]} });
}

TEST_F(AStaticNumaHypergraph, RestoresARemovedHyperedge4) {
  hypergraph.removeEdge(GLOBAL_EDGE_ID(hypergraph, 3));
  hypergraph.restoreEdge(GLOBAL_EDGE_ID(hypergraph, 3), 3);
  verifyIncidentNets(id[2], { GLOBAL_EDGE_ID(hypergraph, 0),
                              GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(id[5], { GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(id[6], { GLOBAL_EDGE_ID(hypergraph, 2),
                              GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyPins({ GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[2], id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, ComparesStatsIfCopiedParallel) {
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(AStaticNumaHypergraph, ComparesStatsIfCopiedSequential) {
  NumaHyperGraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.initialNumNodes(), copy_hg.initialNumNodes());
  ASSERT_EQ(hypergraph.initialNumEdges(), copy_hg.initialNumEdges());
  ASSERT_EQ(hypergraph.initialNumPins(), copy_hg.initialNumPins());
  ASSERT_EQ(hypergraph.initialTotalVertexDegree(), copy_hg.initialTotalVertexDegree());
  ASSERT_EQ(hypergraph.totalWeight(), copy_hg.totalWeight());
  ASSERT_EQ(hypergraph.maxEdgeSize(), copy_hg.maxEdgeSize());
}

TEST_F(AStaticNumaHypergraph, ComparesIncidentNetsIfCopiedParallel) {
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  verifyIncidentNets(copy_hg, id[0], { GLOBAL_EDGE_ID(hypergraph, 0),
                                       GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(copy_hg, id[1], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(copy_hg, id[2], { GLOBAL_EDGE_ID(hypergraph, 0),
                                       GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(copy_hg, id[3], { GLOBAL_EDGE_ID(hypergraph, 1),
                                       GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(copy_hg, id[4], { GLOBAL_EDGE_ID(hypergraph, 1),
                                       GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(copy_hg, id[5], { GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(copy_hg, id[6], { GLOBAL_EDGE_ID(hypergraph, 2),
                                       GLOBAL_EDGE_ID(hypergraph, 3) });
}

TEST_F(AStaticNumaHypergraph, ComparesIncidentNetsIfCopiedSequential) {
  NumaHyperGraph copy_hg = hypergraph.copy();
  verifyIncidentNets(copy_hg, id[0], { GLOBAL_EDGE_ID(hypergraph, 0),
                                       GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(copy_hg, id[1], { GLOBAL_EDGE_ID(hypergraph, 1) });
  verifyIncidentNets(copy_hg, id[2], { GLOBAL_EDGE_ID(hypergraph, 0),
                                       GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(copy_hg, id[3], { GLOBAL_EDGE_ID(hypergraph, 1),
                                       GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(copy_hg, id[4], { GLOBAL_EDGE_ID(hypergraph, 1),
                                       GLOBAL_EDGE_ID(hypergraph, 2) });
  verifyIncidentNets(copy_hg, id[5], { GLOBAL_EDGE_ID(hypergraph, 3) });
  verifyIncidentNets(copy_hg, id[6], { GLOBAL_EDGE_ID(hypergraph, 2),
                                       GLOBAL_EDGE_ID(hypergraph, 3) });
}

TEST_F(AStaticNumaHypergraph, ComparesPinsOfHyperedgesIfCopiedParallel) {
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  verifyPins(copy_hg, { GLOBAL_EDGE_ID(hypergraph, 0),
                        GLOBAL_EDGE_ID(hypergraph, 1),
                        GLOBAL_EDGE_ID(hypergraph, 2),
                        GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1], id[3], id[4]},
      {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, ComparesPinsOfHyperedgesIfCopiedSequential) {
  NumaHyperGraph copy_hg = hypergraph.copy();
  verifyPins(copy_hg, { GLOBAL_EDGE_ID(hypergraph, 0),
                        GLOBAL_EDGE_ID(hypergraph, 1),
                        GLOBAL_EDGE_ID(hypergraph, 2),
                        GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1], id[3], id[4]},
      {id[3], id[4], id[6]}, {id[2], id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, ComparesCommunityStatsIfCopiedParallel) {
  assignCommunityIds();
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
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

TEST_F(AStaticNumaHypergraph, ComparesCommunityStatsIfCopiedSequential) {
  assignCommunityIds();
  NumaHyperGraph copy_hg = hypergraph.copy();
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

TEST_F(AStaticNumaHypergraph, ComparesCommunityIdsIfCopiedParallel) {
  assignCommunityIds();
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.communityID(id[0]), copy_hg.communityID(id[0]));
  ASSERT_EQ(hypergraph.communityID(id[1]), copy_hg.communityID(id[1]));
  ASSERT_EQ(hypergraph.communityID(id[2]), copy_hg.communityID(id[2]));
  ASSERT_EQ(hypergraph.communityID(id[3]), copy_hg.communityID(id[3]));
  ASSERT_EQ(hypergraph.communityID(id[4]), copy_hg.communityID(id[4]));
  ASSERT_EQ(hypergraph.communityID(id[5]), copy_hg.communityID(id[5]));
  ASSERT_EQ(hypergraph.communityID(id[6]), copy_hg.communityID(id[6]));
}

TEST_F(AStaticNumaHypergraph, ComparesCommunityIdsIfCopiedSequential) {
  assignCommunityIds();
  NumaHyperGraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.communityID(id[0]), copy_hg.communityID(id[0]));
  ASSERT_EQ(hypergraph.communityID(id[1]), copy_hg.communityID(id[1]));
  ASSERT_EQ(hypergraph.communityID(id[2]), copy_hg.communityID(id[2]));
  ASSERT_EQ(hypergraph.communityID(id[3]), copy_hg.communityID(id[3]));
  ASSERT_EQ(hypergraph.communityID(id[4]), copy_hg.communityID(id[4]));
  ASSERT_EQ(hypergraph.communityID(id[5]), copy_hg.communityID(id[5]));
  ASSERT_EQ(hypergraph.communityID(id[6]), copy_hg.communityID(id[6]));
}

TEST_F(AStaticNumaHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)));
}

TEST_F(AStaticNumaHypergraph, ComparesNumberOfCommunitiesInHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  NumaHyperGraph copy_hg = hypergraph.copy();
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)));
  ASSERT_EQ(hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)),
            copy_hg.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)));
}

TEST_F(AStaticNumaHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedParallel) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  NumaHyperGraph copy_hg = hypergraph.copy(TBB::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { GLOBAL_EDGE_ID(hypergraph, 0),
                           GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1]}, {id[2]} });
  verifyCommunityPins(1, { GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 2) },
    { {id[3], id[4]}, {id[3], id[4]} });
  verifyCommunityPins(2, { GLOBAL_EDGE_ID(hypergraph, 2),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[6]}, {id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, ComparesPinsOfCommunityHyperedgesIfCopiedSequential) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBB::GLOBAL_TASK_GROUP);
  NumaHyperGraph copy_hg = hypergraph.copy();
  verifyCommunityPins(0, { GLOBAL_EDGE_ID(hypergraph, 0),
                           GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1]}, {id[2]} });
  verifyCommunityPins(1, { GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 2) },
    { {id[3], id[4]}, {id[3], id[4]} });
  verifyCommunityPins(2, { GLOBAL_EDGE_ID(hypergraph, 2),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[6]}, {id[5], id[6]} });
}

TEST_F(AStaticNumaHypergraph, ContractsCommunities1) {
  parallel::scalable_vector<HypernodeID> c_mapping = {1, 4, 1, 5, 5, 4, 5};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0),
    GLOBAL_ID(c_hypergraph, 1), GLOBAL_ID(c_hypergraph, 2) };

  // Verify Mapping
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[0]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[1]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[2]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[3]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[4]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[5]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[6]));

  // Verify Stats
  ASSERT_EQ(3, c_hypergraph.initialNumNodes());
  ASSERT_EQ(1, c_hypergraph.initialNumEdges());
  ASSERT_EQ(3, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(3, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[1]));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(id[2]));

  // Verify Hyperedge Weights
  ASSERT_EQ(2, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(0)));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, id[0], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyIncidentNets(c_hypergraph, id[1], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyIncidentNets(c_hypergraph, id[2], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyPins(c_hypergraph, { GLOBAL_EDGE_ID(c_hypergraph, 0) },
    { {id[0], id[1], id[2]} });
}

TEST_F(AStaticNumaHypergraph, ContractsCommunities2) {
  parallel::scalable_vector<HypernodeID> c_mapping = {1, 4, 1, 5, 5, 6, 5};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0),
    GLOBAL_ID(c_hypergraph, 1), GLOBAL_ID(c_hypergraph, 2), GLOBAL_ID(c_hypergraph, 3) };

  // Verify Mapping
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[0]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[1]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[2]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[3]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[4]));
  ASSERT_EQ(3, c_hypergraph.originalNodeID(c_mapping[5]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[6]));

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(2, c_hypergraph.initialNumEdges());
  ASSERT_EQ(6, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(3, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(id[1]));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(id[2]));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(id[3]));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(0)));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(1)));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, id[0],
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyIncidentNets(c_hypergraph, id[1],
    { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyIncidentNets(c_hypergraph, id[2],
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyIncidentNets(c_hypergraph, id[3],
    { GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyPins(c_hypergraph,
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) },
    { {id[0], id[1], id[2]}, {id[0], id[2], id[3]} });
}

TEST_F(AStaticNumaHypergraph, ContractsCommunities3) {
  parallel::scalable_vector<HypernodeID> c_mapping = {2, 2, 0, 5, 5, 1, 1};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0),
    GLOBAL_ID(c_hypergraph, 1), GLOBAL_ID(c_hypergraph, 2), GLOBAL_ID(c_hypergraph, 3) };

  // Verify Mapping
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[0]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[1]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[2]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[3]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[4]));
  ASSERT_EQ(3, c_hypergraph.originalNodeID(c_mapping[5]));
  ASSERT_EQ(3, c_hypergraph.originalNodeID(c_mapping[6]));

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(4, c_hypergraph.initialNumEdges());
  ASSERT_EQ(8, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(id[1]));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[2]));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[3]));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(0)));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(1)));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(2)));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(3)));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, id[0],
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyIncidentNets(c_hypergraph, id[1],
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 3) });
  verifyIncidentNets(c_hypergraph, id[2],
    { GLOBAL_EDGE_ID(c_hypergraph, 1), GLOBAL_EDGE_ID(c_hypergraph, 2) });
  verifyIncidentNets(c_hypergraph, id[3],
    { GLOBAL_EDGE_ID(c_hypergraph, 2), GLOBAL_EDGE_ID(c_hypergraph, 3) });
  verifyPins(c_hypergraph,
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1),
      GLOBAL_EDGE_ID(c_hypergraph, 2), GLOBAL_EDGE_ID(c_hypergraph, 3) },
    { { id[0], id[1] }, { id[0], id[2] }, { id[2], id[3] }, { id[1], id[3] } });
}

TEST_F(AStaticNumaHypergraph, ContractsCommunitiesWithDisabledHypernodes) {
  hypergraph.disableHypernode(id[0]);
  hypergraph.disableHypernode(id[6]);

  parallel::scalable_vector<HypernodeID> c_mapping = {0, 1, 1, 2, 2, 2, 6};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0), GLOBAL_ID(c_hypergraph, 1) };

  // Verify Mapping
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[1]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[2]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[3]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[4]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[5]));

  // Verify Stats
  ASSERT_EQ(2, c_hypergraph.initialNumNodes());
  ASSERT_EQ(1, c_hypergraph.initialNumEdges());
  ASSERT_EQ(2, c_hypergraph.initialNumPins());
  ASSERT_EQ(5, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(3, c_hypergraph.nodeWeight(id[1]));

  // Verify Hyperedge Weights
  ASSERT_EQ(2, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(0)));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, id[0], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyIncidentNets(c_hypergraph, id[1], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyPins(c_hypergraph, { GLOBAL_EDGE_ID(c_hypergraph, 0) },
    { {id[0], id[1]} });
}

TEST_F(AStaticNumaHypergraph, ContractsCommunitiesWithDisabledHyperedges) {
  hypergraph.disableHyperedge(GLOBAL_EDGE_ID(hypergraph, 3));

  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 0, 1, 1, 2, 3};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0), GLOBAL_ID(c_hypergraph, 1),
    GLOBAL_ID(c_hypergraph, 2), GLOBAL_ID(c_hypergraph, 3) };

  // Verify Mapping
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[0]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[1]));
  ASSERT_EQ(0, c_hypergraph.originalNodeID(c_mapping[2]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[3]));
  ASSERT_EQ(1, c_hypergraph.originalNodeID(c_mapping[4]));
  ASSERT_EQ(2, c_hypergraph.originalNodeID(c_mapping[5]));
  ASSERT_EQ(3, c_hypergraph.originalNodeID(c_mapping[6]));

  // Verify Stats
  ASSERT_EQ(4, c_hypergraph.initialNumNodes());
  ASSERT_EQ(2, c_hypergraph.initialNumEdges());
  ASSERT_EQ(4, c_hypergraph.initialNumPins());
  ASSERT_EQ(7, c_hypergraph.totalWeight());
  ASSERT_EQ(2, c_hypergraph.maxEdgeSize());

  // Verify Vertex Weights
  ASSERT_EQ(3, c_hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, c_hypergraph.nodeWeight(id[1]));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(id[2]));
  ASSERT_EQ(1, c_hypergraph.nodeWeight(id[3]));

  // Verify Hyperedge Weights
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(0)));
  ASSERT_EQ(1, c_hypergraph.edgeWeight(c_hypergraph.globalEdgeID(1)));

  // Verify Hypergraph Structure
  verifyIncidentNets(c_hypergraph, id[0], { GLOBAL_EDGE_ID(c_hypergraph, 0) });
  verifyIncidentNets(c_hypergraph, id[1],
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyIncidentNets(c_hypergraph, id[2], { });
  verifyIncidentNets(c_hypergraph, id[3], { GLOBAL_EDGE_ID(c_hypergraph, 1) });
  verifyPins(c_hypergraph,
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1) },
    { {id[0], id[1]}, {id[1], id[3]} });
}

TEST_F(AStaticNumaHypergraph, ContractCommunitiesIfCommunityInformationAreAvailable) {
  assignCommunityIds();
  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 1, 2, 2, 3, 3};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0),
    GLOBAL_ID(c_hypergraph, 1), GLOBAL_ID(c_hypergraph, 2), GLOBAL_ID(c_hypergraph, 3) };

  // Verify Community Ids
  ASSERT_EQ(0, c_hypergraph.communityID(id[0]));
  ASSERT_EQ(0, c_hypergraph.communityID(id[1]));
  ASSERT_EQ(1, c_hypergraph.communityID(id[2]));
  ASSERT_EQ(2, c_hypergraph.communityID(id[3]));

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

TEST_F(AStaticNumaHypergraph, ContractCommunitiesIfCommunityHyperedgesAreAvailable) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  parallel::scalable_vector<HypernodeID> c_mapping = {0, 0, 1, 2, 2, 3, 3};
  NumaHyperGraph c_hypergraph = hypergraph.contract(
    c_mapping, TBBNumaArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(c_hypergraph, 0),
    GLOBAL_ID(c_hypergraph, 1), GLOBAL_ID(c_hypergraph, 2), GLOBAL_ID(c_hypergraph, 3) };

  // Verify Hypergraph Structure
  verifyCommunityPins(c_hypergraph, 0,
    { GLOBAL_EDGE_ID(c_hypergraph, 0), GLOBAL_EDGE_ID(c_hypergraph, 1),
      GLOBAL_EDGE_ID(c_hypergraph, 3) },
    { {id[0], id[1]}, {id[0]}, {id[1]} });
  verifyCommunityPins(c_hypergraph, 1,
    { GLOBAL_EDGE_ID(c_hypergraph, 1), GLOBAL_EDGE_ID(c_hypergraph, 2) },
    { {id[2]}, {id[2]} });
  verifyCommunityPins(c_hypergraph, 2,
    { GLOBAL_EDGE_ID(c_hypergraph, 2), GLOBAL_EDGE_ID(c_hypergraph, 3) },
    { {id[3]}, {id[3]} });
}

}
} // namespace mt_kahypar