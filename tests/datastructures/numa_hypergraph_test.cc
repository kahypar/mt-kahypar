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
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/datastructures/numa_hypergraph_factory.h"

#include "tests/datastructures/hypergraph_fixtures.h"

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
  hypergraph.doParallelForAllNodes(TBBNumaArena::GLOBAL_TASK_GROUP,
    [&](const HypernodeID hn) {
      visited[hypergraph.originalNodeID(hn)] = true;
    });

  for ( size_t i = 0; i < visited.size(); ++i ) {
    ASSERT_TRUE(visited[i]) << i;
  }
}

TEST_F(AStaticNumaHypergraph, IteratesParallelOverAllEdges) {
  std::vector<uint8_t> visited(4, false);
  hypergraph.doParallelForAllEdges(TBBNumaArena::GLOBAL_TASK_GROUP,
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
  hypergraph.updateTotalWeight(TBBNumaArena::GLOBAL_TASK_GROUP);
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

TEST_F(AStaticNumaHypergraph, MarksVerticesWithADegreeGreaterOrEqualThanTwoAsHighDegree) {
  hypergraph.markAllHighDegreeVertices(TBBNumaArena::GLOBAL_TASK_GROUP, 2UL);
  ASSERT_TRUE(hypergraph.isHighDegreeVertex(id[0]));
  ASSERT_FALSE(hypergraph.isHighDegreeVertex(id[1]));
  ASSERT_TRUE(hypergraph.isHighDegreeVertex(id[2]));
  ASSERT_TRUE(hypergraph.isHighDegreeVertex(id[3]));
  ASSERT_TRUE(hypergraph.isHighDegreeVertex(id[4]));
  ASSERT_FALSE(hypergraph.isHighDegreeVertex(id[5]));
  ASSERT_TRUE(hypergraph.isHighDegreeVertex(id[6]));
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

TEST_F(AStaticNumaHypergraph, VerifiesNumberOfCommunitiesInHyperedges) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 0)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 1)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 2)));
  ASSERT_EQ(2, hypergraph.numCommunitiesInHyperedge(GLOBAL_EDGE_ID(hypergraph, 3)));
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 0), 0));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 0), 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 1), 0));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 1), 0));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 3), 0));
  ASSERT_EQ(1, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 3), 0));
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 1), 1));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 1), 1));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 2), 1));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 2), 1));
}

TEST_F(AStaticNumaHypergraph, VerifiesCommunityHyperedgeInternals3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 2), 2));
  ASSERT_EQ(1, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 2), 2));
  ASSERT_EQ(1, hypergraph.edgeWeight(GLOBAL_EDGE_ID(hypergraph, 3), 2));
  ASSERT_EQ(2, hypergraph.edgeSize(GLOBAL_EDGE_ID(hypergraph, 3), 2));
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge1) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(0, { GLOBAL_EDGE_ID(hypergraph, 0),
                           GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[0], id[2]}, {id[0], id[1]}, {id[2]} });
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge2) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(1, { GLOBAL_EDGE_ID(hypergraph, 1),
                           GLOBAL_EDGE_ID(hypergraph, 2) },
    { {id[3], id[4]}, {id[3], id[4]} });
}

TEST_F(AStaticNumaHypergraph, CanIterateOverThePinsOfASpecificCommunityInAHyperedge3) {
  assignCommunityIds();
  hypergraph.initializeCommunityHyperedges(TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyCommunityPins(2, { GLOBAL_EDGE_ID(hypergraph, 2),
                           GLOBAL_EDGE_ID(hypergraph, 3) },
    { {id[6]}, {id[5], id[6]} });
}

}
} // namespace mt_kahypar