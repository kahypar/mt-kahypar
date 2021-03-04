/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#include "mt-kahypar/partition/refinement/ilp/ilp_hypergraph.h"

using ::testing::Test;

namespace mt_kahypar {

class AILPHypergraph : public Test {
public:

  AILPHypergraph() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
        7, 4, {{0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6}})),
    phg() {
    phg = PartitionedHypergraph(2, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
    phg.setOnlyNodePart(0, 0);
    phg.setOnlyNodePart(1, 0);
    phg.setOnlyNodePart(2, 0);
    phg.setOnlyNodePart(3, 0);
    phg.setOnlyNodePart(4, 1);
    phg.setOnlyNodePart(5, 1);
    phg.setOnlyNodePart(6, 1);
    phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  Hypergraph hg;
  PartitionedHypergraph phg;
};

TEST_F(AILPHypergraph, ChecksNumberOfNodesAndEdges) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);
  ASSERT_EQ(3, ilp_hg.numNodes());
  ASSERT_EQ(2, ilp_hg.numEdges());
}

TEST_F(AILPHypergraph, IteratesOverNodes) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  HypernodeID expected_hn = 0;
  for ( const HypernodeID& hn : ilp_hg.nodes() ) {
    ASSERT_EQ(expected_hn++, hn);
  }
  ASSERT_EQ(3, expected_hn);
}

TEST_F(AILPHypergraph, IteratesOverEdges) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  HyperedgeID expected_he = 0;
  for ( const HyperedgeID& he : ilp_hg.edges() ) {
    ASSERT_EQ(expected_he++, he);
  }
  ASSERT_EQ(2, expected_he);
}

TEST_F(AILPHypergraph, IteratesOverPins1) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  const std::vector<HypernodeID> expected_pins = {3, 0, 1, 2};
  for ( const HypernodeID& pin : ilp_hg.pins(0) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, IteratesOverPins2) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  const std::vector<HypernodeID> expected_pins = {1, 2, 4};
  for ( const HypernodeID& pin : ilp_hg.pins(1) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, IteratesOverPins3) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  const std::vector<HypernodeID> expected_pins = {2, 0, 1, 3};
  for ( const HypernodeID& pin : ilp_hg.pins(0) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, IteratesOverPins4) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  const std::vector<HypernodeID> expected_pins = {1, 3};
  for ( const HypernodeID& pin : ilp_hg.pins(1) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, IteratesOverPinsTwoTimes1) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  std::vector<HypernodeID> expected_pins = {2, 0, 1, 3};
  for ( const HypernodeID& pin : ilp_hg.pins(0) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);

  idx = 0;
  expected_pins = {1, 3};
  for ( const HypernodeID& pin : ilp_hg.pins(1) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, IteratesOverPinsTwoTimes2) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  size_t idx = 0;
  std::vector<HypernodeID> expected_pins = {2, 0, 1, 3};
  for ( const HypernodeID& pin : ilp_hg.pins(0) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);

  idx = 0;
  for ( const HypernodeID& pin : ilp_hg.pins(0) ) {
    ASSERT_EQ(expected_pins[idx++], pin);
  }
  ASSERT_EQ(expected_pins.size(), idx);
}

TEST_F(AILPHypergraph, VerifyNodeWeights1) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(1, ilp_hg.nodeWeight(0));
  ASSERT_EQ(1, ilp_hg.nodeWeight(1));
  ASSERT_EQ(1, ilp_hg.nodeWeight(2));
}

TEST_F(AILPHypergraph, VerifySuperVertexWeights1) {
  vec<HypernodeID> nodes = { 1, 3, 4 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(2, ilp_hg.superVertexWeight(0));
  ASSERT_EQ(2, ilp_hg.superVertexWeight(1));
}


TEST_F(AILPHypergraph, VerifyNodeWeights2) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(1, ilp_hg.nodeWeight(0));
  ASSERT_EQ(1, ilp_hg.nodeWeight(1));
}

TEST_F(AILPHypergraph, VerifySuperVertexWeights2) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(2, ilp_hg.superVertexWeight(0));
  ASSERT_EQ(3, ilp_hg.superVertexWeight(1));
}

TEST_F(AILPHypergraph, VerifyEdgeWeights1) {
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(1, ilp_hg.edgeWeight(0));
  ASSERT_EQ(1, ilp_hg.edgeWeight(1));
}

TEST_F(AILPHypergraph, VerifyEdgeWeights2) {
  vec<HypernodeID> nodes = { 0, 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_EQ(1, ilp_hg.edgeWeight(0));
  ASSERT_EQ(1, ilp_hg.edgeWeight(1));
  ASSERT_EQ(1, ilp_hg.edgeWeight(2));
}

TEST_F(AILPHypergraph, CheckIfEdgeContainsAPinInBlock) {
  vec<HypernodeID> nodes = { 0, 1, 3 };
  ILPHypergraph ilp_hg(phg, nodes);

  ASSERT_TRUE(ilp_hg.containsPinInPart(0, 0));
  ASSERT_FALSE(ilp_hg.containsPinInPart(0, 1));
  ASSERT_TRUE(ilp_hg.containsPinInPart(1, 0));
  ASSERT_TRUE(ilp_hg.containsPinInPart(1, 1));
  ASSERT_TRUE(ilp_hg.containsPinInPart(2, 0));
  ASSERT_TRUE(ilp_hg.containsPinInPart(2, 1));
}

TEST_F(AILPHypergraph, MapsBlockIdsCorrectly1) {
  PartitionedHypergraph phg2(3, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  phg2.setOnlyNodePart(0, 0);
  phg2.setOnlyNodePart(1, 0);
  phg2.setOnlyNodePart(2, 2);
  phg2.setOnlyNodePart(3, 1);
  phg2.setOnlyNodePart(4, 1);
  phg2.setOnlyNodePart(5, 2);
  phg2.setOnlyNodePart(6, 1);
  phg2.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  vec<HypernodeID> nodes = { 1, 3 };
  ILPHypergraph ilp_hg(phg2, nodes);

  ASSERT_EQ(2, ilp_hg.k());
  ASSERT_EQ(0, ilp_hg.partID(0));
  ASSERT_EQ(1, ilp_hg.partID(1));
  ASSERT_EQ(0, ilp_hg.superVertexBlock(2));
  ASSERT_EQ(1, ilp_hg.superVertexBlock(3));
}

TEST_F(AILPHypergraph, MapsBlockIdsCorrectly2) {
  PartitionedHypergraph phg2(3, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  phg2.setOnlyNodePart(0, 0);
  phg2.setOnlyNodePart(1, 0);
  phg2.setOnlyNodePart(2, 2);
  phg2.setOnlyNodePart(3, 1);
  phg2.setOnlyNodePart(4, 1);
  phg2.setOnlyNodePart(5, 2);
  phg2.setOnlyNodePart(6, 2);
  phg2.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  vec<HypernodeID> nodes = { 2, 5 };
  ILPHypergraph ilp_hg(phg2, nodes);

  ASSERT_EQ(2, ilp_hg.k());
  ASSERT_EQ(1, ilp_hg.partID(0));
  ASSERT_EQ(1, ilp_hg.partID(1));
  ASSERT_EQ(0, ilp_hg.superVertexBlock(2));
  ASSERT_EQ(1, ilp_hg.superVertexBlock(3));
}

TEST_F(AILPHypergraph, MapsBlockIdsCorrectly3) {
  PartitionedHypergraph phg2(3, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
  phg2.setOnlyNodePart(0, 0);
  phg2.setOnlyNodePart(1, 0);
  phg2.setOnlyNodePart(2, 2);
  phg2.setOnlyNodePart(3, 1);
  phg2.setOnlyNodePart(4, 1);
  phg2.setOnlyNodePart(5, 2);
  phg2.setOnlyNodePart(6, 1);
  phg2.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  vec<HypernodeID> nodes = { 0, 1, 3 };
  ILPHypergraph ilp_hg(phg2, nodes);

  ASSERT_EQ(3, ilp_hg.k());
  ASSERT_EQ(0, ilp_hg.partID(0));
  ASSERT_EQ(0, ilp_hg.partID(1));
  ASSERT_EQ(1, ilp_hg.partID(2));
  ASSERT_EQ(0, ilp_hg.superVertexBlock(3));
  ASSERT_EQ(1, ilp_hg.superVertexBlock(4));
  ASSERT_EQ(2, ilp_hg.superVertexBlock(5));
}

} // namespace