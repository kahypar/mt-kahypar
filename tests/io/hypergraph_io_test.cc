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
#include "mt-kahypar/io/hypergraph_io.h"

using ::testing::Test;

namespace mt_kahypar {
namespace io {

class AHypergraphReader : public Test {

 public:
  AHypergraphReader() :
    hypergraph() { }

  void readHypergraph(const std::string& filename) {
    hypergraph = readHypergraphFile(
      filename, TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  void verifyIncidentNets(const std::vector< std::set<HyperedgeID> >& references) {
    ASSERT(hypergraph.initialNumNodes() == references.size());
    for (HypernodeID hn = 0; hn < hypergraph.initialNumNodes(); ++hn) {
      const std::set<HyperedgeID>& reference = references[hn];
      size_t count = 0;
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        ASSERT_TRUE(reference.find(he) != reference.end()) << V(hn) << V(he);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  void verifyPins(const std::vector< std::set<HypernodeID> >& references) {
    ASSERT(hypergraph.initialNumEdges() == references.size());
    for (HyperedgeID he = 0; he < hypergraph.initialNumEdges(); ++he) {
      const std::set<HypernodeID>& reference = references[he];
      size_t count = 0;
      for (const HypernodeID& pin : hypergraph.pins(he)) {
        ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  Hypergraph hypergraph;
};

TEST_F(AHypergraphReader, ReadsAnUnweightedHypergraph) {
  this->readHypergraph("../tests/instances/unweighted_hypergraph.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));

  // Verify Edge Weights
  ASSERT_EQ(1, this->hypergraph.edgeWeight(0));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(1));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(2));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(3));
}

TEST_F(AHypergraphReader, ReadsAnHypergraphWithEdgeWeights) {
  this->readHypergraph("../tests/instances/hypergraph_with_edge_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));

  // Verify Edge Weights
  ASSERT_EQ(4, this->hypergraph.edgeWeight(0));
  ASSERT_EQ(2, this->hypergraph.edgeWeight(1));
  ASSERT_EQ(3, this->hypergraph.edgeWeight(2));
  ASSERT_EQ(8, this->hypergraph.edgeWeight(3));
}

TEST_F(AHypergraphReader, ReadsAnHypergraphWithNodeWeights) {
  this->readHypergraph("../tests/instances/hypergraph_with_node_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(5, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(4, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(9, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(6));

  // Verify Edge Weights
  ASSERT_EQ(1, this->hypergraph.edgeWeight(0));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(1));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(2));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(3));
}

TEST_F(AHypergraphReader, ReadsAnHypergraphWithNodeAndEdgeWeights) {
  this->readHypergraph("../tests/instances/hypergraph_with_node_and_edge_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(5, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(4, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(9, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(6));

  // Verify Edge Weights
  ASSERT_EQ(4, this->hypergraph.edgeWeight(0));
  ASSERT_EQ(2, this->hypergraph.edgeWeight(1));
  ASSERT_EQ(3, this->hypergraph.edgeWeight(2));
  ASSERT_EQ(8, this->hypergraph.edgeWeight(3));
}
}  // namespace io
}  // namespace mt_kahypar
