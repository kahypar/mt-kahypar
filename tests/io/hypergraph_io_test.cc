/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
    hypergraph = readHypergraphFile(filename);
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

#ifdef USE_GRAPH_PARTITIONER
  void verifyIDs(const std::vector<std::pair<HyperedgeID, HyperedgeID>>& edge_pairs, HyperedgeID expected_max) {
    HyperedgeID max_id = 0;
    for (auto [forward, backward] : edge_pairs) {
      HyperedgeID f_id = hypergraph.uniqueEdgeID(forward);
      HyperedgeID b_id = hypergraph.uniqueEdgeID(backward);
      ASSERT_EQ(f_id, b_id);
      max_id = std::max(max_id, std::max(f_id, b_id));
    }
    ASSERT_EQ(expected_max, max_id);
  }
#endif

  Hypergraph hypergraph;
};

#ifndef USE_GRAPH_PARTITIONER
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

TEST_F(AHypergraphReader, ReadsAMetisGraph) {
  this->hypergraph = readGraphFile("../tests/instances/unweighted_graph.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 2, 3, 4 }, { 1, 3, 5, 6 }, { 4, 6, 7, 8 },
      { 0, 5, 9 }, { 7, 9, 10 }, { 8, 10 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 },
    { 1, 3 }, { 2, 4 }, { 2, 3 }, { 3, 5 },
    { 3, 6 }, { 4, 5 }, { 5, 6 }
  });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithNodeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_node_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 2, 3, 4 }, { 1, 3, 5, 6 }, { 4, 6, 7, 8 },
      { 0, 5, 9 }, { 7, 9, 10 }, { 8, 10 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 },
    { 1, 3 }, { 2, 4 }, { 2, 3 }, { 3, 5 },
    { 3, 6 }, { 4, 5 }, { 5, 6 }
  });

  // Verify Node Weights
  ASSERT_EQ(4, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(5, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(6, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithEdgeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_edge_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 2, 3, 4 }, { 1, 3, 5, 6 }, { 4, 6, 7, 8 },
      { 0, 5, 9 }, { 7, 9, 10 }, { 8, 10 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 },
    { 1, 3 }, { 2, 4 }, { 2, 3 }, { 3, 5 },
    { 3, 6 }, { 4, 5 }, { 5, 6 }
  });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 2, 4} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {1, 3, 6, 7, 9} ) {
    ASSERT_EQ(2, this->hypergraph.edgeWeight(e));
  }
  ASSERT_EQ(3, this->hypergraph.edgeWeight(5));
  ASSERT_EQ(5, this->hypergraph.edgeWeight(8));
  ASSERT_EQ(6, this->hypergraph.edgeWeight(10));
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithNodeAndEdgeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_node_and_edge_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 2, 3, 4 }, { 1, 3, 5, 6 }, { 4, 6, 7, 8 },
      { 0, 5, 9 }, { 7, 9, 10 }, { 8, 10 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 },
    { 1, 3 }, { 2, 4 }, { 2, 3 }, { 3, 5 },
    { 3, 6 }, { 4, 5 }, { 5, 6 }
  });

  // Verify Node Weights
  ASSERT_EQ(4, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(5, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(6, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 2, 4} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {1, 3, 6, 7, 9} ) {
    ASSERT_EQ(2, this->hypergraph.edgeWeight(e));
  }
  ASSERT_EQ(3, this->hypergraph.edgeWeight(5));
  ASSERT_EQ(5, this->hypergraph.edgeWeight(8));
  ASSERT_EQ(6, this->hypergraph.edgeWeight(10));
}
#endif

#ifdef USE_GRAPH_PARTITIONER
TEST_F(AHypergraphReader, ReadsAMetisGraph) {
  this->hypergraph = readGraphFile("../tests/instances/unweighted_graph.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8, 9 }, { 10, 11, 12, 13 },
      { 14, 15, 16 }, { 17, 18, 19 }, { 20, 21 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 0 },
    { 1, 2 }, { 1, 3 }, { 2, 0 }, { 2, 1 },
    { 2, 3 }, { 2, 4 }, { 3, 1 }, { 3, 2 },
    { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 },
    { 4, 5 }, { 5, 3 }, { 5, 4 }, { 5, 6 },
    { 6, 3 }, { 6, 5 }
  });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }

  // Verify IDs
  this->verifyIDs({
      {0, 3}, {1, 6}, {2, 14}, {4, 7},
      {5, 10}, {8, 11}, {9, 15}, {12, 17},
      {13, 20}, {16, 18}, {19, 21},
    }, 10
  );
}

TEST_F(AHypergraphReader, ReadsAMetisGraphNoNewline) {
  this->hypergraph = readGraphFile("../tests/instances/graph_no_newline.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8, 9 }, { 10, 11, 12, 13 },
      { 14, 15, 16 }, { 17, 18, 19 }, { 20, 21 } });

  // Verify Pins
  this->verifyPins({
    { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 0 },
    { 1, 2 }, { 1, 3 }, { 2, 0 }, { 2, 1 },
    { 2, 3 }, { 2, 4 }, { 3, 1 }, { 3, 2 },
    { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 },
    { 4, 5 }, { 5, 3 }, { 5, 4 }, { 5, 6 },
    { 6, 3 }, { 6, 5 }
  });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }

  // Verify IDs
  this->verifyIDs({
      {0, 3}, {1, 6}, {2, 14}, {4, 7},
      {5, 10}, {8, 11}, {9, 15}, {12, 17},
      {13, 20}, {16, 18}, {19, 21},
    }, 10
  );
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithNodeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_node_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8, 9 }, { 10, 11, 12, 13 },
      { 14, 15, 16 }, { 17, 18, 19 }, { 20, 21 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 0 },
    { 1, 2 }, { 1, 3 }, { 2, 0 }, { 2, 1 },
    { 2, 3 }, { 2, 4 }, { 3, 1 }, { 3, 2 },
    { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 },
    { 4, 5 }, { 5, 3 }, { 5, 4 }, { 5, 6 },
    { 6, 3 }, { 6, 5 }
  });

  // Verify Node Weights
  ASSERT_EQ(4, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(5, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(6, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }

  // Verify IDs
  this->verifyIDs({
      {0, 3}, {1, 6}, {2, 14}, {4, 7},
      {5, 10}, {8, 11}, {9, 15}, {12, 17},
      {13, 20}, {16, 18}, {19, 21},
    }, 10
  );
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithEdgeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_edge_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8, 9 }, { 10, 11, 12, 13 },
      { 14, 15, 16 }, { 17, 18, 19 }, { 20, 21 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 0 },
    { 1, 2 }, { 1, 3 }, { 2, 0 }, { 2, 1 },
    { 2, 3 }, { 2, 4 }, { 3, 1 }, { 3, 2 },
    { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 },
    { 4, 5 }, { 5, 3 }, { 5, 4 }, { 5, 6 },
    { 6, 3 }, { 6, 5 }
  });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 3, 2, 14, 5, 10} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {1, 6, 4, 7, 8, 11, 12, 17, 16, 18} ) {
    ASSERT_EQ(2, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {9, 15} ) {
    ASSERT_EQ(3, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {13, 20} ) {
    ASSERT_EQ(5, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {19, 21} ) {
    ASSERT_EQ(6, this->hypergraph.edgeWeight(e));
  }

  // Verify IDs
  this->verifyIDs({
      {0, 3}, {1, 6}, {2, 14}, {4, 7},
      {5, 10}, {8, 11}, {9, 15}, {12, 17},
      {13, 20}, {16, 18}, {19, 21},
    }, 10
  );
}

TEST_F(AHypergraphReader, ReadsAMetisGraphWithNodeAndEdgeWeights) {
  this->hypergraph = readGraphFile("../tests/instances/graph_with_node_and_edge_weights.graph", true);

  // Verify Incident Edges
  this->verifyIncidentNets(
    { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8, 9 }, { 10, 11, 12, 13 },
      { 14, 15, 16 }, { 17, 18, 19 }, { 20, 21 }, {} });

  // Verify Pins
  this->verifyPins({
    { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 0 },
    { 1, 2 }, { 1, 3 }, { 2, 0 }, { 2, 1 },
    { 2, 3 }, { 2, 4 }, { 3, 1 }, { 3, 2 },
    { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 },
    { 4, 5 }, { 5, 3 }, { 5, 4 }, { 5, 6 },
    { 6, 3 }, { 6, 5 }
  });

  // Verify Node Weights
  ASSERT_EQ(4, this->hypergraph.nodeWeight(0));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(1));
  ASSERT_EQ(5, this->hypergraph.nodeWeight(2));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(3));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(4));
  ASSERT_EQ(6, this->hypergraph.nodeWeight(5));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(6));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(7));

  // Verify Edge Weights
  for ( HyperedgeID e : {0, 3, 2, 14, 5, 10} ) {
    ASSERT_EQ(1, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {1, 6, 4, 7, 8, 11, 12, 17, 16, 18} ) {
    ASSERT_EQ(2, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {9, 15} ) {
    ASSERT_EQ(3, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {13, 20} ) {
    ASSERT_EQ(5, this->hypergraph.edgeWeight(e));
  }
  for ( HyperedgeID e : {19, 21} ) {
    ASSERT_EQ(6, this->hypergraph.edgeWeight(e));
  }

  // Verify IDs
  this->verifyIDs({
      {0, 3}, {1, 6}, {2, 14}, {4, 7},
      {5, 10}, {8, 11}, {9, 15}, {12, 17},
      {13, 20}, {16, 18}, {19, 21},
    }, 10
  );
}
#endif

}  // namespace io
}  // namespace mt_kahypar
