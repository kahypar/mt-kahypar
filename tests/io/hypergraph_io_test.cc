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
#include "mt-kahypar/io/tmp_hypergraph_io.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace tmp_io {

template< typename HyperGraph,
          typename HyperGraphFactory,
          typename TBBArena>
struct HypergraphTypeTraits {
  using Hypergraph = HyperGraph;
  using HypergraphFactory = HyperGraphFactory;
  using TBB = TBBArena;
};

template<typename TypeTraits>
class AHypergraphReader : public Test {

 using Hypergraph = typename TypeTraits::Hypergraph;
 using HypergraphFactory = typename TypeTraits::HypergraphFactory;
 using TBB = typename TypeTraits::TBB;

 public:
  AHypergraphReader() :
    hypergraph(),
    node_id(),
    edge_id() { }

  static void SetUpTestSuite() {
    TBB::instance(HardwareTopology::instance().num_cpus());
  }

  void readHypergraph(const std::string& filename) {
    hypergraph = readHypergraphFile<Hypergraph, HypergraphFactory>(
      filename, TBB::GLOBAL_TASK_GROUP);

    node_id.resize(hypergraph.initialNumNodes());
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      node_id[hypergraph.originalNodeID(hn)] = hn;
    }

    edge_id.resize(hypergraph.initialNumEdges());
    for ( const HyperedgeID& he : hypergraph.edges() ) {
      edge_id[hypergraph.originalEdgeID(he)] = he;
    }
  }

  void verifyIncidentNets(const std::vector< std::set<HyperedgeID> >& references) {
    ASSERT(hypergraph.initialNumNodes() == references.size());
    for (HypernodeID id = 0; id < hypergraph.initialNumNodes(); ++id) {
      const HypernodeID hn = hypergraph.globalNodeID(id);
      const std::set<HyperedgeID>& reference = references[id];
      size_t count = 0;
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        const HyperedgeID original_he = hypergraph.originalEdgeID(he);
        ASSERT_TRUE(reference.find(original_he) != reference.end()) << V(hn) << V(he);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  void verifyPins(const std::vector< std::set<HypernodeID> >& references) {
    ASSERT(hypergraph.initialNumEdges() == references.size());
    for (HyperedgeID id = 0; id < hypergraph.initialNumEdges(); ++id) {
      const HyperedgeID he = hypergraph.globalEdgeID(id);
      const std::set<HypernodeID>& reference = references[id];
      size_t count = 0;
      for (const HypernodeID& pin : hypergraph.pins(he)) {
        const HypernodeID original_pin = hypergraph.originalNodeID(pin);
        ASSERT_TRUE(reference.find(original_pin) != reference.end()) << V(he) << V(pin);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
  }

  Hypergraph hypergraph;
  std::vector<HypernodeID> node_id;
  std::vector<HypernodeID> edge_id;
};

// Mocking Numa Architecture (=> 2 NUMA Nodes)
using TypeTraits = ds::TestTypeTraits<2>;
using HwTopology = typename TypeTraits::HwTopology;
using TBB = typename TypeTraits::TBB;

// Define NUMA Hypergraph and Factory
using StaticHypergraph = ds::StaticHypergraph;
using StaticHypergraphFactory = ds::StaticHypergraphFactory;
using NumaHypergraph = ds::NumaHypergraph<StaticHypergraph, HwTopology, TBB>;
using NumaHypergraphFactory = ds::NumaHypergraphFactory<
  StaticHypergraph, StaticHypergraphFactory, HwTopology, TBB>;

typedef ::testing::Types<HypergraphTypeTraits<
                          StaticHypergraph,
                          StaticHypergraphFactory,
                          TBBNumaArena>,
                        HypergraphTypeTraits<
                          NumaHypergraph,
                          NumaHypergraphFactory,
                          TBB>> HypergraphTestTypes;


TYPED_TEST_CASE(AHypergraphReader, HypergraphTestTypes);

TYPED_TEST(AHypergraphReader, ReadsAnUnweightedHypergraph) {
  this->readHypergraph("test_instances/unweighted_hypergraph.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[0]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[1]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[2]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[3]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[4]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[5]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[6]));

  // Verify Edge Weights
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[0]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[1]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[2]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[3]));
}

TYPED_TEST(AHypergraphReader, ReadsAnHypergraphWithEdgeWeights) {
  this->readHypergraph("test_instances/hypergraph_with_edge_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[0]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[1]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[2]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[3]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[4]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[5]));
  ASSERT_EQ(1, this->hypergraph.nodeWeight(this->node_id[6]));

  // Verify Edge Weights
  ASSERT_EQ(4, this->hypergraph.edgeWeight(this->edge_id[0]));
  ASSERT_EQ(2, this->hypergraph.edgeWeight(this->edge_id[1]));
  ASSERT_EQ(3, this->hypergraph.edgeWeight(this->edge_id[2]));
  ASSERT_EQ(8, this->hypergraph.edgeWeight(this->edge_id[3]));
}

TYPED_TEST(AHypergraphReader, ReadsAnHypergraphWithNodeWeights) {
  this->readHypergraph("test_instances/hypergraph_with_node_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(5, this->hypergraph.nodeWeight(this->node_id[0]));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(this->node_id[1]));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(this->node_id[2]));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(this->node_id[3]));
  ASSERT_EQ(4, this->hypergraph.nodeWeight(this->node_id[4]));
  ASSERT_EQ(9, this->hypergraph.nodeWeight(this->node_id[5]));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(this->node_id[6]));

  // Verify Edge Weights
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[0]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[1]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[2]));
  ASSERT_EQ(1, this->hypergraph.edgeWeight(this->edge_id[3]));
}

TYPED_TEST(AHypergraphReader, ReadsAnHypergraphWithNodeAndEdgeWeights) {
  this->readHypergraph("test_instances/hypergraph_with_node_and_edge_weights.hgr");

  // Verify Incident Nets
  this->verifyIncidentNets(
    { { 0, 1 }, { 1 }, { 0, 3 }, { 1, 2 },
      {1, 2}, { 3 }, { 2, 3 } });

  // Verify Pins
  this->verifyPins({ { 0, 2 }, { 0, 1, 3, 4 },
    { 3, 4, 6 }, { 2, 5, 6 } });

  // Verify Node Weights
  ASSERT_EQ(5, this->hypergraph.nodeWeight(this->node_id[0]));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(this->node_id[1]));
  ASSERT_EQ(2, this->hypergraph.nodeWeight(this->node_id[2]));
  ASSERT_EQ(3, this->hypergraph.nodeWeight(this->node_id[3]));
  ASSERT_EQ(4, this->hypergraph.nodeWeight(this->node_id[4]));
  ASSERT_EQ(9, this->hypergraph.nodeWeight(this->node_id[5]));
  ASSERT_EQ(8, this->hypergraph.nodeWeight(this->node_id[6]));

  // Verify Edge Weights
  ASSERT_EQ(4, this->hypergraph.edgeWeight(this->edge_id[0]));
  ASSERT_EQ(2, this->hypergraph.edgeWeight(this->edge_id[1]));
  ASSERT_EQ(3, this->hypergraph.edgeWeight(this->edge_id[2]));
  ASSERT_EQ(8, this->hypergraph.edgeWeight(this->edge_id[3]));
}
}  // namespace tmp_io
}  // namespace mt_kahypar
