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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/mt_kahypar.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/partitioner.h"

using ::testing::Test;

namespace mt_kahypar {

class APartitioner : public Test {

static size_t num_threads;

 public:
  APartitioner() :
    hypergraph(),
    context() {
    context.partition.graph_filename = "test_instances/ibm01.hgr";
    context.partition.graph_community_filename = "test_instances/ibm01.hgr.community";
    context.partition.mode = Mode::direct_kway;
    context.partition.objective = Objective::km1;
    context.partition.epsilon = 0.03;
    context.partition.k = 2;
    context.partition.verbose_output = false;

    // Coarsening
    context.coarsening.algorithm = CoarseningAlgorithm::community_coarsener;
    context.coarsening.contraction_limit_multiplier = 20;
    context.coarsening.max_allowed_weight_multiplier = 1;
    context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
    context.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
    context.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;

    // Initial Partitioning
    context.initial_partitioning.runs = 1;
    context.initial_partitioning.context_file = "test_instances/fast_initial_partitioning.ini";

    // Label Propagation
    context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::label_propagation_km1;
    context.refinement.label_propagation.execution_policy = ExecutionType::exponential;

    // Shared Memory
    context.shared_memory.num_threads = num_threads;
    context.shared_memory.use_community_redistribution = true;
    context.shared_memory.assignment_strategy = CommunityAssignmentStrategy::bin_packing;
    context.shared_memory.assignment_objective = CommunityAssignmentObjective::pin_objective;

    // Read hypergraph
    hypergraph = io::readHypergraphFile("test_instances/ibm01.hgr",
      context.partition.k, InitialHyperedgeDistribution::equally);
  }

  static void SetUpTestSuite() {
    TBBNumaArena::instance(num_threads);
  }

  static void TearDownTestSuite() {
    TBBNumaArena::instance().terminate();
  }

  Hypergraph hypergraph;
  Context context;
};

size_t APartitioner::num_threads = std::thread::hardware_concurrency();

void verifyThatHypergraphsAreEquivalent(const Hypergraph& hypergraph,
                                        const Hypergraph& reference) {

  // Verify equivallence of hypernodes and incident nets
  for ( const HypernodeID& hn : reference.nodes() ) {
    const HypernodeID original_id = reference.originalNodeID(hn);
    const HypernodeID u = hypergraph.globalNodeID(original_id);
    ASSERT_TRUE(hypergraph.nodeIsEnabled(u));

    std::set<HyperedgeID> incident_nets;
    for ( const HyperedgeID& he : reference.incidentEdges(hn) ) {
      const HyperedgeID original_edge_id = reference.originalEdgeID(he);
      incident_nets.insert(original_edge_id);
    }

    size_t num_incident_nets = 0;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(u) ) {
      const HyperedgeID original_edge_id = hypergraph.originalEdgeID(he);
      ASSERT_TRUE(incident_nets.find(original_edge_id) != incident_nets.end()) << V(u) << V(original_edge_id);
      ++num_incident_nets;
    }
    ASSERT_EQ(num_incident_nets, incident_nets.size());
  }

  // Verify equivallence of hyperedges and pins
  for ( const HyperedgeID& he : reference.edges() ) {
    const HyperedgeID original_id = reference.originalEdgeID(he);
    const HyperedgeID e = hypergraph.globalEdgeID(original_id);
    ASSERT_TRUE(hypergraph.edgeIsEnabled(e));

    std::set<HypernodeID> pins;
    for ( const HypernodeID& pin : reference.pins(he) ) {
      const HypernodeID original_pin_id = reference.originalNodeID(pin);
      pins.insert(original_pin_id);
    }

    size_t num_pins = 0;
    for ( const HypernodeID& pin : hypergraph.pins(e) ) {
      const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
      ASSERT_TRUE(pins.find(original_pin_id) != pins.end()) << V(e) << V(original_pin_id);
      ++num_pins;
    }
    ASSERT_EQ(num_pins, pins.size());
  }
}



TEST_F(APartitioner, AssignsEachVertexAPartID) {
  partition::Partitioner().partition(hypergraph, context);

  size_t num_nodes = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_NE(-1, hypergraph.partID(hn));
    ASSERT_LE(hypergraph.partID(hn), context.partition.k);
    ++num_nodes;
  }
  ASSERT_EQ(hypergraph.initialNumNodes(), num_nodes);
}

TEST_F(APartitioner, ComputesCorrectBlockWeightsAndPartSizes) {
  partition::Partitioner().partition(hypergraph, context);

  std::vector<HypernodeWeight> weights(hypergraph.k(), 0);
  std::vector<size_t> sizes(hypergraph.k(), 0);
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    weights[hypergraph.partID(hn)] += hypergraph.nodeWeight(hn);
    ++sizes[hypergraph.partID(hn)];
  }

  for ( PartitionID k = 0; k < hypergraph.k(); ++k ) {
    ASSERT_EQ(weights[k], hypergraph.partWeight(k));
    ASSERT_EQ(sizes[k], hypergraph.partSize(k));
  }
}

TEST_F(APartitioner, ComputesCorrectPinCountsInPartValuesAndConnectivitySets) {
  partition::Partitioner().partition(hypergraph, context);

  for ( const HyperedgeID& he : hypergraph.edges() ) {
    std::vector<HypernodeID> pin_count_in_part(context.partition.k, 0);
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      ++pin_count_in_part[hypergraph.partID(pin)];
    }

    PartitionID connectivity = 0;
    for ( PartitionID k = 0; k < context.partition.k; ++k ) {
      if ( pin_count_in_part[k] > 0 ) {
        ++connectivity;
      }
      ASSERT_EQ(pin_count_in_part[k], hypergraph.pinCountInPart(he, k));
    }

    ASSERT_EQ(connectivity, hypergraph.connectivity(he));
    size_t connectivity_2 = 0;
    for ( const PartitionID& id : hypergraph.connectivitySet(he) ) {
      ++connectivity_2;
      ASSERT_GT(pin_count_in_part[id], 0);
    }
    ASSERT_EQ(connectivity, connectivity_2);
  }
}

TEST_F(APartitioner, ComputesCorrectBorderNodes) {
  partition::Partitioner().partition(hypergraph, context);

  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    bool is_border_node = false;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      if ( hypergraph.connectivity(he) > 1 ) {
        is_border_node = true;
        break;
      }
    }
    ASSERT_EQ(is_border_node, hypergraph.isBorderNode(hn)) << V(hn);
  }
}

TEST_F(APartitioner, IsEqualWithInputHypergraphAfterPartitioning) {
  Hypergraph reference = io::readHypergraphFile("test_instances/ibm01.hgr",
    context.partition.k, InitialHyperedgeDistribution::equally);
  partition::Partitioner().partition(hypergraph, context);
  verifyThatHypergraphsAreEquivalent(hypergraph, reference);
}

} // namespace mt_kahypar
