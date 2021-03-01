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

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"

using ::testing::Test;

namespace mt_kahypar {

template< PartitionID k,
          kahypar::Objective objective,
          bool use_community_detection,
          InitialPartitioningMode initial_partitioning_mode,
          LabelPropagationAlgorithm lp_algorithm >
struct TestConfig {
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = objective;
  static constexpr bool USE_COMMUNITY_DETECTION = use_community_detection;
  static constexpr InitialPartitioningMode INITIAL_PARTITIONING_MODE = initial_partitioning_mode;
  static constexpr LabelPropagationAlgorithm LP_ALGORITHM = lp_algorithm;
};

template< typename Config >
class MtKaHyPar : public Test {
  static size_t num_threads;

 public:
  MtKaHyPar() :
    context() {

    if ( context.partition.paradigm == Paradigm::multilevel ) {
      parseIniToContext(context, "../config/fast_preset.ini");
    } else {
      parseIniToContext(context, "../config/strong_preset.ini");
    }

    context.partition.graph_filename = "../tests/instances/ibm01.hgr";
    context.partition.k = Config::K;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.epsilon = 0.03;
    context.partition.seed = 42;
    context.partition.verbose_output = true;
    context.partition.show_detailed_timings = true;
    context.preprocessing.use_community_detection = Config::USE_COMMUNITY_DETECTION;
    context.initial_partitioning.mode = Config::INITIAL_PARTITIONING_MODE;
    context.initial_partitioning.refinement.label_propagation.algorithm = Config::LP_ALGORITHM;
    context.refinement.label_propagation.algorithm = Config::LP_ALGORITHM;
    context.shared_memory.num_threads = num_threads;
  }

  static void SetUpTestSuite() {
    TBBNumaArena::instance(num_threads);
  }

  Context context;
};

template< typename Config >
size_t MtKaHyPar<Config>::num_threads = HardwareTopology::instance().num_cpus();

void verifyThatHypergraphsAreEquivalent(const Hypergraph& hypergraph,
                                        const Hypergraph& reference) {
  // Verify equivalence of hypernodes and incident nets
  for (const HypernodeID& hn : reference.nodes()) {
    ASSERT_TRUE(hypergraph.nodeIsEnabled(hn));

    std::set<HyperedgeID> incident_nets;
    for (const HyperedgeID& he : reference.incidentEdges(hn)) {
      incident_nets.insert(he);
    }

    size_t num_incident_nets = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      ASSERT_TRUE(incident_nets.find(he) != incident_nets.end()) << V(hn) << V(he);
      ++num_incident_nets;
    }
    ASSERT_EQ(num_incident_nets, incident_nets.size());
  }

  // Verify equivalence of hyperedges and pins
  for (const HyperedgeID& he : reference.edges()) {
    ASSERT_TRUE(hypergraph.edgeIsEnabled(he));

    std::set<HypernodeID> pins;
    for (const HypernodeID& pin : reference.pins(he)) {
      pins.insert(pin);
    }

    size_t num_pins = 0;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      ASSERT_TRUE(pins.find(pin) != pins.end()) << V(he) << V(pin);
      ++num_pins;
    }
    ASSERT_EQ(num_pins, pins.size());
  }
}

template< PartitionID k,
          bool use_community_detection,
          InitialPartitioningMode initial_partitioning_mode>
using MultiLevelKm1Config = TestConfig<k,
                                       kahypar::Objective::km1,
                                       use_community_detection,
                                       initial_partitioning_mode,
                                       LabelPropagationAlgorithm::label_propagation_km1>;

typedef ::testing::Types<MultiLevelKm1Config<2,
                                             false,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<4,
                                             false,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<8,
                                             false,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<2,
                                             false,
                                             InitialPartitioningMode::recursive_bisection>,
                         MultiLevelKm1Config<4,
                                             false,
                                             InitialPartitioningMode::recursive_bisection>,
                         MultiLevelKm1Config<8,
                                             false,
                                             InitialPartitioningMode::recursive_bisection>,
                         MultiLevelKm1Config<2,
                                             true,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<4,
                                             true,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<8,
                                             true,
                                             InitialPartitioningMode::recursive>,
                         MultiLevelKm1Config<2,
                                             true,
                                             InitialPartitioningMode::recursive_bisection>,
                         MultiLevelKm1Config<4,
                                             true,
                                             InitialPartitioningMode::recursive_bisection>,
                         MultiLevelKm1Config<8,
                                             true,
                                             InitialPartitioningMode::recursive_bisection> > TestConfigs;

TYPED_TEST_CASE(MtKaHyPar, TestConfigs);

void partitionHypergraph(Hypergraph& hypergraph, Context& context) {
  if (context.partition.objective != kahypar::Objective::km1) {
    context.refinement.fm.algorithm = mt_kahypar::FMAlgorithm::do_nothing;
    context.refinement.global_fm.use_global_fm = false;

    context.initial_partitioning.refinement.fm.algorithm = mt_kahypar::FMAlgorithm::do_nothing;
    context.initial_partitioning.refinement.global_fm.use_global_fm = false;
  }
  mt_kahypar::register_memory_pool(hypergraph, context);

  // Partition Hypergraph
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  PartitionedHypergraph partitioned_hypergraph = partition(hypergraph, context);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds(end - start);
  mt_kahypar::io::printPartitioningResults(partitioned_hypergraph, context, elapsed_seconds);

  // Verify that partitioned hypergraph is
  // equivalent with input hypergraph
  Hypergraph reference = io::readHypergraphFile(
    context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);
  verifyThatHypergraphsAreEquivalent(hypergraph, reference);

  HypernodeID num_hypernodes = 0;
  std::vector<HypernodeWeight> weights(partitioned_hypergraph.k(), 0);
  std::vector<size_t> sizes(partitioned_hypergraph.k(), 0);
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    // Make sure that each hypernode has a block id between [0,k)
    PartitionID part_id = partitioned_hypergraph.partID(hn);
    ASSERT_NE(-1, part_id);
    ASSERT_LE(part_id, context.partition.k);
    weights[part_id] += partitioned_hypergraph.nodeWeight(hn);
    ++sizes[part_id];
    ++num_hypernodes;

    // Verify that border nodes are stored correctly
    bool is_border_node = false;
    for ( const HyperedgeID& he : partitioned_hypergraph.incidentEdges(hn) ) {
      if ( partitioned_hypergraph.connectivity(he) > 1 ) {
        is_border_node = true;
        break;
      }
    }
    ASSERT_EQ(is_border_node, partitioned_hypergraph.isBorderNode(hn));
  }
  // Verify that hypergraph is fully uncontracted
  ASSERT_EQ(partitioned_hypergraph.initialNumNodes(), num_hypernodes);

  // Verify block weights and sizes
  for ( PartitionID k = 0; k < partitioned_hypergraph.k(); ++k ) {
    ASSERT_EQ(weights[k], partitioned_hypergraph.partWeight(k));
  }

  // Verify connectivity (sets) and pin count in part
  size_t num_hyperedges = 0;
  for ( const HyperedgeID& he : partitioned_hypergraph.edges() ) {
    std::vector<HypernodeID> pin_count_in_part(partitioned_hypergraph.k(), 0);
    PartitionID connectivity = 0;
    for ( const HypernodeID& pin : partitioned_hypergraph.pins(he) ) {
      PartitionID part_id = partitioned_hypergraph.partID(pin);
      if ( pin_count_in_part[part_id] == 0 ) {
        ++connectivity;
      }
      ++pin_count_in_part[part_id];
    }
    ASSERT_EQ(connectivity, partitioned_hypergraph.connectivity(he));

    for ( PartitionID k = 0; k < partitioned_hypergraph.k(); ++k ) {
      ASSERT_EQ(pin_count_in_part[k], partitioned_hypergraph.pinCountInPart(he, k));
    }

    PartitionID connectivity_2 = 0;
    for ( const PartitionID& k : partitioned_hypergraph.connectivitySet(he) ) {
      ++connectivity_2;
      ASSERT_GT(pin_count_in_part[k], 0);
    }
    ASSERT_EQ(connectivity, connectivity_2);

    ++num_hyperedges;
  }
  ASSERT_EQ(hypergraph.initialNumEdges(), num_hyperedges);

  mt_kahypar::parallel::MemoryPool::instance().free_memory_chunks();
}

TYPED_TEST(MtKaHyPar, PartitionsAVLSIInstance) {
  // Read Hypergraph
  this->context.partition.graph_filename = "../tests/instances/ibm01.hgr";
  Hypergraph hypergraph = io::readHypergraphFile(
    this->context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);

  partitionHypergraph(hypergraph, this->context);
}

TYPED_TEST(MtKaHyPar, PartitionsASparseMatrixInstance) {
  // Read Hypergraph
  this->context.partition.graph_filename = "../tests/instances/powersim.mtx.hgr";
  Hypergraph hypergraph = io::readHypergraphFile(
    this->context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);

  partitionHypergraph(hypergraph, this->context);
}

TYPED_TEST(MtKaHyPar, PartitionsASATInstance) {
  // Read Hypergraph
  this->context.partition.graph_filename = "../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr";
  Hypergraph hypergraph = io::readHypergraphFile(
    this->context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);

  partitionHypergraph(hypergraph, this->context);
}

TYPED_TEST(MtKaHyPar, PartitionsAVLSIInstanceWithIndividualPartWeights) {
  // Read Hypergraph
  this->context.partition.graph_filename = "../tests/instances/ibm01.hgr";
  Hypergraph hypergraph = io::readHypergraphFile(
    this->context.partition.graph_filename, TBBNumaArena::GLOBAL_TASK_GROUP);

  // setup individual part weights
  HypernodeID k = this->context.partition.k;
  double max_part_weight_balanced = (1 + this->context.partition.epsilon) *
                                             ceil(hypergraph.totalWeight() / static_cast<double>(k));
  double weight_delta = 2 * max_part_weight_balanced / (k + 1);
  this->context.partition.use_individual_part_weights = true;
  this->context.partition.max_part_weights.clear();
  for ( HypernodeID n = 1; n <= k; ++n ) {
    this->context.partition.max_part_weights.push_back(n * weight_delta);
  }

  partitionHypergraph(hypergraph, this->context);
}


}  // namespace mt_kahypar
