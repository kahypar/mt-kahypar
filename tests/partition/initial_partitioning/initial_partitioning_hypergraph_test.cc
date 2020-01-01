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

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_hypergraph.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
class AInitialPartitioningHypergraph : public ds::AHypergraph<2> {
 private:
  using Base = ds::AHypergraph<2>;

 public:
  using Base::TypeTraits;
  using Base::TBBArena;
  using Base::TestHypergraph;
  using InitialPartitioningHypergraph = InitialPartitioningHypergraphT<TypeTraits>;

  AInitialPartitioningHypergraph() :
    Base(),
    hypergraph(construct_hypergraph(7,
                                    { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } },
                                    { 0, 0, 0, 1, 1, 1, 1 },
                                    { 0, 0, 1, 1 },
                                    { 0, 0, 1, 1, 2, 3, 2 })),
    context() {
    context.partition.k = 2;
    context.partition.epsilon = 0.2;
    context.partition.objective = kahypar::Objective::km1;
    // Max Part Weight = 4
    context.setupPartWeights(hypergraph.totalWeight());
    utils::Timer::instance().disable();
  }

  static void TearDownTestSuite() {
    TBBArena::instance().terminate();
  }

  TestHypergraph hypergraph;
  Context context;
};

TEST_F(AInitialPartitioningHypergraph, HasUnequalGlobalAndLocalHypergraphs) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  ASSERT_NE(&ip_hypergraph.global_hypergraph(), &ip_hypergraph.local_hypergraph());
}

TEST_F(AInitialPartitioningHypergraph, LocalAndGlobalHypergraphHaveSameNumberOfNodesEdgesAndPins) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  ASSERT_EQ(ip_hypergraph.global_hypergraph().initialNumNodes(),
            ip_hypergraph.local_hypergraph().initialNumNodes());
  ASSERT_EQ(ip_hypergraph.global_hypergraph().initialNumEdges(),
            ip_hypergraph.local_hypergraph().initialNumEdges());
  ASSERT_EQ(ip_hypergraph.global_hypergraph().initialNumPins(),
            ip_hypergraph.local_hypergraph().initialNumPins());
}

TEST_F(AInitialPartitioningHypergraph, ReturnsAnUnassignedLocalHypernode1) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();
  ASSERT_EQ(-1, local_hg.partID(ip_hypergraph.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningHypergraph, ReturnsAnUnassignedLocalHypernode2) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 2;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_hypergraph.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningHypergraph, ReturnsAnUnassignedLocalHypernode3) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 4;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_hypergraph.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningHypergraph, ReturnsAnUnassignedLocalHypernode4) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 6;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_hypergraph.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningHypergraph, ReturnsInvalidHypernodeIfAllHypernodesAreAssigned) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();

  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
  }

  ASSERT_EQ(std::numeric_limits<HypernodeID>::max(),
            ip_hypergraph.get_unassigned_hypernode());
}

TEST_F(AInitialPartitioningHypergraph, ReturnsValidUnassignedHypernodeIfPartitionIsResetted) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
  ip_hypergraph.reset_unassigned_hypernodes();

  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
  }

  ASSERT_EQ(std::numeric_limits<HypernodeID>::max(),
            ip_hypergraph.get_unassigned_hypernode());

  local_hg.resetPartition();
  ip_hypergraph.reset_unassigned_hypernodes();
  ASSERT_EQ(-1, local_hg.partID(ip_hypergraph.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningHypergraph, AppliesPartitionToHypergraph) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();

  // Cut = 2
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionToHypergraph) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();

  // Cut = 3
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  // Cut = 2
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionWithImbalancedPartitionToHypergraph1) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();

  // Cut = 1, but imbalanced
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  // Cut = 2
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionWithImbalancedPartitionToHypergraph2) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();

  // Cut = 1, but imbalanced
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  // Cut = 2, also imbalanced, but has better balance than first cut
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(1, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionWithImbalancedPartitionToHypergraph3) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();

  // Cut = 2, but imbalanced
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  // Cut = 1, also imbalanced but better metric
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
  local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
  ip_hypergraph.commit();

  ip_hypergraph.apply();

  ASSERT_EQ(1, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(1, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(1, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionToHypergraphInParallel1) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 3
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 2
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  });

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionToHypergraphInParallel2) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 1, but imbalanced
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 2, but balanced
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  });

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(0, hypergraph.partID(id[2]));
  ASSERT_EQ(1, hypergraph.partID(id[3]));
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningHypergraph, AppliesBestPartitionToHypergraphInParallel3) {
  InitialPartitioningHypergraph ip_hypergraph(
    hypergraph, context, TBBArena::GLOBAL_TASK_GROUP);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5),
                                  GLOBAL_ID(hypergraph, 6) };

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 3
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    TestHypergraph& local_hg = ip_hypergraph.local_hypergraph();
    // Cut = 2
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[0]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[1]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[2]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[3]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[4]), 0);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[5]), 1);
    local_hg.setNodePart(ip_hypergraph.map_hypernode_to_local_hypergraph(id[6]), 1);
    ip_hypergraph.commit();
  });

  ip_hypergraph.apply();

  ASSERT_EQ(2, metrics::objective(hypergraph, context.partition.objective));
  ASSERT_EQ(0, hypergraph.partID(id[0]));
  ASSERT_EQ(0, hypergraph.partID(id[1]));
  ASSERT_EQ(1, hypergraph.partID(id[2]));
  ASSERT_EQ(0, hypergraph.partID(id[3]));
  ASSERT_EQ(0, hypergraph.partID(id[4]));
  ASSERT_EQ(1, hypergraph.partID(id[5]));
  ASSERT_EQ(1, hypergraph.partID(id[6]));
}


}  // namespace mt_kahypar
