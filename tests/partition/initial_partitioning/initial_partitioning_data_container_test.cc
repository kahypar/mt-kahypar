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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"

using ::testing::Test;

namespace mt_kahypar {

using TestTypeTraits = ds::TestTypeTraits<2>;
using HyperGraph = typename TestTypeTraits::HyperGraph;
using HyperGraphFactory = typename TestTypeTraits::HyperGraphFactory;
using PartitionedHyperGraph = typename TestTypeTraits::PartitionedHyperGraph<>;
using PartitionedHyperGraphWithoutBorderVertices = typename TestTypeTraits::PartitionedHyperGraph<false>;
using TBB = typename TestTypeTraits::TBB;

class AInitialPartitioningDataContainer :
  public ds::HypergraphFixture<HyperGraph, HyperGraphFactory, TBB> {
 private:
  using Base = ds::HypergraphFixture<HyperGraph, HyperGraphFactory, TBB>;

 public:
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TestTypeTraits>;

  AInitialPartitioningDataContainer() :
    Base(),
    context() {
    context.partition.k = 2;
    context.partition.epsilon = 0.2;
    context.partition.objective = kahypar::Objective::km1;
    // Max Part Weight = 4
    context.setupPartWeights(hypergraph.totalWeight());
    utils::Timer::instance().disable();
  }

  static void TearDownTestSuite() {
    TBB::instance().terminate();
  }

  using Base::hypergraph;
  using Base::id;
  Context context;
};

TEST_F(AInitialPartitioningDataContainer, ReturnsAnUnassignedLocalHypernode1) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();
  ASSERT_EQ(-1, local_hg.partID(ip_data.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningDataContainer, ReturnsAnUnassignedLocalHypernode2) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 2;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_data.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningDataContainer, ReturnsAnUnassignedLocalHypernode3) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 4;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_data.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningDataContainer, ReturnsAnUnassignedLocalHypernode4) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();

  size_t num_hypernodes_to_assign = 6;
  size_t assigned_hypernodes = 0;
  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
    ++assigned_hypernodes;
    if ( assigned_hypernodes == num_hypernodes_to_assign ) {
      break;
    }
  }

  ASSERT_EQ(-1, local_hg.partID(ip_data.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningDataContainer, ReturnsInvalidHypernodeIfAllHypernodesAreAssigned) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();

  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
  }

  ASSERT_EQ(std::numeric_limits<HypernodeID>::max(),
            ip_data.get_unassigned_hypernode());
}

TEST_F(AInitialPartitioningDataContainer, ReturnsValidUnassignedHypernodeIfPartitionIsResetted) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
  ip_data.reset_unassigned_hypernodes();

  for ( const HypernodeID& hn : local_hg.nodes() ) {
    local_hg.setNodePart(hn, 0);
  }

  ASSERT_EQ(std::numeric_limits<HypernodeID>::max(),
            ip_data.get_unassigned_hypernode());

  local_hg.resetPartition();
  ip_data.reset_unassigned_hypernodes();
  ASSERT_EQ(-1, local_hg.partID(ip_data.get_unassigned_hypernode()));
}

TEST_F(AInitialPartitioningDataContainer, AppliesPartitionToHypergraph) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();

  // Cut = 2
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 0);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionToHypergraph) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();

  // Cut = 3
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 0);
  local_hg.setNodePart(id[3], 0);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  // Cut = 2
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 0);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionWithImbalancedPartitionToHypergraph1) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();

  // Cut = 1, but imbalanced
  local_hg.setNodePart(id[0], 1);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 1);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  // Cut = 2
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 0);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionWithImbalancedPartitionToHypergraph2) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();

  // Cut = 1, but imbalanced
  local_hg.setNodePart(id[0], 1);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 1);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  // Cut = 2, also imbalanced but better balance
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 1);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionWithImbalancedPartitionToHypergraph3) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);
  PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();

  // Cut = 1
  local_hg.setNodePart(id[0], 1);
  local_hg.setNodePart(id[1], 0);
  local_hg.setNodePart(id[2], 1);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  // Cut = 2
  local_hg.setNodePart(id[0], 0);
  local_hg.setNodePart(id[1], 1);
  local_hg.setNodePart(id[2], 1);
  local_hg.setNodePart(id[3], 1);
  local_hg.setNodePart(id[4], 1);
  local_hg.setNodePart(id[5], 1);
  local_hg.setNodePart(id[6], 1);
  ip_data.commit(InitialPartitioningAlgorithm::random);

  ip_data.apply();

  ASSERT_EQ(1, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionToHypergraphInParallel1) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 3
    local_hg.setNodePart(id[0], 0);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 0);
    local_hg.setNodePart(id[3], 0);
    local_hg.setNodePart(id[4], 1);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 2
    local_hg.setNodePart(id[0], 0);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 0);
    local_hg.setNodePart(id[3], 1);
    local_hg.setNodePart(id[4], 1);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  });

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionToHypergraphInParallel2) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 1, but imbalanced
    local_hg.setNodePart(id[0], 1);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 1);
    local_hg.setNodePart(id[3], 1);
    local_hg.setNodePart(id[4], 1);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 2, but balanced
    local_hg.setNodePart(id[0], 0);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 0);
    local_hg.setNodePart(id[3], 1);
    local_hg.setNodePart(id[4], 1);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  });

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}

TEST_F(AInitialPartitioningDataContainer, AppliesBestPartitionToHypergraphInParallel3) {
  PartitionedHyperGraph partitioned_hypergraph(context.partition.k, hypergraph);
  InitialPartitioningDataContainer ip_data(
    partitioned_hypergraph, context, TBB::GLOBAL_TASK_GROUP);

  std::atomic<size_t> cnt(0);
  tbb::parallel_invoke([&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 3
    local_hg.setNodePart(id[0], 0);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 0);
    local_hg.setNodePart(id[3], 0);
    local_hg.setNodePart(id[4], 1);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  }, [&] {
    ++cnt;
    while(cnt < 2) { }
    PartitionedHyperGraphWithoutBorderVertices& local_hg = ip_data.local_partitioned_hypergraph();
    // Cut = 2
    local_hg.setNodePart(id[0], 0);
    local_hg.setNodePart(id[1], 0);
    local_hg.setNodePart(id[2], 1);
    local_hg.setNodePart(id[3], 0);
    local_hg.setNodePart(id[4], 0);
    local_hg.setNodePart(id[5], 1);
    local_hg.setNodePart(id[6], 1);
    ip_data.commit(InitialPartitioningAlgorithm::random);
  });

  ip_data.apply();

  ASSERT_EQ(2, metrics::objective(partitioned_hypergraph, context.partition.objective));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[0]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[1]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[2]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[3]));
  ASSERT_EQ(0, partitioned_hypergraph.partID(id[4]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[5]));
  ASSERT_EQ(1, partitioned_hypergraph.partID(id[6]));
}


}  // namespace mt_kahypar
