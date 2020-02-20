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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/community_redistributor.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/policies/community_assignment_objective.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

template< int NUM_NUMA_NODES >
class ARedistributor : public Test {

 public:
  using TypeTraits = TestTypeTraits<NUM_NUMA_NODES>;
  using HyperGraph = typename TypeTraits::HyperGraph;
  using Factory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using Redistributor = mt_kahypar::preprocessing::CommunityRedistributorT<TypeTraits>;

  ARedistributor() :
    hypergraph(Factory::construct(TBB::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    context(),
    id() {
    id.resize(7);
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      id[hypergraph.originalNodeID(hn)] = hn;
    }

    hypergraph.setCommunityID(id[0], 0);
    hypergraph.setCommunityID(id[1], 0);
    hypergraph.setCommunityID(id[2], 1);
    hypergraph.setCommunityID(id[3], 1);
    hypergraph.setCommunityID(id[4], 2);
    hypergraph.setCommunityID(id[5], 3);
    hypergraph.setCommunityID(id[6], 2);
    hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);
  }

  static void SetUpTestSuite() {
    TBB::instance(HwTopology::instance().num_cpus());
  }

  static void TearDownTestSuite() {
    TBB::instance().terminate();
  }

  HyperGraph hypergraph;
  Context context;
  std::vector<HypernodeID> id;
};

using ARedistributorOnTwoNumaNodes = ARedistributor<2>;

TEST_F(ARedistributorOnTwoNumaNodes, RedistributesCommunities) {
  HyperGraph r_hypergraph =
    Redistributor::redistribute(TBB::GLOBAL_TASK_GROUP, hypergraph, { 0, 0, 1, 1 });

  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(1)));
  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(2)));
  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(3)));
  ASSERT_EQ(1, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(5)));
  ASSERT_EQ(1, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(6)));
}

#define SYSTEM_HAS_MORE_THAN_FOUR_CORES false

#if SYSTEM_HAS_NORE_THAN_FOUR_CORES

using ARedistributorOnFourNumaNodes = ARedistributor<4>;

TEST_F(ARedistributorOnFourNumaNodes, RedistributesCommunities) {
  HyperGraph r_hypergraph =
    Redistributor::redistribute(TBB::GLOBAL_TASK_GROUP, hypergraph, { 0, 1, 2, 3 });

  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(1)));
  ASSERT_EQ(1, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(2)));
  ASSERT_EQ(1, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(3)));
  ASSERT_EQ(2, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(4)));
  ASSERT_EQ(3, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, common::get_numa_node_of_vertex(r_hypergraph.globalNodeID(6)));
}
#endif
}  // namespace ds
}  // namespace mt_kahypar
