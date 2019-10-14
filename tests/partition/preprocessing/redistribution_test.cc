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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/community_redistributor.h"
#include "mt-kahypar/partition/preprocessing/policies/community_assignment_objective.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

class ARedistributorOnTwoNumaNodes : public AHypergraph<2> {

 private:
  using Base = AHypergraph<2>;

 public:
  using Base::TBBArena;
  using Base::TestStreamingHypergraph;
  using Base::TestHypergraph;
  using Redistributor = mt_kahypar::preprocessing::CommunityRedistributorT<TestTypeTraits<2>>;

  ARedistributorOnTwoNumaNodes() :
    Base(),
    hypergraph(construct_hypergraph( 7,
      { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
      { 0, 1, 0, 1, 0, 1, 0 },
      { 0, 0, 1, 1 },
      { 0, 0, 1, 1, 2, 3, 2 } )),
    context() { }

  static void TearDownTestSuite() {
    TBBArena::instance().terminate();
  }

  TestHypergraph hypergraph;
  Context context;
};

TEST_F(ARedistributorOnTwoNumaNodes, RedistributesCommunities) {
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(0)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(1)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(2)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(3)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(5)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(6)));

  TestHypergraph r_hypergraph = Redistributor::redistribute(hypergraph, 2, {0, 0, 1, 1});

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(1)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(2)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(3)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(5)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(6)));
}

class ARedistributorOnFourNumaNodes : public AHypergraph<4> {

 private:
  using Base = AHypergraph<4>;

 public:
  using Base::TestStreamingHypergraph;
  using Base::TestHypergraph;
  using Redistributor = mt_kahypar::preprocessing::CommunityRedistributorT<TestTypeTraits<4>>;

  ARedistributorOnFourNumaNodes() :
    Base(),
    hypergraph(construct_hypergraph( 7,
      { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
      { 0, 1, 2, 3, 0, 1, 2 },
      { 0, 1, 2, 3 },
      { 0, 0, 1, 1, 2, 3, 2 } )),
    context() { }

  TestHypergraph hypergraph;
  Context context;
};

TEST_F(ARedistributorOnFourNumaNodes, RedistributesCommunities) {
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(0)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(1)));
  ASSERT_EQ(2, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(2)));
  ASSERT_EQ(3, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(3)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, TestStreamingHypergraph::get_numa_node_of_vertex(hypergraph.globalNodeID(6)));

  TestHypergraph r_hypergraph = Redistributor::redistribute(hypergraph, 2, {0, 1, 2, 3});

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(1)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(2)));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(3)));
  ASSERT_EQ(2, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(4)));
  ASSERT_EQ(3, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, TestStreamingHypergraph::get_numa_node_of_vertex(r_hypergraph.globalNodeID(6)));
}

} // namespace ds
} // namespace mt_kahypar
