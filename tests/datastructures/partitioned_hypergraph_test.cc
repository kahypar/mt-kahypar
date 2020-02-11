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

#include <atomic>

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

template< typename PartitionedHG,
          typename HG,
          typename HGFactory,
          typename TBBArena>
struct PartitionedHypergraphTypeTraits {
  using PartitionedHyperGraph = PartitionedHG;
  using Hypergraph = HG;
  using Factory = HGFactory;
  using TBB = TBBArena;
};

template<typename TypeTraits>
class APartitionedHypergraph : public Test {

 using PartitionedHyperGraph = typename TypeTraits::PartitionedHyperGraph;
 using Hypergraph = typename TypeTraits::Hypergraph;
 using Factory = typename TypeTraits::Factory;
 using TBB = typename TypeTraits::TBB;

 public:
  APartitionedHypergraph() :
    hypergraph(Factory::construct(TBB::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    partitioned_hypergraph(3, hypergraph),
    id() {
    id.resize(7);
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      id[hypergraph.originalNodeID(hn)] = hn;
    }

    partitioned_hypergraph.setNodePart(id[0], 0);
    partitioned_hypergraph.setNodePart(id[1], 0);
    partitioned_hypergraph.setNodePart(id[2], 0);
    partitioned_hypergraph.setNodePart(id[3], 1);
    partitioned_hypergraph.setNodePart(id[4], 1);
    partitioned_hypergraph.setNodePart(id[5], 2);
    partitioned_hypergraph.setNodePart(id[6], 2);
    partitioned_hypergraph.initializeNumCutHyperedges(TBB::GLOBAL_TASK_GROUP);
  }

  static void SetUpTestSuite() {
    TBB::instance(HardwareTopology::instance().num_cpus());
  }

  void verifyPartitionPinCounts(const HyperedgeID he,
                                const std::vector<HypernodeID>& expected_pin_counts) {
    ASSERT(expected_pin_counts.size() == static_cast<size_t>(partitioned_hypergraph.k()));
    for (PartitionID block = 0; block < 3; ++block) {
      ASSERT_EQ(expected_pin_counts[block], partitioned_hypergraph.pinCountInPart(he, block)) << V(he) << V(block);
    }
  }

  void verifyConnectivitySet(const HyperedgeID he,
                             const std::set<PartitionID>& connectivity_set) {
    ASSERT_EQ(connectivity_set.size(), partitioned_hypergraph.connectivity(he)) << V(he);
    PartitionID connectivity = 0;
    for (const PartitionID& block : partitioned_hypergraph.connectivitySet(he)) {
      ASSERT_TRUE(connectivity_set.find(block) != connectivity_set.end()) << V(he) << V(block);
      ++connectivity;
    }
    ASSERT_EQ(connectivity_set.size(), connectivity) << V(he);
  }

  Hypergraph hypergraph;
  PartitionedHyperGraph partitioned_hypergraph;
  std::vector<HypernodeID> id;
};

template <class F1, class F2>
void executeConcurrent(const F1& f1, const F2& f2) {
  std::atomic<int> cnt(0);
  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < 2) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < 2) { }
    f2();
  });
}

typedef ::testing::Types<PartitionedHypergraphTypeTraits<
                          PartitionedHypergraph<StaticHypergraph, TBBNumaArena, true>,
                          StaticHypergraph,
                          StaticHypergraphFactory,
                          TBBNumaArena>,
                        PartitionedHypergraphTypeTraits<
                          PartitionedHypergraph<StaticHypergraph, TBBNumaArena, false>,
                          StaticHypergraph,
                          StaticHypergraphFactory,
                          TBBNumaArena>> PartitionedHypergraphTestTypes;

TYPED_TEST_CASE(APartitionedHypergraph, PartitionedHypergraphTestTypes);

TYPED_TEST(APartitionedHypergraph, HasCorrectPartWeightAndSizes) {
  ASSERT_EQ(3, this->partitioned_hypergraph.partWeight(0));
  ASSERT_EQ(3, this->partitioned_hypergraph.partSize(0));
  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(1));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(1));
  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(2));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(2));
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartWeightsIfOnlyOneThreadPerformsModifications) {
  ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));

  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(0));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(0));
  ASSERT_EQ(3, this->partitioned_hypergraph.partWeight(1));
  ASSERT_EQ(3, this->partitioned_hypergraph.partSize(1));
  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(2));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(2));
}


TYPED_TEST(APartitionedHypergraph, PerformsTwoConcurrentMovesWhereOnlyOneSucceeds) {
  std::array<bool, 2> success;
  executeConcurrent([&] {
    success[0] = this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1);
  }, [&] {
    success[1] = this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 2);
  });

  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(0));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(0));
  if ( success[0] ) {
    ASSERT_FALSE(success[1]);
    ASSERT_EQ(3, this->partitioned_hypergraph.partWeight(1));
    ASSERT_EQ(3, this->partitioned_hypergraph.partSize(1));
    ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(2));
    ASSERT_EQ(2, this->partitioned_hypergraph.partSize(2));
  } else {
    ASSERT_TRUE(success[1]);
    ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(1));
    ASSERT_EQ(2, this->partitioned_hypergraph.partSize(1));
    ASSERT_EQ(3, this->partitioned_hypergraph.partWeight(2));
    ASSERT_EQ(3, this->partitioned_hypergraph.partSize(2));
  }
}

TYPED_TEST(APartitionedHypergraph, PerformsConcurrentMovesWhereAllSucceed) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 2));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 2));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 1));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 2));
  });

  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(0));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(0));
  ASSERT_EQ(2, this->partitioned_hypergraph.partWeight(1));
  ASSERT_EQ(2, this->partitioned_hypergraph.partSize(1));
  ASSERT_EQ(3, this->partitioned_hypergraph.partWeight(2));
  ASSERT_EQ(3, this->partitioned_hypergraph.partSize(2));
}

TYPED_TEST(APartitionedHypergraph, HasCorrectInitialPartitionPinCounts) {
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 2, 0, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 2, 2, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 2, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent1) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[1], 0, 2));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 1, 1, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 3, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 2, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent2) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 2));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 0));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 2, 0, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 2, 1, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 1, 1, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 2, 0, 1 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent3) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 2));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 2));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 2, 0, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 2, 0, 2 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 0, 3 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent4) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 2));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 0));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 1, 0, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 2, 2, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 2, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent5) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 1));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 1, 1, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 1, 3, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 3, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 1, 1 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectPartitionPinCountsIfAllNodesMovesConcurrent) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 2));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 1));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[1], 0, 2));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 1));
  });

  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 0), { 0, 1, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 1), { 2, 1, 1 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 2), { 2, 1, 0 });
  this->verifyPartitionPinCounts(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 2, 1 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent1) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 0, 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent2) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 2));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 0, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 1, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent3) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 1));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 1, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 1, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent4) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 0 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent5) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[1], 0, 2));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 2));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 0 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 1, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 1, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 2 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectConnectivitySetIfAllNodesMovesConcurrent) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[0], 0, 1));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[1], 0, 2));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 1));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 1));
  });

  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 0), { 1 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 1), { 0, 1, 2 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 2), { 0 });
  this->verifyConnectivitySet(GLOBAL_EDGE_ID(this->hypergraph, 3), { 0, 1 });
}

TYPED_TEST(APartitionedHypergraph, HasCorrectInitialBorderNodes) {
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[0]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[1]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[2]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[3]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[4]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[5]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[6]));

  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[0]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[1]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[2]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[3]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[4]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[5]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[6]));
}

TYPED_TEST(APartitionedHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently1) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
  });

  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[0]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[1]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[2]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[3]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[4]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[5]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[6]));

  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[0]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[1]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[2]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[3]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[4]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[5]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[6]));
}

TYPED_TEST(APartitionedHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently2) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[1], 0, 1));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 1));
  });

  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[0]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[1]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[2]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[3]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[4]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[5]));
  ASSERT_TRUE(this->partitioned_hypergraph.isBorderNode(this->id[6]));

  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[0]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[1]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[2]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[3]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[4]));
  ASSERT_EQ(1, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[5]));
  ASSERT_EQ(2, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[6]));
}

TYPED_TEST(APartitionedHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently3) {
  executeConcurrent([&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[3], 1, 0));
  }, [&] {
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[5], 2, 0));
    ASSERT_TRUE(this->partitioned_hypergraph.changeNodePart(this->id[4], 1, 0));
  });

  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[0]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[1]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[2]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[3]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[4]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[5]));
  ASSERT_FALSE(this->partitioned_hypergraph.isBorderNode(this->id[6]));

  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[0]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[1]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[2]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[3]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[4]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[5]));
  ASSERT_EQ(0, this->partitioned_hypergraph.numIncidentCutHyperedges(this->id[6]));
}

}  // namespace ds
}  // namespace mt_kahypar