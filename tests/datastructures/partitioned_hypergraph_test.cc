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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/datastructures/numa_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/numa_partitioned_hypergraph.h"

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
 using Factory = typename TypeTraits::Factory;

 public:
 using Hypergraph = typename TypeTraits::Hypergraph;
 using TBB = typename TypeTraits::TBB;

  APartitionedHypergraph() :
    hypergraph(Factory::construct(TBB::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    partitioned_hypergraph(3, TBB::GLOBAL_TASK_GROUP, hypergraph),
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

  void verifyPins(const Hypergraph& hg,
                  const std::vector<HyperedgeID> hyperedges,
                  const std::vector< std::set<HypernodeID> >& references,
                  bool log = false) {
    ASSERT(hyperedges.size() == references.size());
    for (size_t i = 0; i < hyperedges.size(); ++i) {
      const HyperedgeID he = hyperedges[i];
      const std::set<HypernodeID>& reference = references[i];
      size_t count = 0;
      for (const HypernodeID& pin : hg.pins(he)) {
        if (log) LOG << V(he) << V(pin);
        ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
        count++;
      }
      ASSERT_EQ(count, reference.size());
    }
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

// Mocking Numa Architecture (=> 2 NUMA Nodes)
using TypeTraits = TestTypeTraits<2>;
using HwTopology = typename TypeTraits::HwTopology;
using TBB = typename TypeTraits::TBB;

// Define NUMA Hypergraph and Factory
using NumaHyperGraph = NumaHypergraph<StaticHypergraph, HwTopology, TBB>;
using NumaHyperGraphFactory = NumaHypergraphFactory<
  StaticHypergraph, StaticHypergraphFactory, HwTopology, TBB>;

typedef ::testing::Types<PartitionedHypergraphTypeTraits<
                          PartitionedHypergraph<StaticHypergraph, StaticHypergraphFactory, true>,
                          StaticHypergraph,
                          StaticHypergraphFactory,
                          TBBNumaArena>,
                        PartitionedHypergraphTypeTraits<
                          PartitionedHypergraph<StaticHypergraph, StaticHypergraphFactory, false>,
                          StaticHypergraph,
                          StaticHypergraphFactory,
                          TBBNumaArena>,
                        PartitionedHypergraphTypeTraits<
                          NumaPartitionedHypergraph<NumaHyperGraph, NumaHyperGraphFactory, true>,
                          NumaHyperGraph,
                          NumaHyperGraphFactory,
                          TBB>,
                        PartitionedHypergraphTypeTraits<
                          NumaPartitionedHypergraph<NumaHyperGraph, NumaHyperGraphFactory, false>,
                          NumaHyperGraph,
                          NumaHyperGraphFactory,
                          TBB>> PartitionedHypergraphTestTypes;

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

TYPED_TEST(APartitionedHypergraph, ExtractBlockZeroWithCutNetSplitting) {
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 0, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(3, hg.initialNumNodes());
  ASSERT_EQ(2, hg.initialNumEdges());
  ASSERT_EQ(4, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[0]),
    map_from_original_to_extracted_hg(this->id[1]),
    map_from_original_to_extracted_hg(this->id[2])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0), hg.globalEdgeID(1) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[0]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[1]), common::get_numa_node_of_vertex(node_id[1]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[2]), common::get_numa_node_of_vertex(node_id[2]));

  this->verifyPins(hg, {edge_id[0], edge_id[1]},
    { {node_id[0], node_id[2]}, {node_id[0], node_id[1]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockOneWithCutNetSplitting) {
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 1, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(2, hg.initialNumNodes());
  ASSERT_EQ(2, hg.initialNumEdges());
  ASSERT_EQ(4, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[3]),
    map_from_original_to_extracted_hg(this->id[4])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0), hg.globalEdgeID(1) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[3]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[4]), common::get_numa_node_of_vertex(node_id[1]));

  this->verifyPins(hg, {edge_id[0], edge_id[1]},
    { {node_id[0], node_id[1]}, {node_id[0], node_id[1]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockTwoWithCutNetSplitting) {
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 2, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(2, hg.initialNumNodes());
  ASSERT_EQ(1, hg.initialNumEdges());
  ASSERT_EQ(2, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[5]),
    map_from_original_to_extracted_hg(this->id[6])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[5]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[6]), common::get_numa_node_of_vertex(node_id[1]));

  this->verifyPins(hg, {edge_id[0]},
    { {node_id[0], node_id[1]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockZeroWithCutNetRemoval) {
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 0, false);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(3, hg.initialNumNodes());
  ASSERT_EQ(1, hg.initialNumEdges());
  ASSERT_EQ(2, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[0]),
    map_from_original_to_extracted_hg(this->id[1]),
    map_from_original_to_extracted_hg(this->id[2])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[0]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[1]), common::get_numa_node_of_vertex(node_id[1]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[2]), common::get_numa_node_of_vertex(node_id[2]));

  this->verifyPins(hg, {edge_id[0]},
    { {node_id[0], node_id[2]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockOneWithCutNetRemoval) {
  this->partitioned_hypergraph.changeNodePart(this->id[6], 2, 1);
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 1, false);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(3, hg.initialNumNodes());
  ASSERT_EQ(1, hg.initialNumEdges());
  ASSERT_EQ(3, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[3]),
    map_from_original_to_extracted_hg(this->id[4]),
    map_from_original_to_extracted_hg(this->id[6])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[3]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[4]), common::get_numa_node_of_vertex(node_id[1]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[6]), common::get_numa_node_of_vertex(node_id[2]));

  this->verifyPins(hg, {edge_id[0]},
    { {node_id[0], node_id[1], node_id[2]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockTwoWithCutNetRemoval) {
  this->partitioned_hypergraph.changeNodePart(this->id[2], 0, 2);
  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 2, false);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  ASSERT_EQ(3, hg.initialNumNodes());
  ASSERT_EQ(1, hg.initialNumEdges());
  ASSERT_EQ(3, hg.initialNumPins());

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };
  parallel::scalable_vector<HypernodeID> node_id = {
    map_from_original_to_extracted_hg(this->id[2]),
    map_from_original_to_extracted_hg(this->id[5]),
    map_from_original_to_extracted_hg(this->id[6])
  };
  parallel::scalable_vector<HypernodeID> edge_id = {
    hg.globalEdgeID(0) };

  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[2]), common::get_numa_node_of_vertex(node_id[0]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[5]), common::get_numa_node_of_vertex(node_id[1]));
  ASSERT_EQ(common::get_numa_node_of_vertex(this->id[6]), common::get_numa_node_of_vertex(node_id[2]));

  this->verifyPins(hg, {edge_id[0]},
    { {node_id[0], node_id[1], node_id[2]} });
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockZeroWithCommunityInformation) {
  this->hypergraph.setCommunityID(this->id[0], 0);
  this->hypergraph.setCommunityID(this->id[1], 1);
  this->hypergraph.setCommunityID(this->id[2], 0);
  this->hypergraph.setCommunityID(this->id[3], 2);
  this->hypergraph.setCommunityID(this->id[4], 3);
  this->hypergraph.setCommunityID(this->id[5], 4);
  this->hypergraph.setCommunityID(this->id[6], 5);
  this->hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 0, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };

  ASSERT_EQ(0, hg.communityID(map_from_original_to_extracted_hg(this->id[0])));
  ASSERT_EQ(1, hg.communityID(map_from_original_to_extracted_hg(this->id[1])));
  ASSERT_EQ(0, hg.communityID(map_from_original_to_extracted_hg(this->id[2])));
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockOneWithCommunityInformation) {
  this->hypergraph.setCommunityID(this->id[0], 0);
  this->hypergraph.setCommunityID(this->id[1], 1);
  this->hypergraph.setCommunityID(this->id[2], 0);
  this->hypergraph.setCommunityID(this->id[3], 2);
  this->hypergraph.setCommunityID(this->id[4], 3);
  this->hypergraph.setCommunityID(this->id[5], 4);
  this->hypergraph.setCommunityID(this->id[6], 5);
  this->hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 1, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };

  ASSERT_EQ(2, hg.communityID(map_from_original_to_extracted_hg(this->id[3])));
  ASSERT_EQ(3, hg.communityID(map_from_original_to_extracted_hg(this->id[4])));
}

TYPED_TEST(APartitionedHypergraph, ExtractBlockTwoWithCommunityInformation) {
  this->hypergraph.setCommunityID(this->id[0], 0);
  this->hypergraph.setCommunityID(this->id[1], 1);
  this->hypergraph.setCommunityID(this->id[2], 0);
  this->hypergraph.setCommunityID(this->id[3], 2);
  this->hypergraph.setCommunityID(this->id[4], 3);
  this->hypergraph.setCommunityID(this->id[5], 4);
  this->hypergraph.setCommunityID(this->id[6], 5);
  this->hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);

  auto extracted_hg = this->partitioned_hypergraph.extract(TBB::GLOBAL_TASK_GROUP, 2, true);
  auto& hg = extracted_hg.first;
  auto& hn_mapping = extracted_hg.second;

  auto map_from_original_to_extracted_hg = [&](const HypernodeID hn) {
    return hg.globalNodeID(hn_mapping[this->hypergraph.originalNodeID(hn)]);
  };

  ASSERT_EQ(4, hg.communityID(map_from_original_to_extracted_hg(this->id[5])));
  ASSERT_EQ(5, hg.communityID(map_from_original_to_extracted_hg(this->id[6])));
}

}  // namespace ds
}  // namespace mt_kahypar