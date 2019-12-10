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
#include <atomic>

#include "gmock/gmock.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {
using AConcurrentHypergraph = AHypergraph<2>;
using TestHypergraph = typename AConcurrentHypergraph::TestHypergraph;
using TestStreamingHypergraph = typename AConcurrentHypergraph::TestStreamingHypergraph;
using Memento = typename TestStreamingHypergraph::Memento;

void assignPartitionIDs(TestHypergraph& hypergraph) {
  hypergraph.setNodePart(hypergraph.globalNodeID(0), 0);
  hypergraph.setNodePart(hypergraph.globalNodeID(1), 0);
  hypergraph.setNodePart(hypergraph.globalNodeID(2), 0);
  hypergraph.setNodePart(hypergraph.globalNodeID(3), 1);
  hypergraph.setNodePart(hypergraph.globalNodeID(4), 1);
  hypergraph.setNodePart(hypergraph.globalNodeID(5), 2);
  hypergraph.setNodePart(hypergraph.globalNodeID(6), 2);
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();
}

TestHypergraph construct_test_hypergraph_without_partition(const AConcurrentHypergraph& test) {
  TestHypergraph hypergraph = test.construct_hypergraph(7,
                                                        { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } },
                                                        { 0, 0, 0, 1, 1, 1, 1 },
                                                        { 0, 0, 1, 1 },
                                                        { 0, 0, 0, 1, 1, 2, 1 }, 3);
  return hypergraph;
}


TestHypergraph construct_test_hypergraph(const AConcurrentHypergraph& test) {
  TestHypergraph hypergraph = test.construct_hypergraph(7,
                                                        { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } },
                                                        { 0, 0, 0, 1, 1, 1, 1 },
                                                        { 0, 0, 1, 1 },
                                                        { 0, 0, 0, 1, 1, 2, 1 }, 3);
  assignPartitionIDs(hypergraph);
  return hypergraph;
}

template <class F, class K>
void executeConcurrent(F f1, K f2) {
  std::atomic<int> cnt(0);
  tbb::task_group group;

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f1();
      });

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f2();
      });

  group.wait();
}

void verifyPinIterators(const TestHypergraph& hypergraph,
                        const std::vector<HyperedgeID> hyperedges,
                        const std::vector<std::set<HypernodeID> >& references,
                        bool log = false) {
  ASSERT(hyperedges.size() == references.size());
  for (size_t i = 0; i < hyperedges.size(); ++i) {
    const HyperedgeID he = hyperedges[i];
    const std::set<HypernodeID>& reference = references[i];
    size_t count = 0;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      if (log) LOG << V(he) << V(pin);
      ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
      count++;
    }
    ASSERT_EQ(count, reference.size());
  }
}

TEST_F(AConcurrentHypergraph, HasCorrectLocalPartWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  executeConcurrent([&] {
        ASSERT_EQ(3, hypergraph.localPartWeight(0));
        ASSERT_EQ(3, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(2, hypergraph.localPartWeight(2));
        ASSERT_EQ(2, hypergraph.localPartSize(2));
      }, [&] {
        ASSERT_EQ(3, hypergraph.localPartWeight(0));
        ASSERT_EQ(3, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(2, hypergraph.localPartWeight(2));
        ASSERT_EQ(2, hypergraph.localPartSize(2));
      });
}

TEST_F(AConcurrentHypergraph, HasCorrectLocalPartWeightsIfOnlyOneThreadPerformsModificationsBefore) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  hypergraph.changeNodePart(0, 0, 1);
  hypergraph.updateGlobalPartInfos();

  executeConcurrent([&] {
        ASSERT_EQ(2, hypergraph.localPartWeight(0));
        ASSERT_EQ(2, hypergraph.localPartSize(0));
        ASSERT_EQ(3, hypergraph.localPartWeight(1));
        ASSERT_EQ(3, hypergraph.localPartSize(1));
        ASSERT_EQ(2, hypergraph.localPartWeight(2));
        ASSERT_EQ(2, hypergraph.localPartSize(2));
      }, [&] {
        ASSERT_EQ(2, hypergraph.localPartWeight(0));
        ASSERT_EQ(2, hypergraph.localPartSize(0));
        ASSERT_EQ(3, hypergraph.localPartWeight(1));
        ASSERT_EQ(3, hypergraph.localPartSize(1));
        ASSERT_EQ(2, hypergraph.localPartWeight(2));
        ASSERT_EQ(2, hypergraph.localPartSize(2));
      });
}

TEST_F(AConcurrentHypergraph, PerformsTwoConcurrentMovesWhereOnlyOneSucceeds) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  std::array<bool, 2> success;
  executeConcurrent([&] {
        success[0] = hypergraph.changeNodePart(0, 0, 1);

        if (success[0]) {
          ASSERT_EQ(2, hypergraph.localPartWeight(0));
          ASSERT_EQ(2, hypergraph.localPartSize(0));
          ASSERT_EQ(3, hypergraph.localPartWeight(1));
          ASSERT_EQ(3, hypergraph.localPartSize(1));
        }
      }, [&] {
        success[1] = hypergraph.changeNodePart(0, 0, 2);

        if (success[1]) {
          ASSERT_EQ(2, hypergraph.localPartWeight(0));
          ASSERT_EQ(2, hypergraph.localPartSize(0));
          ASSERT_EQ(3, hypergraph.localPartWeight(2));
          ASSERT_EQ(3, hypergraph.localPartSize(2));
        }
      });

  ASSERT_TRUE((success[0] && !success[1]) || (!success[0] && success[1]));

  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(2, hypergraph.partSize(0));
  if (success[0]) {
    ASSERT_EQ(3, hypergraph.partWeight(1));
    ASSERT_EQ(3, hypergraph.partSize(1));
    ASSERT_EQ(2, hypergraph.partWeight(2));
    ASSERT_EQ(2, hypergraph.partSize(2));
  } else {
    ASSERT_EQ(2, hypergraph.partWeight(1));
    ASSERT_EQ(2, hypergraph.partSize(1));
    ASSERT_EQ(3, hypergraph.partWeight(2));
    ASSERT_EQ(3, hypergraph.partSize(2));
  }
}

TEST_F(AConcurrentHypergraph, PerformsConcurrentMovesWhereAllSucceed) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 2));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 2));

        ASSERT_EQ(1, hypergraph.localPartWeight(0));
        ASSERT_EQ(1, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(4, hypergraph.localPartWeight(2));
        ASSERT_EQ(4, hypergraph.localPartSize(2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 1));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 2));

        ASSERT_EQ(4, hypergraph.localPartWeight(0));
        ASSERT_EQ(4, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(1, hypergraph.localPartWeight(2));
        ASSERT_EQ(1, hypergraph.localPartSize(2));
      });

  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(2, hypergraph.partSize(0));
  ASSERT_EQ(2, hypergraph.partWeight(1));
  ASSERT_EQ(2, hypergraph.partSize(1));
  ASSERT_EQ(3, hypergraph.partWeight(2));
  ASSERT_EQ(3, hypergraph.partSize(2));
}

TEST_F(AConcurrentHypergraph, PerformsConcurrentMovesAndUpdatesLocalPartInfos) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  std::atomic<size_t> cnt(0);
  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));  // Move 1

        cnt++;
        while (cnt < 3) { }

        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 2));  // Move 5

        hypergraph.updateLocalPartInfos();  // Move 1, 2, 3, 4, 5 are applied
        ASSERT_EQ(3, hypergraph.localPartWeight(0));
        ASSERT_EQ(3, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(2, hypergraph.localPartWeight(2));
        ASSERT_EQ(2, hypergraph.localPartSize(2));

        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 2));  // Move 6

        ASSERT_EQ(2, hypergraph.localPartWeight(0));
        ASSERT_EQ(2, hypergraph.localPartSize(0));
        ASSERT_EQ(2, hypergraph.localPartWeight(1));
        ASSERT_EQ(2, hypergraph.localPartSize(1));
        ASSERT_EQ(3, hypergraph.localPartWeight(2));
        ASSERT_EQ(3, hypergraph.localPartSize(2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 1));  // Move 2
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));  // Move 3

        cnt++;
        while (cnt < 2) { }

        hypergraph.updateLocalPartInfos();  // Move 1, 2, 3 are applied
        ASSERT_EQ(3, hypergraph.localPartWeight(0));
        ASSERT_EQ(3, hypergraph.localPartSize(0));
        ASSERT_EQ(4, hypergraph.localPartWeight(1));
        ASSERT_EQ(4, hypergraph.localPartSize(1));
        ASSERT_EQ(0, hypergraph.localPartWeight(2));
        ASSERT_EQ(0, hypergraph.localPartSize(2));

        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 2));  // Move 4

        cnt++;

        ASSERT_EQ(3, hypergraph.localPartWeight(0));
        ASSERT_EQ(3, hypergraph.localPartSize(0));
        ASSERT_EQ(3, hypergraph.localPartWeight(1));
        ASSERT_EQ(3, hypergraph.localPartSize(1));
        ASSERT_EQ(1, hypergraph.localPartWeight(2));
        ASSERT_EQ(1, hypergraph.localPartSize(2));
      });

  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(2, hypergraph.partSize(0));
  ASSERT_EQ(2, hypergraph.partWeight(1));
  ASSERT_EQ(2, hypergraph.partSize(1));
  ASSERT_EQ(3, hypergraph.partWeight(2));
  ASSERT_EQ(3, hypergraph.partSize(2));
}

void verifyPartitionPinCounts(TestHypergraph& hypergraph,
                              const HyperedgeID he,
                              const std::vector<HypernodeID>& expected_pin_counts) {
  ASSERT(expected_pin_counts.size() == 3);
  for (PartitionID k = 0; k < 3; ++k) {
    ASSERT_EQ(expected_pin_counts[k], hypergraph.pinCountInPart(he, k)) << V(he) << V(k);
  }
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCounts) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 2, 0, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 2, 2, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 0, 2, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 1, 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(1), 0, 2));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 1, 1, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 0, 3, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 0, 2, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 1, 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 2, 0, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 2, 1, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 1, 1, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 2, 0, 1 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 2));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 2, 0, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 2, 0, 2 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 0, 0, 3 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 1, 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 0));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 1, 0, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 2, 2, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 0, 2, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 1, 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfTwoNodesMovesConcurrent5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 1));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 1, 1, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 1, 3, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 0, 3, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 1, 1, 1 });
}

TEST_F(AConcurrentHypergraph, HasCorrectPartitionPinCountsIfAllNodesMovesConcurrent) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 2));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 1));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(1), 0, 2));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 1));
      });

  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(0), { 0, 1, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(1), { 2, 1, 1 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(2), { 2, 1, 0 });
  verifyPartitionPinCounts(hypergraph, hypergraph.globalEdgeID(3), { 0, 2, 1 });
}

void verifyConnectivitySet(TestHypergraph& hypergraph,
                           const HyperedgeID he,
                           const std::set<PartitionID>& connectivity_set) {
  ASSERT_EQ(connectivity_set.size(), hypergraph.connectivity(he)) << V(he);
  PartitionID connectivity = 0;
  for (const PartitionID& id : hypergraph.connectivitySet(he)) {
    ASSERT_TRUE(connectivity_set.find(id) != connectivity_set.end()) << V(he) << V(id);
    ++connectivity;
  }
  ASSERT_EQ(connectivity_set.size(), connectivity) << V(he);
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 0, 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0, 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 0, 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 2));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 0, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0, 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 1, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 1));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0, 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 1, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 1, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 0 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 0, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfTwoNodesMovesConcurrent5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(1), 0, 2));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 2));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 0 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0, 1, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 1, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 0, 2 });
}

TEST_F(AConcurrentHypergraph, HasCorrectConnectivitySetIfAllNodesMovesConcurrent) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(1), 0, 2));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 1));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 1));
      });

  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(0), { 1 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(1), { 0, 1, 2 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(2), { 0 });
  verifyConnectivitySet(hypergraph, hypergraph.globalEdgeID(3), { 0, 1 });
}

TEST_F(AConcurrentHypergraph, HasCorrectBorderNodes) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(0)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(1)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(2)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(3)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(4)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(5)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(6)));

  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(0)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(1)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(2)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(3)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(6)));
}

TEST_F(AConcurrentHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
      });

  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(0)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(1)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(2)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(3)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(4)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(5)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(6)));

  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(1)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(2)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(3)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(6)));
}

TEST_F(AConcurrentHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(1), 0, 1));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(2), 0, 1));
      });

  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(0)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(1)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(2)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(3)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(4)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(5)));
  ASSERT_TRUE(hypergraph.isBorderNode(hypergraph.globalNodeID(6)));

  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(0)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(1)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(2)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(3)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(4)));
  ASSERT_EQ(1, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(5)));
  ASSERT_EQ(2, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(6)));
}

TEST_F(AConcurrentHypergraph, HasCorrectBorderNodesIfNodesAreMovingConcurrently3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  executeConcurrent([&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 2, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0));
      }, [&] {
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(5), 2, 0));
        ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(4), 1, 0));
      });

  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(0)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(1)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(2)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(3)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(4)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(5)));
  ASSERT_FALSE(hypergraph.isBorderNode(hypergraph.globalNodeID(6)));

  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(0)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(1)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(2)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(3)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(4)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(5)));
  ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hypergraph.globalNodeID(6)));
}

TEST_F(AConcurrentHypergraph, UncontractsABatchOfContractionsConcurrently1)  {
  TestHypergraph hypergraph = construct_test_hypergraph_without_partition(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  kahypar::ds::FastResetFlagArray<> batch_hypernodes(hypergraph.initialNumNodes());

  std::vector<Memento> batch;
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(0), hypergraph.globalNodeID(2)));
  batch_hypernodes.set(2, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(3), hypergraph.globalNodeID(4)));
  batch_hypernodes.set(4, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(5), hypergraph.globalNodeID(6)));
  batch_hypernodes.set(6, true);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0] }, { id[0], id[1], id[3] }, { id[3], id[5] }, { id[0], id[5] } });

  hypergraph.buildContractionHierarchy(batch);
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[3], 0);
  hypergraph.setNodePart(id[5], 1);
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();

  std::reverse(batch.begin(), batch.end());
  for ( const Memento& memento : batch ) {
    hypergraph.preprocessMementoForBatchUncontraction(
      memento, parallel_he_representative, batch_hypernodes);
  }
  hypergraph.uncontract(batch, parallel_he_representative, batch_hypernodes);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AConcurrentHypergraph, UncontractsABatchOfContractionsConcurrently2)  {
  TestHypergraph hypergraph = construct_test_hypergraph_without_partition(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  kahypar::ds::FastResetFlagArray<> batch_hypernodes(hypergraph.initialNumNodes());

  std::vector<Memento> batch;
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(0), hypergraph.globalNodeID(3)));
  batch_hypernodes.set(3, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(1), hypergraph.globalNodeID(4)));
  batch_hypernodes.set(4, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(2), hypergraph.globalNodeID(5)));
  batch_hypernodes.set(5, true);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0], id[2] }, { id[0], id[1] }, { id[0], id[1], id[6] }, { id[2], id[6] } });

  hypergraph.buildContractionHierarchy(batch);
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 0);
  hypergraph.setNodePart(id[6], 1);
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();

  std::reverse(batch.begin(), batch.end());
  for ( const Memento& memento : batch ) {
    hypergraph.preprocessMementoForBatchUncontraction(
      memento, parallel_he_representative, batch_hypernodes);
  }
  hypergraph.uncontract(batch, parallel_he_representative, batch_hypernodes);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AConcurrentHypergraph, UncontractsABatchOfContractionsConcurrently3)  {
  TestHypergraph hypergraph = construct_test_hypergraph_without_partition(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  kahypar::ds::FastResetFlagArray<> batch_hypernodes(hypergraph.initialNumNodes());

  std::vector<Memento> batch;
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(0), hypergraph.globalNodeID(5)));
  batch_hypernodes.set(5, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(3), hypergraph.globalNodeID(2)));
  batch_hypernodes.set(2, true);
  batch.emplace_back(hypergraph.contract(hypergraph.globalNodeID(1), hypergraph.globalNodeID(6)));
  batch_hypernodes.set(6, true);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0], id[3] }, { id[0], id[1], id[3], id[4] }, { id[1], id[3], id[4] }, { id[0], id[1], id[3] } });

  hypergraph.buildContractionHierarchy(batch);
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[3], 0);
  hypergraph.setNodePart(id[4], 1);
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();

  std::reverse(batch.begin(), batch.end());
  for ( const Memento& memento : batch ) {
    hypergraph.preprocessMementoForBatchUncontraction(
      memento, parallel_he_representative, batch_hypernodes);
  }
  hypergraph.uncontract(batch, parallel_he_representative, batch_hypernodes);

  verifyPinIterators(hypergraph, { hypergraph.globalEdgeID(0), hypergraph.globalEdgeID(1),
    hypergraph.globalEdgeID(2), hypergraph.globalEdgeID(3) },
    { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}
}  // namespace ds
}  // namespace mt_kahypar
