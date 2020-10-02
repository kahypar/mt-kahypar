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

#include "tests/datastructures/hypergraph_fixtures.h"

#include <boost/range/irange.hpp>
#include <mt-kahypar/partition/refinement/policies/gain_policy.h>
#include "gmock/gmock.h"

#include "tbb/blocked_range.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {
template <PartitionID k,
          kahypar::Objective objective>
struct TestConfig {
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = objective;
};

template <typename Config>
class AConcurrentHypergraph : public Test {

 public:
  AConcurrentHypergraph() :
    k(Config::K),
    objective(Config::OBJECTIVE),
    underlying_hypergraph(),
    hypergraph()
  {
    int cpu_id = sched_getcpu();
    underlying_hypergraph = io::readHypergraphFile(
      "../tests/instances/contracted_ibm01.hgr", TBBNumaArena::GLOBAL_TASK_GROUP);
    hypergraph = mt_kahypar::PartitionedHypergraph(k, TBBNumaArena::GLOBAL_TASK_GROUP, underlying_hypergraph);
    for (const HypernodeID& hn : hypergraph.nodes()) {
      PartitionID id = utils::Randomize::instance().getRandomInt(0, k - 1, cpu_id);
      hypergraph.setNodePart(hn, id);
    }

    hypergraph.initializeGainCache();
  }

  static void SetUpTestSuite() {
    utils::Randomize::instance().setSeed(0);
  }

  PartitionID k;
  kahypar::Objective objective;
  Hypergraph underlying_hypergraph;
  mt_kahypar::PartitionedHypergraph hypergraph;
};

typedef ::testing::Types<TestConfig<2, kahypar::Objective::cut>,
                         TestConfig<4, kahypar::Objective::cut>,
                         TestConfig<8, kahypar::Objective::cut>,
                         TestConfig<16, kahypar::Objective::cut>,
                         TestConfig<32, kahypar::Objective::cut>,
                         TestConfig<64, kahypar::Objective::cut>,
                         TestConfig<128, kahypar::Objective::cut>,
                         TestConfig<2, kahypar::Objective::km1>,
                         TestConfig<4, kahypar::Objective::km1>,
                         TestConfig<8, kahypar::Objective::km1>,
                         TestConfig<16, kahypar::Objective::km1>,
                         TestConfig<32, kahypar::Objective::km1>,
                         TestConfig<64, kahypar::Objective::km1>,
                         TestConfig<128, kahypar::Objective::km1>> TestConfigs;

TYPED_TEST_CASE(AConcurrentHypergraph, TestConfigs);

template<typename HyperGraph>
void moveAllNodesOfHypergraphRandom(HyperGraph& hypergraph,
                                    const PartitionID k,
                                    const kahypar::Objective objective,
                                    const bool show_timings) {

  tbb::enumerable_thread_specific<HyperedgeWeight> deltas(0);

  auto objective_delta = [&](const HyperedgeID he,
                             const HyperedgeWeight edge_weight,
                             const HypernodeID edge_size,
                             const HypernodeID pin_count_in_from_part_after,
                             const HypernodeID pin_count_in_to_part_after) {
                           if (objective == kahypar::Objective::km1) {
                             deltas.local() += km1Delta(
                               he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
                           } else if (objective == kahypar::Objective::cut) {
                             deltas.local() += cutDelta(
                               he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
                           }
                         };

  HyperedgeWeight metric_before = metrics::objective(hypergraph, objective);
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  tbb::parallel_for(ID(0), hypergraph.initialNumNodes(), [&](const HypernodeID& hn) {
    int cpu_id = sched_getcpu();
    const PartitionID from = hypergraph.partID(hn);
    PartitionID to = -1;
    while (to == -1 || to == from) {
      to = utils::Randomize::instance().getRandomInt(0, k - 1, cpu_id);
    }
    ASSERT((to >= 0 && to < k) && to != from);
    hypergraph.changeNodePart(hn, from, to, objective_delta);
  } );

  hypergraph.recomputePartWeights();

  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  double timing = std::chrono::duration<double>(end - start).count();

  HyperedgeWeight delta = 0;
  for (const HyperedgeWeight& local_delta : deltas) {
    delta += local_delta;
  }

  HyperedgeWeight metric_after = metrics::objective(hypergraph, objective);
  ASSERT_EQ(metric_after, metric_before + delta) << V(metric_before) << V(delta);
  if (show_timings) {
    LOG << V(k) << V(objective) << V(metric_before) << V(delta) << V(metric_after) << V(timing);
  }
}

template<typename HyperGraph>
void verifyBlockWeightsAndSizes(HyperGraph& hypergraph,
                                const PartitionID k) {
  std::vector<HypernodeWeight> block_weight(k, 0);
  for (const HypernodeID& hn : hypergraph.nodes()) {
    block_weight[hypergraph.partID(hn)] += hypergraph.nodeWeight(hn);
  }

  for (PartitionID i = 0; i < k; ++i) {
    ASSERT_EQ(block_weight[i], hypergraph.partWeight(i));
  }
}

template<typename HyperGraph>
void verifyPinCountsInParts(HyperGraph& hypergraph,
                            const PartitionID k) {
  for (const HyperedgeID& he : hypergraph.edges()) {
    std::vector<HypernodeID> pin_count_in_part(k, 0);
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      ++pin_count_in_part[hypergraph.partID(pin)];
    }

    for (PartitionID i = 0; i < k; ++i) {
      ASSERT_EQ(pin_count_in_part[i], hypergraph.pinCountInPart(he, i));
    }
  }
}

template<typename HyperGraph>
void verifyConnectivitySet(HyperGraph& hypergraph,
                           const PartitionID k) {
  for (const HyperedgeID& he : hypergraph.edges()) {
    std::vector<HypernodeID> pin_count_in_part(k, 0);
    std::set<PartitionID> recomputed_connectivity_set;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      PartitionID id = hypergraph.partID(pin);
      HypernodeID pin_count_before = pin_count_in_part[id]++;
      if (pin_count_before == 0) {
        recomputed_connectivity_set.insert(id);
      }
    }

    std::set<PartitionID> connectivity_set;
    for (const PartitionID id : hypergraph.connectivitySet(he)) {
      ASSERT(id < k && id >= 0);
      connectivity_set.insert(id);
    }

    ASSERT_EQ(hypergraph.connectivity(he), connectivity_set.size());
    ASSERT_EQ(connectivity_set.size(), recomputed_connectivity_set.size());
    ASSERT_EQ(connectivity_set, recomputed_connectivity_set);
  }
}

template<typename HyperGraph>
void verifyBorderNodes(HyperGraph& hypergraph) {
  for (const HypernodeID& hn : hypergraph.nodes()) {
    bool is_border_node = false;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if (hypergraph.connectivity(he) > 1) {
        is_border_node = true;
        break;
      }
    }
    ASSERT_EQ(is_border_node, hypergraph.isBorderNode(hn));
  }
}

TYPED_TEST(AConcurrentHypergraph, VerifyBlockWeightsSmokeTest) {
  moveAllNodesOfHypergraphRandom(this->hypergraph, this->k, this->objective, false);
  verifyBlockWeightsAndSizes(this->hypergraph, this->k);
}

TYPED_TEST(AConcurrentHypergraph, VerifyPinCountsInPartsSmokeTest) {
  moveAllNodesOfHypergraphRandom(this->hypergraph, this->k, this->objective, false);
  verifyPinCountsInParts(this->hypergraph, this->k);
}

TYPED_TEST(AConcurrentHypergraph, VerifyConnectivitySetSmokeTest) {
  moveAllNodesOfHypergraphRandom(this->hypergraph, this->k, this->objective, false);
  verifyConnectivitySet(this->hypergraph, this->k);
}

TYPED_TEST(AConcurrentHypergraph, VerifyBorderNodesSmokeTest) {
  moveAllNodesOfHypergraphRandom(this->hypergraph, this->k, this->objective, false);
  verifyBorderNodes(this->hypergraph);
}


}  // namespace ds
}  // namespace mt_kahypar
