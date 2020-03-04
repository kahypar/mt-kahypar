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
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"

using ::testing::Test;

namespace mt_kahypar {

using TypeTraits = ds::TestTypeTraits<2>;
using HyperGraph = typename TypeTraits::HyperGraph;
using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
using TBB = typename TypeTraits::TBB;

template <template <typename> class GainPolicy, PartitionID K>
class AGainPolicy : public Test {
 public:
  using GainCalculator = GainPolicy<PartitionedHyperGraph>;

  AGainPolicy() :
    hg(HyperGraphFactory::construct(TBB::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    context(),
    gain(nullptr) {
    context.partition.k = K;
    context.partition.max_part_weights.assign(K, std::numeric_limits<HypernodeWeight>::max());
    gain = std::make_unique<GainCalculator>(context, true  /* disable randomization */);
    hypergraph = PartitionedHyperGraph(K, TBB::GLOBAL_TASK_GROUP, hg);
  }

  void assignPartitionIDs(const std::vector<PartitionID>& part_ids) {
    HypernodeID hn = 0;
    for (const PartitionID& part : part_ids) {
      ASSERT(part < K);
      hypergraph.setNodePart(hypergraph.globalNodeID(hn++), part);
    }
    hypergraph.initializeNumCutHyperedges();
  }

  HyperGraph hg;
  PartitionedHyperGraph hypergraph;
  Context context;
  std::unique_ptr<GainCalculator> gain;
};

using AKm1PolicyK2 = AGainPolicy<Km1Policy, 2>;

TEST_F(AKm1PolicyK2, ComputesCorrectMoveGainForVertex1) {
  assignPartitionIDs({ 1, 0, 0, 0, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(0));
  ASSERT_EQ(1, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(-2, move.gain);
}

TEST_F(AKm1PolicyK2, ComputesCorrectObjectiveDelta1) {
  assignPartitionIDs({ 1, 0, 0, 0, 0, 1, 1 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 1, 0,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-2, gain->delta());
}

TEST_F(AKm1PolicyK2, ComputesCorrectMoveGainForVertex2) {
  assignPartitionIDs({ 0, 0, 0, 1, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(3));
  ASSERT_EQ(1, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(AKm1PolicyK2, ComputesCorrectObjectiveDelta2) {
  assignPartitionIDs({ 0, 0, 0, 1, 0, 1, 1 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(AKm1PolicyK2, ComputesCorrectMoveGainForVertex3) {
  assignPartitionIDs({ 0, 0, 0, 0, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(4));
  ASSERT_EQ(0, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(0, move.gain);
}

using ACutPolicyK2 = AGainPolicy<CutPolicy, 2>;

TEST_F(ACutPolicyK2, ComputesCorrectMoveGainForVertex1) {
  assignPartitionIDs({ 1, 0, 0, 0, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(0));
  ASSERT_EQ(1, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(-2, move.gain);
}

TEST_F(ACutPolicyK2, ComputesCorrectObjectiveDelta1) {
  assignPartitionIDs({ 1, 0, 0, 0, 0, 1, 1 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 1, 0,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-2, gain->delta());
}

TEST_F(ACutPolicyK2, ComputesCorrectMoveGainForVertex2) {
  assignPartitionIDs({ 0, 0, 0, 1, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(3));
  ASSERT_EQ(1, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(ACutPolicyK2, ComputesCorrectObjectiveDelta2) {
  assignPartitionIDs({ 0, 0, 0, 1, 0, 1, 1 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(3), 1, 0,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(ACutPolicyK2, ComputesCorrectMoveGainForVertex3) {
  assignPartitionIDs({ 0, 0, 0, 0, 0, 1, 1 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(4));
  ASSERT_EQ(0, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(0, move.gain);
}

using AKm1PolicyK4 = AGainPolicy<Km1Policy, 4>;

TEST_F(AKm1PolicyK4, ComputesCorrectMoveGainForVertex1) {
  assignPartitionIDs({ 0, 1, 2, 3, 3, 1, 2 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(0));
  ASSERT_EQ(0, move.from);
  ASSERT_EQ(1, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(AKm1PolicyK4, ComputesCorrectObjectiveDelta1) {
  assignPartitionIDs({ 0, 1, 2, 3, 3, 1, 2 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 1,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(AKm1PolicyK4, ComputesCorrectMoveGainForVertex2) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(6));
  ASSERT_EQ(3, move.from);
  ASSERT_EQ(0, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(AKm1PolicyK4, ComputesCorrectObjectiveDelta2) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 3, 0,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(AKm1PolicyK4, ComputesCorrectMoveGainForVertex3) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(3));
  ASSERT_EQ(2, move.from);
  ASSERT_EQ(2, move.to);
  ASSERT_EQ(0, move.gain);
}

using ACutPolicyK4 = AGainPolicy<CutPolicy, 4>;

TEST_F(ACutPolicyK4, ComputesCorrectMoveGainForVertex1) {
  assignPartitionIDs({ 0, 1, 2, 3, 3, 1, 2 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(0));
  ASSERT_EQ(0, move.from);
  ASSERT_EQ(2, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(ACutPolicyK4, ComputesCorrectObjectiveDelta1) {
  assignPartitionIDs({ 0, 1, 2, 3, 3, 1, 2 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(0), 0, 2,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(ACutPolicyK4, ComputesCorrectMoveGainForVertex2) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(6));
  ASSERT_EQ(3, move.from);
  ASSERT_EQ(2, move.to);
  ASSERT_EQ(-1, move.gain);
}

TEST_F(ACutPolicyK4, ComputesCorrectObjectiveDelta2) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  ASSERT_TRUE(hypergraph.changeNodePart(hypergraph.globalNodeID(6), 3, 2,
                                        [&](const HyperedgeID he,
                                            const HyperedgeWeight edge_weight,
                                            const HypernodeID edge_size,
                                            const HypernodeID pin_count_in_from_part_after,
                                            const HypernodeID pin_count_in_to_part_after) {
      gain->computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    }));
  ASSERT_EQ(-1, gain->delta());
}

TEST_F(ACutPolicyK4, ComputesCorrectMoveGainForVertex3) {
  assignPartitionIDs({ 0, 3, 1, 2, 2, 0, 3 });
  Move move = gain->computeMaxGainMove(hypergraph, hypergraph.globalNodeID(3));
  ASSERT_EQ(2, move.from);
  ASSERT_EQ(2, move.to);
  ASSERT_EQ(0, move.gain);
}
}  // namespace mt_kahypar
