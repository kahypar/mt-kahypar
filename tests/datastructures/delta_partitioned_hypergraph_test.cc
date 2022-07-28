/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <atomic>
#include <mt-kahypar/parallel/tbb_initializer.h>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

class ADeltaPartitionedHypergraph : public Test {

 using DeltaPartitionedHyperGraph = DeltaPartitionedHypergraph<mt_kahypar::PartitionedHypergraph>;

 public:

  ADeltaPartitionedHypergraph() :
    hg(mt_kahypar::HypergraphFactory::construct(
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    phg(3, hg, parallel_tag_t()),
    delta_phg(),
    context() {
    phg.setOnlyNodePart(0, 0);
    phg.setOnlyNodePart(1, 0);
    phg.setOnlyNodePart(2, 0);
    phg.setOnlyNodePart(3, 1);
    phg.setOnlyNodePart(4, 1);
    phg.setOnlyNodePart(5, 2);
    phg.setOnlyNodePart(6, 2);
    phg.initializePartition();
    phg.initializeGainCache();

    context.partition.k = 3;
    delta_phg = std::move(DeltaPartitionedHyperGraph(context));
    delta_phg.setPartitionedHypergraph(&phg);
  }

  void verifyPinCounts(const HyperedgeID he,
                       const std::vector<HypernodeID>& expected_pin_counts) {
    ASSERT(expected_pin_counts.size() == static_cast<size_t>(phg.k()));
    for (PartitionID block = 0; block < 3; ++block) {
      ASSERT_EQ(expected_pin_counts[block], delta_phg.pinCountInPart(he, block)) << V(he) << V(block);
    }
  }

  void verifymoveToBenefit(const HypernodeID hn,
                           const std::vector<HypernodeID>& expected_penalties) {
    ASSERT(expected_penalties.size() == static_cast<size_t>(phg.k()));
    for (PartitionID block = 0; block < 3; ++block) {
      ASSERT_EQ(expected_penalties[block], delta_phg.moveToBenefit(hn, block)) << V(hn) << V(block);
    }
  }

  Hypergraph hg;
  mt_kahypar::PartitionedHypergraph phg;
  DeltaPartitionedHyperGraph delta_phg;
  Context context;
};

TEST_F(ADeltaPartitionedHypergraph, VerifiesInitialPinCounts) {
  verifyPinCounts(0, { 2, 0, 0 });
  verifyPinCounts(1, { 2, 2, 0 });
  verifyPinCounts(2, { 0, 2, 1 });
  verifyPinCounts(3, { 1, 0, 2 });
}

TEST_F(ADeltaPartitionedHypergraph, VerifyInitialmoveFromPenaltys) {
  ASSERT_EQ(2, delta_phg.moveFromPenalty(0));
  ASSERT_EQ(1, delta_phg.moveFromPenalty(1));
  ASSERT_EQ(1, delta_phg.moveFromPenalty(2));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(3));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(4));
  ASSERT_EQ(1, delta_phg.moveFromPenalty(5));
  ASSERT_EQ(1, delta_phg.moveFromPenalty(6));
}

TEST_F(ADeltaPartitionedHypergraph, VerifyInitialMoveToPenalties) {
  verifymoveToBenefit(0, { 2, 1, 0 });
  verifymoveToBenefit(1, { 1, 1, 0 });
  verifymoveToBenefit(2, { 2, 0, 1 });
  verifymoveToBenefit(3, { 1, 2, 1 });
  verifymoveToBenefit(4, { 1, 2, 1 });
  verifymoveToBenefit(5, { 1, 0, 1 });
  verifymoveToBenefit(6, { 1, 1, 2 });
}
TEST_F(ADeltaPartitionedHypergraph, MovesAVertex1) {
  delta_phg.changeNodePartWithGainCacheUpdate(1, 0, 1, 1000);
  ASSERT_EQ(0, phg.partID(1));
  ASSERT_EQ(1, delta_phg.partID(1));

  // Verify Pin Counts
  verifyPinCounts(1, { 1, 3, 0 });

  // Verify Move From Benefit
  ASSERT_EQ(1, delta_phg.moveFromPenalty(0));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(3));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(4));

  // Verify Move To Penalty
  verifymoveToBenefit(0, { 2, 1, 0 });
  verifymoveToBenefit(3, { 1, 2, 1 });
  verifymoveToBenefit(4, { 1, 2, 1 });
}

TEST_F(ADeltaPartitionedHypergraph, MovesAVertex2) {
  delta_phg.changeNodePartWithGainCacheUpdate(6, 2, 1, 1000);
  ASSERT_EQ(2, phg.partID(6));
  ASSERT_EQ(1, delta_phg.partID(6));

  // Verify Pin Counts
  verifyPinCounts(2, { 0, 3, 0 });
  verifyPinCounts(3, { 1, 1, 1 });

  // Verify Move From Benefit
  ASSERT_EQ(1, delta_phg.moveFromPenalty(2));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(3));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(4));
  ASSERT_EQ(0, delta_phg.moveFromPenalty(5));

  // Verify Move To Penalty
  verifymoveToBenefit(2, { 2, 1, 1 });
  verifymoveToBenefit(3, { 1, 2, 0 });
  verifymoveToBenefit(4, { 1, 2, 0 });
  verifymoveToBenefit(5, { 1, 1, 1 });
}

TEST_F(ADeltaPartitionedHypergraph, MovesSeveralVertices) {
  delta_phg.changeNodePartWithGainCacheUpdate(6, 2, 1, 1000);
  delta_phg.changeNodePartWithGainCacheUpdate(2, 0, 1, 1000);
  delta_phg.changeNodePartWithGainCacheUpdate(5, 2, 1, 1000);
  ASSERT_EQ(0, phg.partID(2));
  ASSERT_EQ(2, phg.partID(5));
  ASSERT_EQ(2, phg.partID(6));
  ASSERT_EQ(1, delta_phg.partID(2));
  ASSERT_EQ(1, delta_phg.partID(5));
  ASSERT_EQ(1, delta_phg.partID(6));

  // Verify Pin Counts
  verifyPinCounts(0, { 1, 1, 0 });
  verifyPinCounts(1, { 2, 2, 0 });
  verifyPinCounts(2, { 0, 3, 0 });
  verifyPinCounts(3, { 0, 3, 0 });

  // Verify Move From Benefit
  ASSERT_EQ(1, delta_phg.moveFromPenalty(0));
  ASSERT_EQ(1, delta_phg.moveFromPenalty(1));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(3));
  ASSERT_EQ(2, delta_phg.moveFromPenalty(4));

  // Verify Move To Penalty
  verifymoveToBenefit(0, { 2, 2, 0 });
  verifymoveToBenefit(1, { 1, 1, 0 });
  verifymoveToBenefit(2, { 1, 2, 0 });
  verifymoveToBenefit(3, { 1, 2, 0 });
  verifymoveToBenefit(4, { 1, 2, 0 });
  verifymoveToBenefit(5, { 0, 1, 0 });
  verifymoveToBenefit(6, { 0, 2, 0 });
}

} // namespace ds
} // namespace mt_kahypar