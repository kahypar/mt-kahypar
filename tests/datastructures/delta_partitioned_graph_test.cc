/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <atomic>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/partitioned_graph.h"
#include "mt-kahypar/datastructures/delta_partitioned_graph.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

class ADeltaPartitionedGraph : public Test {

 using DeltaPartitionedGraph = ds::DeltaPartitionedGraph<mt_kahypar::PartitionedHypergraph>;

 public:

  ADeltaPartitionedGraph() :
    hg(mt_kahypar::HypergraphFactory::construct(7 , 6,
      { {1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6} }, nullptr, nullptr, true)),
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
    delta_phg = std::move(DeltaPartitionedGraph(context));
    delta_phg.setPartitionedHypergraph(&phg);
  }

  void verifyPinCounts(const HyperedgeID he,
                       const std::vector<HypernodeID>& expected_pin_counts) {
    ASSERT(expected_pin_counts.size() == static_cast<size_t>(phg.k()));
    for (PartitionID block = 0; block < 3; ++block) {
      ASSERT_EQ(expected_pin_counts[block], delta_phg.pinCountInPart(he, block)) << V(he) << V(block);
    }
  }

  void verifyKm1Gain(const HypernodeID hn,
                           const std::vector<HyperedgeWeight>& expected_penalties) {
    ASSERT(expected_penalties.size() == static_cast<size_t>(phg.k()));
    for (PartitionID block = 0; block < 3; ++block) {
      if (block != delta_phg.partID(hn)) {
        ASSERT_EQ(expected_penalties[block], delta_phg.km1Gain(hn, delta_phg.partID(hn), block)) << V(hn) << "; " << V(block);
      } else {
        ASSERT_EQ(delta_phg.moveToBenefit(hn, block), delta_phg.moveFromPenalty(hn)) << V(hn) << "; " << V(block);
      }
    }
  }

  Hypergraph hg;
  mt_kahypar::PartitionedHypergraph phg;
  DeltaPartitionedGraph delta_phg;
  Context context;
};

TEST_F(ADeltaPartitionedGraph, VerifiesInitialPinCounts) {
  // edge 1 - 2
  verifyPinCounts(0, { 2, 0, 0 });
  verifyPinCounts(2, { 2, 0, 0 });
  // edge 1 - 4
  verifyPinCounts(1, { 1, 1, 0 });
  verifyPinCounts(5, { 1, 1, 0 });
  // edge 2 - 3
  verifyPinCounts(3, { 1, 1, 0 });
  verifyPinCounts(4, { 1, 1, 0 });
  // edge 4 - 5
  verifyPinCounts(6, { 0, 1, 1 });
  verifyPinCounts(8, { 0, 1, 1 });
  // edge 4 - 6
  verifyPinCounts(7, { 0, 1, 1 });
  verifyPinCounts(10, { 0, 1, 1 });
  // edge 5 - 6
  verifyPinCounts(9, { 0, 0, 2 });
  verifyPinCounts(11, { 0, 0, 2 });
}

TEST_F(ADeltaPartitionedGraph, VerifyInitialKm1Gain) {
  verifyKm1Gain(0, { 0, 0, 0 });
  verifyKm1Gain(1, { 0, 0, -1 });
  verifyKm1Gain(2, { 0, 0, -1 });
  verifyKm1Gain(3, { 1, 0, 0 });
  verifyKm1Gain(4, { 1, 0, 2 });
  verifyKm1Gain(5, { -1, 0, 0 });
  verifyKm1Gain(6, { -1, 0, 0 });
}

TEST_F(ADeltaPartitionedGraph, MovesVertices) {
  delta_phg.changeNodePartWithGainCacheUpdate(1, 0, 1, 1000);
  ASSERT_EQ(0, phg.partID(1));
  ASSERT_EQ(1, delta_phg.partID(1));

  // Verify Pin Counts
  verifyPinCounts(0, { 1, 1, 0 });
  verifyPinCounts(2, { 1, 1, 0 });
  verifyPinCounts(1, { 0, 2, 0 });
  verifyPinCounts(5, { 0, 2, 0 });

  // Verify Move To Penalty
  verifyKm1Gain(1, { 0, 0, -1 });
  verifyKm1Gain(2, { 0, 2, 0 });
  verifyKm1Gain(4, { -1, 0, 1 });

  delta_phg.changeNodePartWithGainCacheUpdate(4, 1, 2, 1000);
  ASSERT_EQ(1, phg.partID(4));
  ASSERT_EQ(2, delta_phg.partID(4));

  // Verify Pin Counts
  verifyPinCounts(1, { 0, 1, 1 });
  verifyPinCounts(5, { 0, 1, 1 });
  verifyPinCounts(6, { 0, 0, 2 });
  verifyPinCounts(8, { 0, 0, 2 });
  verifyPinCounts(7, { 0, 0, 2 });
  verifyPinCounts(10, { 0, 0, 2 });

  // Verify Move To Penalty
  verifyKm1Gain(4, { -2, -1, 0 });
  verifyKm1Gain(1, { 1, 0, 1 });
  verifyKm1Gain(5, { -2, -2, 0 });
  verifyKm1Gain(6, { -2, -2, 0 });
}

} // namespace ds
} // namespace mt_kahypar