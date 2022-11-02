/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/refinement/fm/strategies/km1_gains.h"

using ::testing::Test;

namespace mt_kahypar {


  template<PartitionID K>
  class GainComputerTest : public Test {
  public:

    GainComputerTest() :
            hg(HypergraphFactory::construct(7, 4,
                                            {{0, 2},
                                                   {0, 1, 3, 4},
                                                   {3, 4, 6},
                                                   {2, 5, 6}})),
            context(),
            gain(K) {
      context.partition.k = K;
      context.partition.max_part_weights.assign(K, std::numeric_limits<HypernodeWeight>::max());
      phg = PartitionedHypergraph(K, hg, parallel_tag_t());
    }

    void assignPartitionIDs(const std::vector<PartitionID>& part_ids) {
      HypernodeID hn = 0;
      for (const PartitionID& part : part_ids) {
        ASSERT(part < K);
        phg.setNodePart(hn++, part);
      }
    }

    Hypergraph hg;
    PartitionedHypergraph phg;
    Context context;
    Km1GainComputer gain;
  };

  using Km1GainsK2 = GainComputerTest<2>;

  TEST_F(Km1GainsK2, ComputesCorrectMoveGainForVertex1) {
    assignPartitionIDs({1, 0, 0, 0, 0, 1, 1});
    auto [to, g] = gain.computeBestTargetBlock(phg, 0, context.partition.max_part_weights);
    ASSERT_EQ(0, to);
    ASSERT_EQ(2, g);
  }


  TEST_F(Km1GainsK2, ComputesCorrectMoveGainForVertex2) {
    assignPartitionIDs({0, 0, 0, 1, 0, 1, 1});
    auto [to, g] = gain.computeBestTargetBlock(phg, 3, context.partition.max_part_weights);
    ASSERT_EQ(0, to);
    ASSERT_EQ(1, g);
  }


  TEST_F(Km1GainsK2, ComputesCorrectMoveGainForVertex3) {
    assignPartitionIDs({0, 0, 0, 0, 0, 1, 1});
    auto [to, g] = gain.computeBestTargetBlock(phg, 4, context.partition.max_part_weights);
    ASSERT_EQ(1, to);
    ASSERT_EQ(-1, g);   // computeBestTarget block will select negative gain moves!
  }



  using Km1GainsK4 = GainComputerTest<4>;

  TEST_F(Km1GainsK4, ComputesCorrectMoveGainForVertex1) {
    assignPartitionIDs({0, 1, 2, 3, 3, 1, 2});
    auto [to, g] = gain.computeBestTargetBlock(phg, 0, context.partition.max_part_weights);

    gain.computeGains(phg, 0);
    ASSERT_EQ(gain.gains[1], 1);
    ASSERT_EQ(gain.gains[2], 1);
    ASSERT_EQ(gain.gains[3], 1);

    ASSERT_EQ(1, to); // all target blocks have equal weight --> take the first
    ASSERT_EQ(1, g);
  }


  TEST_F(Km1GainsK4, ComputesCorrectMoveGainForVertex2) {
    assignPartitionIDs({0, 3, 1, 2, 2, 0, 3});
    auto [to, g] = gain.computeBestTargetBlock(phg, 6, context.partition.max_part_weights);

    gain.computeGains(phg, 6);
    ASSERT_EQ(gain.gains[0], 1);
    ASSERT_EQ(gain.gains[1], 1);
    ASSERT_EQ(gain.gains[2], 1);
    ASSERT_EQ(1, to); // block 1 is lighter than block 0
    ASSERT_EQ(1, g);

    gain.clear();
    std::tie(to, g) = gain.computeBestTargetBlock(phg, 2, context.partition.max_part_weights);

    gain.computeGains(phg, 2);
    ASSERT_EQ(gain.gains[0], 2);
    ASSERT_EQ(gain.gains[2], 0);
    ASSERT_EQ(gain.gains[3], 1);

    ASSERT_EQ(to, 0);
    ASSERT_EQ(g, 2);
  }


  TEST_F(Km1GainsK4, ComputesCorrectMoveGainForVertex3) {
    assignPartitionIDs({0, 3, 1, 2, 2, 0, 3});
    auto [to, g] = gain.computeBestTargetBlock(phg, 3, context.partition.max_part_weights);

    gain.computeGains(phg, 3);
    ASSERT_EQ(gain.gains[0], -1);
    ASSERT_EQ(gain.gains[1], -2);
    ASSERT_EQ(gain.gains[3], 0);

    ASSERT_EQ(3, to);
    ASSERT_EQ(0, g);
  }
} // namespace