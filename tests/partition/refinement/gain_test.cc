/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/fm/strategies/km1_gains.h"

using ::testing::Test;

namespace mt_kahypar {


  template<PartitionID K>
  class GainComputerTest : public Test {
  public:

    GainComputerTest() :
            hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
                                            7, 4,
                                            {{0, 2},
                                                   {0, 1, 3, 4},
                                                   {3, 4, 6},
                                                   {2, 5, 6}})),
            context(),
            gain(K) {
      context.partition.k = K;
      context.partition.max_part_weights.assign(K, std::numeric_limits<HypernodeWeight>::max());
      phg = PartitionedHypergraph(K, TBBNumaArena::GLOBAL_TASK_GROUP, hg);
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

    ASSERT_EQ(gain.gains[1], 1);
    ASSERT_EQ(gain.gains[2], 1);
    ASSERT_EQ(gain.gains[3], 1);

    ASSERT_EQ(1, to); // all target blocks have equal weight --> take the first
    ASSERT_EQ(1, g);
  }


  TEST_F(Km1GainsK4, ComputesCorrectMoveGainForVertex2) {
    assignPartitionIDs({0, 3, 1, 2, 2, 0, 3});
    auto [to, g] = gain.computeBestTargetBlock(phg, 6, context.partition.max_part_weights);
    ASSERT_EQ(gain.gains[0], 1);
    ASSERT_EQ(gain.gains[1], 1);
    ASSERT_EQ(gain.gains[2], 1);
    ASSERT_EQ(1, to); // block 1 is lighter than block 0
    ASSERT_EQ(1, g);

    std::tie(to, g) = gain.computeBestTargetBlock(phg, 2, context.partition.max_part_weights);

    ASSERT_EQ(gain.gains[0], 2);
    ASSERT_EQ(gain.gains[2], 0);
    ASSERT_EQ(gain.gains[3], 1);

    ASSERT_EQ(to, 0);
    ASSERT_EQ(g, 2);
  }


  TEST_F(Km1GainsK4, ComputesCorrectMoveGainForVertex3) {
    assignPartitionIDs({0, 3, 1, 2, 2, 0, 3});
    auto [to, g] = gain.computeBestTargetBlock(phg, 3, context.partition.max_part_weights);

    ASSERT_EQ(gain.gains[0], -1);
    ASSERT_EQ(gain.gains[1], -2);
    ASSERT_EQ(gain.gains[3], 0);

    ASSERT_EQ(3, to);
    ASSERT_EQ(0, g);
  }
} // namespace 