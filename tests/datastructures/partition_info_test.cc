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

#include "mt-kahypar/datastructures/partition_info.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos0) {
  PartitionInfo part_info(1, 2);
  ASSERT_EQ(1, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos1) {
  PartitionInfo part_info(8, 8);
  ASSERT_EQ(1, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos2) {
  PartitionInfo part_info(8, 2);
  ASSERT_EQ(3, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos3) {
  PartitionInfo part_info(8, 4);
  ASSERT_EQ(2, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos4) {
  PartitionInfo part_info(16, 4);
  ASSERT_EQ(2, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos5) {
  PartitionInfo part_info(64, 2);
  ASSERT_EQ(6, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos6) {
  PartitionInfo part_info(64, 4);
  ASSERT_EQ(3, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos7) {
  PartitionInfo part_info(64, 8);
  ASSERT_EQ(2, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos8) {
  PartitionInfo part_info(64, 16);
  ASSERT_EQ(2, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectNumLocalBlockInfos9) {
  PartitionInfo part_info(64, 64);
  ASSERT_EQ(1, part_info.numLocalBlockInfos());
}

TEST(APartitionInfo, HasCorrectBlockSizeIfSetNodePartIsCalled1) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  ASSERT_EQ(1, part_info.partSize(2));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfSetNodePartIsCalled1) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  ASSERT_EQ(3, part_info.partWeight(2));
}

TEST(APartitionInfo, HasCorrectBlockSizeIfSetNodePartIsCalled2) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 3, 2);
  ASSERT_EQ(1, part_info.partSize(2));
  ASSERT_EQ(1, part_info.partSize(3));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfSetNodePartIsCalled2) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 3, 2);
  ASSERT_EQ(3, part_info.partWeight(2));
  ASSERT_EQ(2, part_info.partWeight(3));
}

TEST(APartitionInfo, HasCorrectBlockSizeIfSetNodePartIsCalled3) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  ASSERT_EQ(2, part_info.partSize(2));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfSetNodePartIsCalled3) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  ASSERT_EQ(5, part_info.partWeight(2));
}

TEST(APartitionInfo, HasCorrectBlockSizeIfChangeNodePartIsCalled1) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  ASSERT_EQ(1, part_info.partSize(2));
  ASSERT_EQ(1, part_info.partSize(0));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfChangeNodePartIsCalled1) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  ASSERT_EQ(3, part_info.partWeight(2));
  ASSERT_EQ(2, part_info.partWeight(0));
}

TEST(APartitionInfo, HasCorrectBlockSizeIfChangeNodePartIsCalled2) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  part_info.changeNodePart(0, 2, 1, 3);
  ASSERT_EQ(0, part_info.partSize(2));
  ASSERT_EQ(1, part_info.partSize(1));
  ASSERT_EQ(1, part_info.partSize(0));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfChangeNodePartIsCalled2) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  part_info.changeNodePart(0, 2, 1, 3);
  ASSERT_EQ(0, part_info.partWeight(2));
  ASSERT_EQ(3, part_info.partWeight(1));
  ASSERT_EQ(2, part_info.partWeight(0));
}

TEST(APartitionInfo, HasCorrectBlockSizeIfChangeNodePartIsCalled3) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  part_info.changeNodePart(0, 2, 0, 3);
  ASSERT_EQ(0, part_info.partSize(2));
  ASSERT_EQ(2, part_info.partSize(0));
}

TEST(APartitionInfo, HasCorrectBlockWeightIfChangeNodePartIsCalled3) {
  PartitionInfo part_info(64, 4);
  part_info.setNodePart(0, 2, 3);
  part_info.setNodePart(1, 2, 2);
  part_info.changeNodePart(1, 2, 0, 2);
  part_info.changeNodePart(0, 2, 0, 3);
  ASSERT_EQ(0, part_info.partWeight(2));
  ASSERT_EQ(5, part_info.partWeight(0));
}

using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;
using AtomicSize = parallel::IntegralAtomicWrapper<HypernodeID>;

void smokeTest(const HypernodeID num_hypernodes,
               const size_t num_threads,
               const PartitionID k) {
  std::vector<HypernodeWeight> weights(num_hypernodes, 0);
  for ( HypernodeID hn = 0; hn < num_hypernodes; ++hn ) {
    weights[hn] = utils::Randomize::instance().getRandomInt(1, 10, sched_getcpu());
  }

  PartitionInfo part_info(num_threads, k);
  std::vector<PartitionID> part_id(num_hypernodes, kInvalidPartition);
  std::vector<AtomicSize> part_size(k, AtomicSize(0));
  std::vector<AtomicWeight> part_weight(k, AtomicWeight(0));

  // Set Vertex Partition IDs
  tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID hn) {
    const PartitionID block = utils::Randomize::instance().getRandomInt(0, k-1, sched_getcpu());
    part_id[hn] = block;
    ++part_size[block];
    part_weight[block] += weights[hn];
    part_info.setNodePart(hn, block, weights[hn]);
  });

  for ( PartitionID block = 0; block < k; ++block ) {
    ASSERT_EQ(part_size[block].load(), part_info.partSize(block));
    ASSERT_EQ(part_weight[block].load(), part_info.partWeight(block));
  }

  // Change Vertex Partition IDs
  tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID hn) {
    const PartitionID from = part_id[hn];
    const PartitionID to = utils::Randomize::instance().getRandomInt(0, k-1, sched_getcpu());
    if ( from != to ) {
      --part_size[from];
      part_weight[from] -= weights[hn];
      ++part_size[to];
      part_weight[to] += weights[hn];
      part_info.changeNodePart(hn, from, to, weights[hn]);
    }
  });

  for ( PartitionID block = 0; block < k; ++block ) {
    ASSERT_EQ(part_size[block].load(), part_info.partSize(block));
    ASSERT_EQ(part_weight[block].load(), part_info.partWeight(block));
  }
}

TEST(APartitionInfo, SmokeTest1) {
  smokeTest(ID(100000), 64, 4);
}

TEST(APartitionInfo, SmokeTest2) {
  smokeTest(ID(100000), 64, 16);
}

TEST(APartitionInfo, SmokeTest3) {
  smokeTest(ID(100000), 64, 64);
}

} // namespace
} // namespace
