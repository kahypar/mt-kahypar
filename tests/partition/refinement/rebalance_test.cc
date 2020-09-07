/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include <functional>
#include <random>


#include "gmock/gmock.h"

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"



using ::testing::Test;

namespace mt_kahypar {

TEST(RebalanceTests, HeapSortWithMoveGainComparator) {
  vec<Move> moves;
  vec<Gain> gains = { 52, 12, 72, -154, 2672, -717, 1346, -7111, -113461, 136682, 3833 };

  for (HypernodeID i = 0; i < gains.size(); ++i) {
    moves.push_back(Move{-1, -1 , i, gains[i]});
  }

  std::make_heap(moves.begin(), moves.end(), Km1Rebalancer::MoveGainComparator());
  for (size_t i = 0; i < moves.size(); ++i) {
    std::pop_heap(moves.begin(), moves.end() - i, Km1Rebalancer::MoveGainComparator());
  }

  // assert that moves is sorted descendingly
  ASSERT_TRUE(std::is_sorted(moves.begin(), moves.end(),
                         [](const Move& lhs, const Move& rhs) {
    return lhs.gain > rhs.gain;
  }) );
}

TEST(RebalanceTests, FindsMoves) {
  PartitionID k = 8;
  Context context;
  context.partition.k = k;
  context.partition.epsilon = 0.03;
  TaskGroupID task_group_id = 0;
  Hypergraph hg = io::readHypergraphFile("../tests/instances/contracted_ibm01.hgr", 0, true /* enable stable construction */);
  context.setupPartWeights(hg.totalWeight());
  PartitionedHypergraph phg = PartitionedHypergraph(k, hg);

  HypernodeID nodes_per_part = hg.initialNumNodes() / (k-4);
  ASSERT(hg.initialNumNodes() % (k - 4) == 0);
  for (PartitionID i = 0; i < k - 4; ++i) {
    for (HypernodeID u = i * nodes_per_part; u < (i+1) * nodes_per_part; ++u) {
      phg.setOnlyNodePart(u, i);
    }
  }
  phg.initializePartition(task_group_id);
    phg.initializeGainCache();

  Km1Rebalancer rebalancer(phg, context);
  vec<Move> moves_to_empty_blocks = rebalancer.repairEmptyBlocks();

  ASSERT_EQ(moves_to_empty_blocks.size(), 4);

  for (Move& m : moves_to_empty_blocks) {
    ASSERT_EQ(phg.km1Gain(m.node, m.from, m.to), m.gain);
    Gain recomputed_gain = phg.moveFromBenefitRecomputed(m.node) - phg.moveToPenaltyRecomputed(m.node, m.to);
    if (recomputed_gain == 0) {
      ASSERT_TRUE([&]() {
        for (HyperedgeID e : phg.incidentEdges(m.node)) {
          if (phg.pinCountInPart(e, m.from) != 1) {
            return false;
          }
        }
        return true;
      }());
    }
    ASSERT_EQ(m.gain, recomputed_gain);
    ASSERT_EQ(m.from, phg.partID(m.node));
    ASSERT_EQ(phg.partWeight(m.to), 0);
    ASSERT_GE(m.to, k - 4);
    ASSERT_LT(m.from, k - 4);
    phg.changeNodePartWithGainCacheUpdate(m.node, m.from, m.to);
  }

  moves_to_empty_blocks = rebalancer.repairEmptyBlocks();
  ASSERT_EQ(moves_to_empty_blocks.size(), 0);
}

}