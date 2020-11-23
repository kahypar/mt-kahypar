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

#include "mt-kahypar/macros.h"

#include <mt-kahypar/definitions.h>
#include <mt-kahypar/io/hypergraph_io.h>

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(GainUpdates, Example1) {
  Hypergraph hg = io::readHypergraphFile("../tests/instances/twocenters.hgr", 0);
  PartitionID k = 2;
  PartitionedHypergraph<Hypergraph, HypergraphFactory> phg(k, hg);

  phg.setNodePart(0, 0);
  phg.setNodePart(1, 0);
  for (HypernodeID u = 4; u < 12; ++u) {
    phg.setNodePart(u, 0);
  }

  phg.setNodePart(2, 1);
  phg.setNodePart(3, 1);
  for (HypernodeID u = 12; u < 20; ++u) {
    phg.setNodePart(u, 1);
  }

  ASSERT_EQ(phg.partWeight(0), phg.partWeight(1));
  ASSERT_EQ(phg.partWeight(0), 10);

    phg.initializeGainCache();
  ASSERT_EQ(phg.km1Gain(0, phg.partID(0), 1), -1);
  ASSERT_EQ(phg.moveFromBenefit(0), 1);
  ASSERT_EQ(phg.moveToPenalty(0, 1), 2);

  ASSERT_EQ(phg.km1Gain(2, phg.partID(2), 0), -1);

  ASSERT_EQ(phg.km1Gain(4, phg.partID(4), 1), -1);
  ASSERT_EQ(phg.km1Gain(6, phg.partID(6), 1), -2);

  ASSERT_EQ(phg.km1Gain(12, phg.partID(12), 0), -1);
  ASSERT_EQ(phg.km1Gain(14, phg.partID(14), 0), -2);

    phg.changeNodePartWithGainCacheUpdate(8, 0, 1);

  phg.recomputeMoveFromBenefit(8);  // nodes are allowed to move once before moveFromBenefit must be recomputed
  ASSERT_EQ(phg.km1Gain(8, 1, 0), 2);

  ASSERT_EQ(phg.km1Gain(6, 0, 1), 0);
}


}  // namespace ds
}  // namespace mt_kahypar
