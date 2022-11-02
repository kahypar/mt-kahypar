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
  mt_kahypar::PartitionedHypergraph phg(k, hg);

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
  ASSERT_EQ(phg.moveFromPenalty(0), 2);
  ASSERT_EQ(phg.moveToBenefit(0, 1), 1);

  ASSERT_EQ(phg.km1Gain(2, phg.partID(2), 0), -1);

  ASSERT_EQ(phg.km1Gain(4, phg.partID(4), 1), -1);
  ASSERT_EQ(phg.km1Gain(6, phg.partID(6), 1), -2);

  ASSERT_EQ(phg.km1Gain(12, phg.partID(12), 0), -1);
  ASSERT_EQ(phg.km1Gain(14, phg.partID(14), 0), -2);

    phg.changeNodePartWithGainCacheUpdate(8, 0, 1);

  phg.recomputeMoveFromPenalty(8);  // nodes are allowed to move once before moveFromPenalty must be recomputed
  ASSERT_EQ(phg.km1Gain(8, 1, 0), 2);

  ASSERT_EQ(phg.km1Gain(6, 0, 1), 0);
}


}  // namespace ds
}  // namespace mt_kahypar
