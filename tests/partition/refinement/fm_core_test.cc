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

#include <mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h>
#include <mt-kahypar/io/hypergraph_io.h>

namespace mt_kahypar {
namespace refinement {

using ::testing::Test;
class FMCoreTest : public Test {
public:
  FMCoreTest() {
    hg = io::readHypergraphFile<Hypergraph, HypergraphFactory>("../test_instances/ibm01.hgr", 0);
    phg = PartitionedHypergraph(k, hg);
    HypernodeID nodes_per_part = hg.initialNumNodes() / k;
    for (PartitionID i = 0; i < k; ++i) {
      for (HypernodeID u = i * nodes_per_part; u < (i+1) * nodes_per_part; ++u) {
        phg.setNodePart(u, i);
      }
    }
    phg.initializeGainInformation();

    sharedData = FMSharedData(hg.initialNumNodes(), hg.initialNumEdges(), k, 4);
    sharedData.setRemainingOriginalPins(phg);
  }

  Hypergraph hg;
  PartitionID k = 8;
  PartitionedHypergraph phg;
  FMSharedData sharedData;
};

TEST(FMCoreTest, PQInsertAndUpdate) {

}

}
}