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

    context.partition.k = k;
    context.partition.epsilon = 0.03;
    context.setupPartWeights(hg.totalWeight());

    // ignore balance in this test
    for (PartitionID i = 0; i < k; ++i) {
      context.partition.max_part_weights[i] = hg.totalWeight();
    }
  }

  Hypergraph hg;
  PartitionID k = 8;
  PartitionedHypergraph phg;
  FMSharedData sharedData;
  Context context;
};

void printGains(PartitionedHypergraph& phg, PartitionID k) {
  for (HypernodeID u = 0; u < phg.initialNumNodes(); ++u) {
    std::cout << "u=" << u << "p=" << phg.partID(u) << ". gains=";
    for (PartitionID i = 0; i < k; ++i) {
      if (i != phg.partID(u)) {
        std::cout << phg.km1Gain(u, phg.partID(u), i) << " ";
      }
    }
    std::cout << std::endl;
  }
}

TEST_F(FMCoreTest, PQInsertAndUpdate) {
  //printGains(phg, k);
  LocalizedKWayFM fm(context, hg.initialNumNodes(), &sharedData.vertexPQHandles);
  fm.findMoves(phg, 23, sharedData, 0);
  for (MoveID move_id = 0; move_id < sharedData.moveTracker.numPerformedMoves(); ++move_id) {
    Move& m = sharedData.moveTracker.globalMoveOrder[move_id];
    LOG << V(move_id) << V(m.node) << V(m.from) << V(m.to) << V(m.gain);
  }

}

}
}