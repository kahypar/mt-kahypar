/*******************************************************************************
 * This file is part of MT-KaHyPar.
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

#include <mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h>

namespace mt_kahypar {

class InitialPartitioningBoundaryFM {
public:
  Gain refine(PartitionedHypergraph& phg, const Context& context) {
    phg.initializeGainInformation();
    ASSERT(phg.k() == 2 && context.partition.k == 2);
    FMSharedData sharedData(phg.initialNumNodes(), context);
    LocalizedKWayFM fm(context, phg.initialNumNodes(), sharedData.vertexPQHandles.data());

    vec<HypernodeID> initialNodes;
    for (HypernodeID u : phg.nodes()) {
      if (phg.isBorderNode(u)) {
        initialNodes.push_back(u);
      }
    }

    fm.findMoves(phg, sharedData, initialNodes);

    return fm.stats.estimated_improvement;  // since the call is single-threaded, the estimated improvement is exact
  }
};

}