/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma once

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_contraction.h"
#include "mt-kahypar/partition/preprocessing/community_detection/plm.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
class ParallelModularityLouvain {
 private:
  static constexpr bool debug = false;

  static ds::Clustering localMovingContractRecurse(ds::AdjListGraph& GFine, PLM& mlv, size_t numTasks) {
    ds::Clustering C(GFine.numNodes());

    DBG << "Start Local Moving";
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    bool clustering_changed = mlv.localMoving(GFine, C);
    utils::Timer::instance().stop_timer("local_moving");

/*

        ERROR("Exiting so we only test local moving");
*/

    if (clustering_changed) {
      // contract
      DBG << "Contract";

      utils::Timer::instance().start_timer("contraction", "Contraction");
      ds::AdjListGraph GCoarse = ParallelClusteringContractionAdjList::contract(GFine, C, numTasks);
      utils::Timer::instance().stop_timer("contraction");

#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
      ds::Clustering coarseGraphSingletons(GCoarse.numNodes());
      coarseGraphSingletons.assignSingleton();
      // assert(PLM::doubleMod(GFine, PLM::intraClusterWeights_And_SumOfSquaredClusterVolumes(GFine, C)) == PLM::integerModularityFromScratch(GCoarse, coarseGraphSingletons));
#endif

      ClusteringStatistics::printLocalMovingStats(GFine, C);

      // recurse
      ds::Clustering coarseC = localMovingContractRecurse(GCoarse, mlv, numTasks);

      utils::Timer::instance().start_timer("prolong", "Prolong");
      // prolong clustering
      for (NodeID u : GFine.nodes()) // parallelize
        C[u] = coarseC[C[u]];
      utils::Timer::instance().stop_timer("prolong");
      // assert(PLM::integerModularityFromScratch(GFine, C) == PLM::integerModularityFromScratch(GCoarse, coarseC));
    }

    return C;
  }

 public:
  static ds::Clustering run(ds::AdjListGraph& graph, const Context& context) {
    PLM mlv(context, graph.numNodes());
    ds::Clustering C = localMovingContractRecurse(graph, mlv, context.shared_memory.num_threads);
    ClusteringStatistics::printLocalMovingStats(graph, C);
    return C;
  }
};
}
