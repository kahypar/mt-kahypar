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
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    bool clustering_changed = mlv.localMoving(GFine, C);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().update_timing("local_moving", "Local Moving",
                                                       "perform_community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING,
                                                       std::chrono::duration<double>(end - start).count());


/*

        DBG << "Exiting so we only test local moving";
        std::exit(-1);
*/

    if (clustering_changed) {
      // contract
      DBG << "Contract";

      start = std::chrono::high_resolution_clock::now();
      ds::AdjListGraph GCoarse = ParallelClusteringContractionAdjList::contract(GFine, C, numTasks);
      end = std::chrono::high_resolution_clock::now();
      mt_kahypar::utils::Timer::instance().update_timing("contraction", "Contraction",
                                                         "perform_community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING,
                                                         std::chrono::duration<double>(end - start).count());


#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
      ds::Clustering coarseGraphSingletons(GCoarse.numNodes());
      coarseGraphSingletons.assignSingleton();
      // assert(PLM::doubleMod(GFine, PLM::intraClusterWeights_And_SumOfSquaredClusterVolumes(GFine, C)) == PLM::integerModularityFromScratch(GCoarse, coarseGraphSingletons));
#endif

      ClusteringStatistics::printLocalMovingStats(GFine, C);

      // recurse
      ds::Clustering coarseC = localMovingContractRecurse(GCoarse, mlv, numTasks);

      start = std::chrono::high_resolution_clock::now();
      // prolong clustering
      for (NodeID u : GFine.nodes()) // parallelize
        C[u] = coarseC[C[u]];
      end = std::chrono::high_resolution_clock::now();
      mt_kahypar::utils::Timer::instance().update_timing("prolong", "Prolong",
                                                         "perform_community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING,
                                                         std::chrono::duration<double>(end - start).count());
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
