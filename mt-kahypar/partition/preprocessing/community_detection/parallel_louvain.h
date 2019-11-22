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

namespace mt_kahypar {
class ParallelModularityLouvain {
 private:
  static constexpr bool debug = false;

  static ds::Clustering localMovingContractRecurse(ds::AdjListGraph& GFine, PLM& mlv, size_t numTasks) {
    ds::Clustering C(GFine.numNodes());

    DBG << "Start Local Moving";
    auto t_lm = tbb::tick_count::now();
    bool clustering_changed = mlv.localMoving(GFine, C);
    mlv.tr.report("Local Moving", tbb::tick_count::now() - t_lm);

/*

        DBG << "Exiting so we only test local moving";
        std::exit(-1);
*/

    if (clustering_changed) {
      // contract
      DBG << "Contract";

      auto t_contract = tbb::tick_count::now();
      ds::AdjListGraph GCoarse = ParallelClusteringContractionAdjList::contract(GFine, C, numTasks);
      mlv.tr.report("Contraction", tbb::tick_count::now() - t_contract);

#ifndef NDEBUG
      ds::Clustering coarseGraphSingletons(GCoarse.numNodes());
      coarseGraphSingletons.assignSingleton();
      // assert(PLM::doubleMod(GFine, PLM::intraClusterWeights_And_SumOfSquaredClusterVolumes(GFine, C)) == PLM::integerModularityFromScratch(GCoarse, coarseGraphSingletons));
#endif

      ClusteringStatistics::printLocalMovingStats(GFine, C, mlv.tr);

      // recurse
      ds::Clustering coarseC = localMovingContractRecurse(GCoarse, mlv, numTasks);

      auto t_prolong = tbb::tick_count::now();
      // prolong clustering
      for (NodeID u : GFine.nodes()) // parallelize
        C[u] = coarseC[C[u]];
      mlv.tr.report("Prolong", tbb::tick_count::now() - t_prolong);
      // assert(PLM::integerModularityFromScratch(GFine, C) == PLM::integerModularityFromScratch(GCoarse, coarseC));
    }

    return C;
  }

 public:
  static ds::Clustering run(ds::AdjListGraph& graph, const Context& context) {
    PLM mlv(context, graph.numNodes());
    ds::Clustering C = localMovingContractRecurse(graph, mlv, context.shared_memory.num_threads);
    ClusteringStatistics::printLocalMovingStats(graph, C, mlv.tr);
    return C;
  }
};
}
