/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (communities) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/community_detection/plm.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
class ParallelModularityLouvain {
 private:
  static constexpr bool debug = false;

  static ds::Clustering localMovingContractRecurse(Graph& fine_graph,
                                                   PLM& mlv,
                                                   const size_t num_tasks) {
    DBG << V(fine_graph.numNodes())
        << V(fine_graph.numArcs())
        << V(fine_graph.totalVolume());

    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    ds::Clustering communities(fine_graph.numNodes());
    bool communities_changed = mlv.localMoving(fine_graph, communities);
    utils::Timer::instance().stop_timer("local_moving");

    if (communities_changed) {
      DBG << "Modularity after PLM:" << metrics::modularity(fine_graph, communities);
      utils::Timer::instance().start_timer("contraction", "Contraction");
      // Contract Communities
      Graph coarse_graph = fine_graph.contract(communities);
      ASSERT(coarse_graph.totalVolume() == fine_graph.totalVolume(),
        V(coarse_graph.totalVolume()) << V(fine_graph.totalVolume()));
      utils::Timer::instance().stop_timer("contraction");
      DBG << "Modularity after Contraction:" << metrics::modularity(fine_graph, communities);

      // Recurse on contracted graph
      ds::Clustering coarse_communities =
        localMovingContractRecurse(coarse_graph, mlv, num_tasks);

      utils::Timer::instance().start_timer("prolong", "Prolong");
      // Prolong Clustering
      for (NodeID u : fine_graph.nodes()) {
        ASSERT(communities[u] < static_cast<PartitionID>(coarse_communities.size()));
        communities[u] = coarse_communities[communities[u]];
      }
      utils::Timer::instance().stop_timer("prolong");

      DBG << "Modularity after Prolong:" << metrics::modularity(fine_graph, communities);
    }

    return communities;
  }

 public:
  static ds::Clustering run(Graph& graph,
                            const Context& context,
                            const size_t num_tasks,
                            const bool disable_randomization = false) {
    PLM mlv(context, graph.numNodes(), disable_randomization);
    ds::Clustering communities = localMovingContractRecurse(graph, mlv, num_tasks);
    return communities;
  }
};
}
