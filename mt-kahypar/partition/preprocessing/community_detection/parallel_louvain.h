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
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_contraction.h"
#include "mt-kahypar/partition/preprocessing/community_detection/plm.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
class ParallelModularityLouvain {
 private:
  static constexpr bool debug = false;

  static ds::Clustering localMovingContractRecurse(ds::Graph& fine_graph,
                                                   PLM& mlv,
                                                   const size_t num_tasks) {
    ds::Clustering communities(fine_graph.numNodes());

    DBG << "Start Local Moving";
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    bool communities_changed = mlv.localMoving(fine_graph, communities);
    utils::Timer::instance().stop_timer("local_moving");

    if (communities_changed) {
      utils::Timer::instance().start_timer("contraction", "Contraction");
      // Contract Communities
      ds::Graph coarse_graph = ParallelContraction::contract(fine_graph, communities, num_tasks);
      utils::Timer::instance().stop_timer("contraction");

      // Recurse on contracted graph
      ds::Clustering coarse_communities =
        localMovingContractRecurse(coarse_graph, mlv, num_tasks);

      utils::Timer::instance().start_timer("prolong", "Prolong");
      // Prolong Clustering
      for (NodeID u : fine_graph.nodes()) {
        communities[u] = coarse_communities[communities[u]];
      }
      utils::Timer::instance().stop_timer("prolong");
    }

    return communities;
  }

 public:
  static ds::Clustering run(ds::Graph& graph, const Context& context) {
    PLM mlv(context, graph.numNodes());
    ds::Clustering communities = localMovingContractRecurse(
      graph, mlv, context.shared_memory.num_threads);
    return communities;
  }
};
}
