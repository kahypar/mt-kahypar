/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "parallel_louvain.h"

#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::community_detection {

  ds::Clustering local_moving_contract_recurse(Graph& fine_graph, ParallelLocalMovingModularity& mlv, const Context& context) {
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    ds::Clustering communities(fine_graph.numNodes());
    bool communities_changed = mlv.localMoving(fine_graph, communities);
    utils::Timer::instance().stop_timer("local_moving");

    if (communities_changed) {
      utils::Timer::instance().start_timer("contraction", "Contraction");
      // Contract Communities
      Graph coarse_graph = fine_graph.contract(communities, context.preprocessing.community_detection.low_memory_contraction);
      ASSERT(coarse_graph.totalVolume() == fine_graph.totalVolume());
      utils::Timer::instance().stop_timer("contraction");

      // Recurse on contracted graph
      ds::Clustering coarse_communities = local_moving_contract_recurse(coarse_graph, mlv, context);

      utils::Timer::instance().start_timer("project", "Project");
      // Prolong Clustering
      tbb::parallel_for(0UL, fine_graph.numNodes(), [&](const NodeID u) {
        ASSERT(communities[u] < static_cast<PartitionID>(coarse_communities.size()));
        communities[u] = coarse_communities[communities[u]];
      });
      utils::Timer::instance().stop_timer("project");
    }

    return communities;
  }

  ds::Clustering run_parallel_louvain(Graph& graph, const Context& context, bool disable_randomization) {
    ParallelLocalMovingModularity mlv(context, graph.numNodes(), disable_randomization);
    ds::Clustering communities = local_moving_contract_recurse(graph, mlv, context);
    return communities;
  }
}