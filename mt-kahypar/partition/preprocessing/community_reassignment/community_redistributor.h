/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/parallel_for.h"
#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace preprocessing {
template <typename TypeTraits>
class CommunityRedistributorT {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool enable_heavy_assert = false;

 public:
  CommunityRedistributorT(const CommunityRedistributorT&) = delete;
  CommunityRedistributorT & operator= (const CommunityRedistributorT &) = delete;
  CommunityRedistributorT(CommunityRedistributorT&&) = delete;
  CommunityRedistributorT & operator= (CommunityRedistributorT &&) = delete;

  static HyperGraph redistribute(const TaskGroupID task_group_id,
                                 const HyperGraph& hg,
                                 const parallel::scalable_vector<PartitionID>& community_assignment) {
    // Compute Node Mapping
    utils::Timer::instance().start_timer("compute_node_mapping", "Compute Node Mapping");
    parallel::scalable_vector<int> vertices_to_numa_node(hg.initialNumNodes(), -1);
    tbb::parallel_for(0UL, hg.initialNumNodes(), [&](const HypernodeID& hn) {
          vertices_to_numa_node[hn] = community_assignment[hg.communityID(hg.globalNodeID(hn))];
        });
    utils::Timer::instance().stop_timer("compute_node_mapping");

    utils::Timer::instance().start_timer("redistribute_hypergraph", "Redistribute Hypergraph");
    HyperGraph redistributed_hypergraph = HyperGraphFactory::construct(
      task_group_id, hg, std::move(vertices_to_numa_node));
    utils::Timer::instance().stop_timer("redistribute_hypergraph");

    utils::Timer::instance().start_timer("setup_communities", "Setup Community Structure");
    redistributed_hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      const HypernodeID original_id = redistributed_hypergraph.originalNodeID(hn);
      redistributed_hypergraph.setCommunityID(hn,
        hg.communityID(hg.globalNodeID(original_id)));
    });
    redistributed_hypergraph.initializeCommunities(task_group_id);
    utils::Timer::instance().stop_timer("setup_communities");

    HEAVY_PREPROCESSING_ASSERT([&] {
          for (const HypernodeID& hn : redistributed_hypergraph.nodes()) {
            int node = common::get_numa_node_of_vertex(hn);
            PartitionID community_id = redistributed_hypergraph.communityID(hn);
            if (community_assignment[community_id] != node) {
              LOG << "Hypernode" << hn << "should be on numa node" << community_assignment[community_id]
                  << "but is on node" << node;
              return false;
            }
          }
          return true;
        } (), "There are verticies assigned to wrong numa node");

    return redistributed_hypergraph;
  }

 protected:
  CommunityRedistributorT() = default;
};

using CommunityRedistributor = CommunityRedistributorT<GlobalTypeTraits>;
}  // namespace preprocessing
}  // namespace mt_kahypar
