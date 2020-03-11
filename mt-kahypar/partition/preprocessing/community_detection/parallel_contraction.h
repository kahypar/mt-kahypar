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

#include <tbb/enumerable_thread_specific.h>

#include "kahypar/datastructure/sparse_map.h"

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/parallel/parallel_counting_sort.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
class ParallelContraction {
 private:
  static constexpr bool debug = false;

 public:
  using ArcWeight = ds::Graph::ArcWeight;
  using IncidentClusterWeights = kahypar::ds::SparseMap<PartitionID, ArcWeight>;

  static ds::Graph contract(const ds::Graph& GFine, ds::Clustering& C, size_t numTasks) {

    utils::Timer::instance().start_timer("compactify", "Compactify");
    size_t numClusters = C.compactify();
    utils::Timer::instance().stop_timer("compactify");
    auto nodes = boost::irange(NodeID(0), static_cast<NodeID>(GFine.numNodes()));

    // clang complains about structured bindings. we should revert once it works with clang
    utils::Timer::instance().start_timer("parallel_counting_sort", "Parallel Counting Sort");
    std::vector<NodeID> nodesSortedByCluster; std::vector<uint32_t> clusterBegin;
    std::tie(nodesSortedByCluster, clusterBegin) = parallel::ParallelCountingSort::sort(nodes, numClusters, C, numTasks);
    utils::Timer::instance().stop_timer("parallel_counting_sort");


    utils::Timer::instance().start_timer("construct_coarse_graph", "Construct Coarse Graph");
    tbb::enumerable_thread_specific<IncidentClusterWeights> ets_incidentClusterWeights(numClusters, 0);
    tbb::enumerable_thread_specific<size_t> ets_nArcs(0);
    ds::Graph GCoarse(numClusters);
    tbb::parallel_for_each(GCoarse.nodes(), [&](const NodeID coarseNode) {
        IncidentClusterWeights& incidentClusterWeights = ets_incidentClusterWeights.local();
        ArcWeight coarseNodeVolume = 0;

        /* accumulate weights to incident coarse nodes */
        size_t firstNext = clusterBegin[coarseNode + 1];
        for (size_t i = clusterBegin[coarseNode]; i < firstNext; ++i) {
          const NodeID fineNode = nodesSortedByCluster[i];
          for (auto& arc : GFine.arcsOf(fineNode)) {
            incidentClusterWeights[C[arc.head]] += arc.weight;
          }
          coarseNodeVolume += GFine.nodeVolume(fineNode);
        }
        for (auto& coarseNeighbor : incidentClusterWeights) {
          PartitionID communityID = coarseNeighbor.key;
          if (static_cast<NodeID>(communityID) != coarseNode)
            GCoarse.addDirectedEdge(coarseNode, static_cast<NodeID>(communityID), incidentClusterWeights[communityID]);
        }
        GCoarse.setNodeVolume(coarseNode, coarseNodeVolume);

        size_t& nArcs = ets_nArcs.local();
        nArcs += incidentClusterWeights.size();
        if (incidentClusterWeights.contains(static_cast<PartitionID>(coarseNode)))
          nArcs -= 1;

        incidentClusterWeights.clear();
      });
    utils::Timer::instance().stop_timer("construct_coarse_graph");

    GCoarse.setNumArcs(ets_nArcs.combine(std::plus<size_t>()));
    GCoarse.setTotalVolume(GFine.totalVolume());

    return GCoarse;
  }
};
}
