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

#include "mt-kahypar/datastructures/clearlist.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/parallel/parallel_counting_sort.h"

namespace mt_kahypar {
class ParallelClusteringContractionAdjList {
 private:
  static constexpr bool debug = false;

 public:
  using ArcWeight = ds::AdjListGraph::ArcWeight;
  using IncidentClusterWeights = ds::ClearListMap<PartitionID, ArcWeight>;

  static ds::AdjListGraph contract(const ds::AdjListGraph& GFine, ds::Clustering& C, size_t numTasks) {
    auto t_compactify = tbb::tick_count::now();
    size_t numClusters = C.compactify();
    DBG << "compactify" << (tbb::tick_count::now() - t_compactify).seconds() << "[s]";

    auto t_bucket_assignment = tbb::tick_count::now();
    auto nodes = boost::irange(NodeID(0), static_cast<NodeID>(GFine.numNodes()));
    auto [nodesSortedByCluster, clusterBegin] = parallel::ParallelCountingSort::sort(nodes, numClusters, C, numTasks);
    DBG << "bucket assignment" << (tbb::tick_count::now() - t_bucket_assignment).seconds() << "[s]";

    auto t_edge_accumulation = tbb::tick_count::now();
    tbb::enumerable_thread_specific<IncidentClusterWeights> ets_incidentClusterWeights(numClusters, 0);
    tbb::enumerable_thread_specific<size_t> ets_nArcs(0);

    // measure edge accumulation load balancing
    tbb::enumerable_thread_specific<tbb::tick_count::interval_t> ts_runtime(0.0);
    tbb::enumerable_thread_specific<size_t> ts_deg(0);

    ds::AdjListGraph GCoarse(numClusters);
    tbb::parallel_for_each(GCoarse.nodes(), [&](const NodeID coarseNode) {
        // std::for_each(GCoarse.nodes().begin(), GCoarse.nodes().end(), [&](const NodeID coarseNode) {
        // if calls to .local() are too expensive, we could implement range-based ourselves. or use task_arena::thread_id() and OMP style allocation. I think .local() hashes thread IDs for some reason
        auto t_handle_cluster = tbb::tick_count::now();
        size_t& fineDegree = ts_deg.local();

        IncidentClusterWeights& incidentClusterWeights = ets_incidentClusterWeights.local();
        ArcWeight coarseNodeVolume = 0;

        /* accumulate weights to incident coarse nodes */
        size_t firstNext = clusterBegin[coarseNode + 1];
        for (size_t i = clusterBegin[coarseNode]; i < firstNext; ++i) {
          const NodeID fineNode = nodesSortedByCluster[i];
          fineDegree += GFine.degree(fineNode);
          for (auto& arc : GFine.arcsOf(fineNode)) {
            incidentClusterWeights.add(C[arc.head], arc.weight);
          }
          coarseNodeVolume += GFine.nodeVolume(fineNode);
        }
        for (PartitionID coarseNeighbor : incidentClusterWeights.keys())
          if (static_cast<NodeID>(coarseNeighbor) != coarseNode)
            GCoarse.addHalfEdge(coarseNode, static_cast<NodeID>(coarseNeighbor), incidentClusterWeights[coarseNeighbor]);

        GCoarse.setNodeVolume(coarseNode, coarseNodeVolume);

        size_t& nArcs = ets_nArcs.local();
        nArcs += incidentClusterWeights.keys().size();
        if (incidentClusterWeights.contains(static_cast<PartitionID>(coarseNode)))
          nArcs -= 1;

        incidentClusterWeights.clear();

        ts_runtime.local() += (tbb::tick_count::now() - t_handle_cluster);
      });

    GCoarse.setNumArcs(ets_nArcs.combine(std::plus<size_t>()));
    GCoarse.setTotalVolume(GFine.totalVolume());

    DBG << "edge accumulation time: " << (tbb::tick_count::now() - t_edge_accumulation).seconds() << "[s]";

    std::stringstream os;
    os << "edge accumulation thread degrees ";
    for (auto& d : ts_deg)
      os << d << " ";
    os << "\nedge accumulation thread runtimes ";
    for (auto& t : ts_runtime)
      os << t.seconds() << " ";
    DBG << os.str();

    return GCoarse;
  }
};

// TODO Parallel Contraction for CSR with ParallelCounting Sort. Single Pass Michi Style with indirections.
}
