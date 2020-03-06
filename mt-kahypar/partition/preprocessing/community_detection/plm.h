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

#include <atomic>

#include <tbb/enumerable_thread_specific.h>

#include "kahypar/datastructure/sparse_map.h"

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
class PLM {
 private:
  static constexpr bool advancedGainAdjustment = false;

  using Graph = ds::Graph;
  using ArcWeight = ds::Graph::ArcWeight;
  using AtomicArcWeight = parallel::AtomicWrapper<ArcWeight>;
  using Arc = ds::Graph::Arc;
  using IncidentClusterWeights = kahypar::ds::SparseMap<PartitionID, ArcWeight>;

 public:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = true;

  explicit PLM(const Context& context, size_t numNodes) :
    _context(context),
    _cluster_volumes(numNodes),
    _local_incident_cluster_weight(numNodes, 0) { }

  bool localMoving(Graph& graph, ds::Clustering& communities) {
    _reciprocal_total_volume = 1.0 / graph.totalVolume();
    _vol_multiplier_div_by_node_vol = _reciprocal_total_volume;

    parallel::scalable_vector<NodeID> nodes(graph.numNodes());
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      nodes[u] = u;
      communities[u] = u;
      _cluster_volumes[u].store(graph.nodeVolume(u));
    });

    // local moving
    bool clustering_changed = false;
    size_t number_of_nodes_moved = graph.numNodes();
    for (size_t currentRound = 0;
         number_of_nodes_moved >=
         _context.preprocessing.community_detection.min_eps_improvement * graph.numNodes() &&
         currentRound < _context.preprocessing.community_detection.max_pass_iterations; currentRound++) {

      utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
      utils::Randomize::instance().localizedParallelShuffleVector(
        nodes, 0UL, nodes.size(), _context.shared_memory.shuffle_block_size);
      utils::Timer::instance().stop_timer("random_shuffle");

      tbb::enumerable_thread_specific<size_t> local_number_of_nodes_moved(0);
      auto moveNode =
        [&](const NodeID u) {         // get rid of named lambda after testing?
          IncidentClusterWeights& incident_cluster_weights =
            _local_incident_cluster_weight.local();
          for (Arc& arc : graph.arcsOf(u)) {
            incident_cluster_weights.add(communities[arc.head], arc.weight);
          }

          PartitionID from = communities[u];
          PartitionID bestCluster = communities[u];

          const ArcWeight volume_from = _cluster_volumes[from];
          const ArcWeight volU = graph.nodeVolume(u);
          const ArcWeight weight_from = incident_cluster_weights[from];

          const double volMultiplier = _vol_multiplier_div_by_node_vol * volU;
          double bestGain = weight_from + volMultiplier * (volume_from - volU);
          for (const auto& clusterWeight : incident_cluster_weights) {
            PartitionID to = clusterWeight.key;
            // if from == to, we would have to remove volU from volume_to as well.
            // just skip it. it has (adjusted) gain zero.
            if (from == to) {
              continue;
            }

            const ArcWeight volume_to = _cluster_volumes[to],
              weight_to = incident_cluster_weights.get(to);

            double gain = modularityGain(weight_to, volume_to, volMultiplier);

            if (gain > bestGain) {
              bestCluster = to;
              bestGain = gain;
            }
          }

          HEAVY_PREPROCESSING_ASSERT(verifyGain(
            graph, communities, u, bestCluster, bestGain, incident_cluster_weights));

          incident_cluster_weights.clear();

          if (bestCluster != from) {
            _cluster_volumes[bestCluster] += volU;
            _cluster_volumes[from] -= volU;
            communities[u] = bestCluster;
            ++local_number_of_nodes_moved.local();
          }
        };

      utils::Timer::instance().start_timer("local_moving_round", "Local Moving Round");

      #ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
      std::for_each(nodes.begin(), nodes.end(), moveNode);
      #else
      tbb::parallel_for_each(nodes, moveNode);
      #endif

      number_of_nodes_moved = local_number_of_nodes_moved.combine(std::plus<size_t>());
      clustering_changed |= number_of_nodes_moved > 0;
      utils::Timer::instance().stop_timer("local_moving_round");

      DBG << V(currentRound) << V(number_of_nodes_moved);
    }
    return clustering_changed;
  }

 private:

  inline double modularityGain(const ArcWeight weight_to,
                               const ArcWeight volume_to,
                               const double multiplier) {
    return weight_to - multiplier * volume_to;
  }

  inline long double adjustAdvancedModGain(double gain,
                                           const ArcWeight weight_from,
                                           const ArcWeight volume_from,
                                           const ArcWeight volume_node) {
    return 2.0L * _reciprocal_total_volume *
      (gain - weight_from + _reciprocal_total_volume *
        volume_node * (volume_from - volume_node));
  }

  bool verifyGain(const Graph& graph,
                  ds::Clustering& communities,
                  const NodeID u,
                  const PartitionID to,
                  double gain,
                  const IncidentClusterWeights& icw) {
    const PartitionID from = communities[u];

    long double adjustedGain = adjustAdvancedModGain(
      gain, icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    const double volMultiplier = _vol_multiplier_div_by_node_vol * graph.nodeVolume(u);
    long double adjustedGainRecomputed = adjustAdvancedModGain(
      modularityGain(icw.get(to), _cluster_volumes[to], volMultiplier),
      icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    unused(adjustedGainRecomputed);

    if (from == to) {
      adjustedGainRecomputed = 0.0L;
      adjustedGain = 0.0L;
    }

    ASSERT(adjustedGain == adjustedGainRecomputed);

    auto eq = [&](const long double x, const long double y) {
                static constexpr double eps = 1e-8;
                long double diff = x - y;
                if (std::abs(diff) >= eps) {
                  LOG << V(x) << V(y) << V(diff);
                }
                return std::abs(diff) < eps;
              };

    long double dTotalVolumeSquared = static_cast<long double>(graph.totalVolume()) * static_cast<long double>(graph.totalVolume());

    auto accBeforeMove = intraClusterWeightsAndSumOfSquaredClusterVolumes(graph, communities);
    long double coverageBeforeMove = static_cast<long double>(accBeforeMove.first) / graph.totalVolume();
    long double expectedCoverageBeforeMove = accBeforeMove.second / dTotalVolumeSquared;
    long double modBeforeMove = coverageBeforeMove - expectedCoverageBeforeMove;

    // apply move
    communities[u] = to;
    _cluster_volumes[to] += graph.nodeVolume(u);
    _cluster_volumes[from] -= graph.nodeVolume(u);

    auto accAfterMove = intraClusterWeightsAndSumOfSquaredClusterVolumes(graph, communities);
    long double coverageAfterMove = static_cast<long double>(accAfterMove.first) / graph.totalVolume();
    long double expectedCoverageAfterMove = accAfterMove.second / dTotalVolumeSquared;
    long double modAfterMove = coverageAfterMove - expectedCoverageAfterMove;

    bool comp = eq(modBeforeMove + adjustedGain, modAfterMove);
    ASSERT(comp,
           V(modBeforeMove + adjustedGain) << V(modAfterMove) << V(gain) << V(adjustedGain)
           << V(coverageBeforeMove) << V(expectedCoverageBeforeMove) << V(modBeforeMove)
           << V(coverageAfterMove) << V(expectedCoverageAfterMove) << V(modAfterMove));

    // revert move
    communities[u] = from;
    _cluster_volumes[to] -= graph.nodeVolume(u);
    _cluster_volumes[from] += graph.nodeVolume(u);

    return comp;
  }

  static std::pair<ArcWeight, ArcWeight> intraClusterWeightsAndSumOfSquaredClusterVolumes(const Graph& graph, const ds::Clustering& communities) {
    ArcWeight intraClusterWeights = 0;
    ArcWeight sumOfSquaredClusterVolumes = 0;
    std::vector<ArcWeight> _cluster_volumes(graph.numNodes(), 0);

    for (NodeID u : graph.nodes()) {
      ArcWeight arcVol = 0;
      for (const Arc& arc : graph.arcsOf(u)) {
        if (communities[u] == communities[arc.head])
          intraClusterWeights += arc.weight;
        arcVol += arc.weight;
      }

      ArcWeight selfLoopWeight = graph.nodeVolume(u) - arcVol;          // already accounted for as twice!
      ASSERT(selfLoopWeight >= 0.0);
      intraClusterWeights += selfLoopWeight;
      _cluster_volumes[communities[u]] += graph.nodeVolume(u);
    }

    for (NodeID cluster : graph.nodes()) {      // unused cluster IDs have volume 0
      sumOfSquaredClusterVolumes += _cluster_volumes[cluster] * _cluster_volumes[cluster];
    }
    return std::make_pair(intraClusterWeights, sumOfSquaredClusterVolumes);
  }

  const Context& _context;
  double _reciprocal_total_volume = 0.0;
  double _vol_multiplier_div_by_node_vol = 0.0;
  std::vector<AtomicArcWeight> _cluster_volumes;
  tbb::enumerable_thread_specific<IncidentClusterWeights> _local_incident_cluster_weight;
};
}
