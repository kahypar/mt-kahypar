/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "local_moving_modularity.h"

#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar::metrics {
  double modularity(const Graph& graph, ds::Clustering& communities) {
    ASSERT(graph.canBeUsed());
    ASSERT(graph.numNodes() == communities.size());
    parallel::scalable_vector<parallel::AtomicWrapper<double>> internal_volume(graph.numNodes());
    parallel::scalable_vector<parallel::AtomicWrapper<double>> total_volume(graph.numNodes());
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      const PartitionID community_u = communities[u];
      ASSERT(community_u < static_cast<PartitionID>(graph.numNodes()));
      total_volume[community_u] += graph.nodeVolume(u);
      internal_volume[community_u] += graph.nodeVolume(u);
      for ( const Arc& arc : graph.arcsOf(u) ) {
        const NodeID v = arc.head;
        const PartitionID community_v = communities[v];
        ASSERT(community_v < static_cast<PartitionID>(graph.numNodes()));
        if ( community_u != community_v ) {
          internal_volume[community_u] -= arc.weight;
        }
      }
    });

    tbb::enumerable_thread_specific<double> local_modularity(0.0);
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      if ( total_volume[u].load(std::memory_order_relaxed) > 0.0 ) {
        local_modularity.local() += internal_volume[u] -
                                    (total_volume[u] * total_volume[u]) / graph.totalVolume();
      }
    });
    return local_modularity.combine(std::plus<>()) / graph.totalVolume();
  }
}

namespace mt_kahypar::community_detection {

  bool ParallelLocalMovingModularity::localMoving(Graph& graph, ds::Clustering& communities) {
    ASSERT(graph.canBeUsed());
    _max_degree = graph.max_degree();
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

    if ( graph.numArcs() > 0 ) {
      size_t number_of_nodes_moved = graph.numNodes();
      for (size_t currentRound = 0;
          number_of_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * graph.numNodes()
          && currentRound < _context.preprocessing.community_detection.max_pass_iterations;
          currentRound++) {

        if ( !_disable_randomization ) {
          utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
          utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
          utils::Timer::instance().stop_timer("random_shuffle");
        }

        tbb::enumerable_thread_specific<size_t> local_number_of_nodes_moved(0);
        auto moveNode = [&](const NodeID u) {
          const ArcWeight volU = graph.nodeVolume(u);
          const PartitionID from = communities[u];
          PartitionID best_cluster;

          if ( ratingsFitIntoSmallSparseMap(graph, u) ) {
            best_cluster = computeMaxGainCluster(graph, communities, u, _local_small_incident_cluster_weight.local());
          } else {
            LargeIncidentClusterWeights& large_incident_cluster_weight = _local_large_incident_cluster_weight.local();
            large_incident_cluster_weight.setMaxSize(3UL * std::min(_max_degree, _vertex_degree_sampling_threshold));
            best_cluster = computeMaxGainCluster(graph, communities, u, large_incident_cluster_weight);
          }

          if (best_cluster != from) {
            _cluster_volumes[best_cluster] += volU;
            _cluster_volumes[from] -= volU;
            communities[u] = best_cluster;
            ++local_number_of_nodes_moved.local();
          }
        };

        utils::Timer::instance().start_timer("local_moving_round", "Local Moving Round");

        #ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
        std::for_each(nodes.begin(), nodes.end(), moveNode);
        #else
        tbb::parallel_for(0UL, nodes.size(), [&](const size_t i) {
          moveNode(nodes[i]);
        });
        #endif

        number_of_nodes_moved = local_number_of_nodes_moved.combine(std::plus<>());
        clustering_changed |= number_of_nodes_moved > 0;
        utils::Timer::instance().stop_timer("local_moving_round");

        DBG << "Louvain-Pass #" << currentRound << " - Modularity:" << metrics::modularity(graph, communities);
      }
    }
    return clustering_changed;
  }

  template<typename Map>
  bool ParallelLocalMovingModularity::verifyGain(const Graph& graph,
                  ds::Clustering& communities,
                  const NodeID u,
                  const PartitionID to,
                  double gain,
                  const Map& icw) {
    const PartitionID from = communities[u];

    long double adjustedGain = adjustAdvancedModGain(gain, icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    const double volMultiplier = _vol_multiplier_div_by_node_vol * graph.nodeVolume(u);
    double modGain = modularityGain(icw.get(to), _cluster_volumes[to], volMultiplier);
    long double adjustedGainRecomputed = adjustAdvancedModGain(modGain, icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    unused(adjustedGainRecomputed);

    if (from == to) {
      adjustedGainRecomputed = 0.0L;
      adjustedGain = 0.0L;
    }

    ASSERT(adjustedGain == adjustedGainRecomputed);

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

    const bool result = math::are_almost_equal_ld(modBeforeMove + adjustedGain, modAfterMove, 1e-8);
    ASSERT(result,
           V(modBeforeMove + adjustedGain) << V(modAfterMove) << V(gain) << V(adjustedGain)
                                           << V(coverageBeforeMove) << V(expectedCoverageBeforeMove) << V(modBeforeMove)
                                           << V(coverageAfterMove) << V(expectedCoverageAfterMove) << V(modAfterMove));

    // revert move
    communities[u] = from;
    _cluster_volumes[to] -= graph.nodeVolume(u);
    _cluster_volumes[from] += graph.nodeVolume(u);

    return result;
  }

  std::pair<ArcWeight, ArcWeight> ParallelLocalMovingModularity::intraClusterWeightsAndSumOfSquaredClusterVolumes(
          const Graph& graph, const ds::Clustering& communities) {
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

  void ParallelLocalMovingModularity::initializeClusterVolumes(const Graph& graph, ds::Clustering& communities) {
    _reciprocal_total_volume = 1.0 / graph.totalVolume();
    _vol_multiplier_div_by_node_vol =  _reciprocal_total_volume;
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      const PartitionID community_id = communities[u];
      _cluster_volumes[community_id] += graph.nodeVolume(u);
    });
  }

  ParallelLocalMovingModularity::~ParallelLocalMovingModularity() {
    tbb::parallel_invoke([&] {
      parallel::parallel_free_thread_local_internal_data(
              _local_small_incident_cluster_weight, [&](CacheEfficientIncidentClusterWeights& data) {
                data.freeInternalData();
              });
    }, [&] {
      parallel::parallel_free_thread_local_internal_data(
              _local_large_incident_cluster_weight, [&](LargeIncidentClusterWeights& data) {
                data.freeInternalData();
              });
    }, [&] {
      parallel::free(_cluster_volumes);
    });
  }


}