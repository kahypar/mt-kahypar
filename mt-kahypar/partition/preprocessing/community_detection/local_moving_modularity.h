/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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



#include "mt-kahypar/datastructures/sparse_map.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"

#include "gtest/gtest_prod.h"

namespace mt_kahypar::metrics {
  double modularity(const Graph& graph, ds::Clustering& communities);
}

namespace mt_kahypar::community_detection {


class ParallelLocalMovingModularity {
 private:
  static constexpr bool advancedGainAdjustment = false;

  using AtomicArcWeight = parallel::AtomicWrapper<ArcWeight>;
  using LargeIncidentClusterWeights = ds::FixedSizeSparseMap<PartitionID, ArcWeight>;
  using CacheEfficientIncidentClusterWeights = ds::FixedSizeSparseMap<PartitionID, ArcWeight>;

 public:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  ParallelLocalMovingModularity(const Context& context,
                                         size_t numNodes,
                                         const bool disable_randomization = false) :
    _context(context),
    _max_degree(numNodes),
    _vertex_degree_sampling_threshold(context.preprocessing.community_detection.vertex_degree_sampling_threshold),
    _cluster_volumes(numNodes),
    _local_small_incident_cluster_weight(0),
    _local_large_incident_cluster_weight([&] {
      return construct_large_incident_cluster_weight_map();
    }),
    _disable_randomization(disable_randomization) { }

  ~ParallelLocalMovingModularity();

  bool localMoving(Graph& graph, ds::Clustering& communities);

 private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool ratingsFitIntoSmallSparseMap(const Graph& graph,
                                                                       const HypernodeID u)  {
    static constexpr size_t cache_efficient_map_size = CacheEfficientIncidentClusterWeights::MAP_SIZE / 3UL;
    return std::min(_vertex_degree_sampling_threshold, _max_degree) > cache_efficient_map_size &&
           graph.degree(u) <= cache_efficient_map_size;
  }

  LargeIncidentClusterWeights construct_large_incident_cluster_weight_map() {
    return LargeIncidentClusterWeights(3UL * std::min(_max_degree, _vertex_degree_sampling_threshold), 0);
  }

  // ! Only for testing
  void initializeClusterVolumes(const Graph& graph, ds::Clustering& communities);

  template<typename Map>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PartitionID computeMaxGainCluster(const Graph& graph,
                                                                       ds::Clustering& communities,
                                                                       const NodeID u,
                                                                       Map& incident_cluster_weights) {
    PartitionID from = communities[u];
    PartitionID bestCluster = communities[u];

    for (const Arc& arc : graph.arcsOf(u, _vertex_degree_sampling_threshold)) {
      incident_cluster_weights[communities[arc.head]] += arc.weight;
    }

    const ArcWeight volume_from = _cluster_volumes[from];
    const ArcWeight volU = graph.nodeVolume(u);
    const ArcWeight weight_from = incident_cluster_weights[from];

    const double volMultiplier = _vol_multiplier_div_by_node_vol * volU;
    double bestGain = weight_from - volMultiplier * (volume_from - volU);
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

    HEAVY_PREPROCESSING_ASSERT(verifyGain(graph, communities, u, bestCluster, bestGain, incident_cluster_weights));

    incident_cluster_weights.clear();

    return bestCluster;
  }


  inline double modularityGain(const ArcWeight weight_to,
                               const ArcWeight volume_to,
                               const double multiplier) {
    return weight_to - multiplier * volume_to;
    // missing term is - weight_from + multiplier * (volume_from - volume_node)
  }

  inline long double adjustAdvancedModGain(double gain,
                                           const ArcWeight weight_from,
                                           const ArcWeight volume_from,
                                           const ArcWeight volume_node) const {
    return 2.0L * _reciprocal_total_volume *
      (gain - weight_from + _reciprocal_total_volume *
        volume_node * (volume_from - volume_node));
  }

  template<typename Map>
  bool verifyGain(const Graph& graph, ds::Clustering& communities, NodeID u, PartitionID to, double gain, const Map& icw);

  static std::pair<ArcWeight, ArcWeight> intraClusterWeightsAndSumOfSquaredClusterVolumes(const Graph& graph, const ds::Clustering& communities);

  const Context& _context;
  size_t _max_degree;
  const size_t _vertex_degree_sampling_threshold;
  double _reciprocal_total_volume = 0.0;
  double _vol_multiplier_div_by_node_vol = 0.0;
  parallel::scalable_vector<AtomicArcWeight> _cluster_volumes;
  tbb::enumerable_thread_specific<CacheEfficientIncidentClusterWeights> _local_small_incident_cluster_weight;
  tbb::enumerable_thread_specific<LargeIncidentClusterWeights> _local_large_incident_cluster_weight;
  const bool _disable_randomization;


  FRIEND_TEST(ALouvain, ComputesMaxGainMove1);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove2);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove3);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove4);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove5);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove6);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove7);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove8);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove9);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove10);
};
}
