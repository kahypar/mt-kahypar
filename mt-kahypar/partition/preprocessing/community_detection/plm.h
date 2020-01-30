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

#include <atomic>

#include <tbb/enumerable_thread_specific.h>

#include "kahypar/datastructure/sparse_map.h"

#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/partition/preprocessing/community_detection/clustering_statistics.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
class PLM {
 private:
  static constexpr bool advancedGainAdjustment = false;

  using Graph = ds::AdjListGraph;
  using ArcWeight = ds::AdjListGraph::ArcWeight;
  using AtomicArcWeight = parallel::AtomicWrapper<ArcWeight>;
  using Arc = ds::AdjListGraph::Arc;
  using IncidentClusterWeights = kahypar::ds::SparseMap<PartitionID, ArcWeight>;

 public:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  // multiply with expected coverage part
  // static constexpr double resolutionGamma = 1.0;
  // static constexpr int resolutionGamma = 1;

  explicit PLM(const Context& context, size_t numNodes) :
    _context(context),
    clusterVolumes(numNodes),
    ets_incidentClusterWeights(numNodes, 0) { }

  bool localMoving(Graph& G, ds::Clustering& C) {
    reciprocalTotalVolume = 1.0 / G.totalVolume();
    totalVolume = G.totalVolume();
    volMultiplierDivByNodeVol = reciprocalTotalVolume;          // * resolutionGamma;

    ArcWeight maxAllowedClusterVolume = totalVolume;
    if (_context.preprocessing.community_detection.load_balancing_strategy ==
        CommunityLoadBalancingStrategy::size_constraint) {
      maxAllowedClusterVolume = std::ceil(maxAllowedClusterVolume /
                                          ((double)(_context.preprocessing.community_detection.size_constraint_factor *
                                                    _context.shared_memory.num_threads)));
    }

    parallel::scalable_vector<NodeID> nodes(G.numNodes());
    for (NodeID u : G.nodes()) {
      nodes[u] = u;
      clusterVolumes[u].store(G.nodeVolume(u));
      // clusterVolumes[u] = G.nodeVolume(u);
    }
    /*
    tbb::parallel_for(G.nodesParallelCoarseChunking(), [&](const NodeID u) {		//set coarse grain size
        clusterVolumes[u] = G.nodeVolume(u);
        nodes[u] = u;								//iota
    });
    */
    C.assignSingleton();

    // local moving
    bool clusteringChanged = false;
    size_t nodesMovedThisRound = G.numNodes();
    for (size_t currentRound = 0;
         nodesMovedThisRound >= _context.preprocessing.community_detection.min_eps_improvement * G.numNodes() &&
         currentRound < _context.preprocessing.community_detection.max_pass_iterations; currentRound++) {
      // parallel shuffle starts becoming competitive with sequential shuffle at four cores... :(
      // TODO implement block-based weak shuffling or use the pseudo-random online permutation approach
      utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
      utils::Randomize::instance().localizedParallelShuffleVector(
        nodes, 0UL, nodes.size(), _context.shared_memory.shuffle_block_size);
      utils::Timer::instance().stop_timer("random_shuffle");

      tbb::enumerable_thread_specific<size_t> ets_nodesMovedThisRound(0);

      tbb::enumerable_thread_specific<tbb::tick_count::interval_t> ts_runtime(0.0);
      tbb::enumerable_thread_specific<size_t> ts_deg(0);            // to measure potential imbalance. measurements don't show any imbalance

      auto moveNode =
        [&](const NodeID u) {         // get rid of named lambda after testing?
          auto t_node_move = tbb::tick_count::now();

          IncidentClusterWeights& incidentClusterWeights = ets_incidentClusterWeights.local();
          for (Arc& arc : G.arcsOf(u))
            incidentClusterWeights.add(C[arc.head], arc.weight);

          PartitionID from = C[u];
          PartitionID bestCluster = C[u];

          const ArcWeight volFrom = clusterVolumes[from],
            volU = G.nodeVolume(u),
            weightFrom = incidentClusterWeights[from];

          const double volMultiplier = volMultiplierDivByNodeVol * volU;

          // double bestGain = 0.0;												//basic gain adjustment
          double bestGain = weightFrom + volMultiplier * (volFrom - volU);                // advanced gain adjustment

          for (const auto& clusterWeight : incidentClusterWeights) {
            PartitionID to = clusterWeight.key;
            if (from == to) // if from == to, we would have to remove volU from volTo as well. just skip it. it has (adjusted) gain zero.
              continue;

            const ArcWeight volTo = clusterVolumes[to],
              weightTo = incidentClusterWeights.get(to);

            if (volU + volTo > maxAllowedClusterVolume) {
              continue;
            }

            // double gain = modularityGain(weightFrom, weightTo, volFrom, volTo, volU);
            double gain = modularityGain(weightTo, volTo, volMultiplier);

            if (gain > bestGain) {
              bestCluster = to;
              bestGain = gain;
            }
          }

          HEAVY_PREPROCESSING_ASSERT(verifyGain(G, C, u, bestCluster, bestGain, incidentClusterWeights));

          incidentClusterWeights.clear();

          if (bestCluster != from) {
            clusterVolumes[bestCluster] += volU;
            clusterVolumes[from] -= volU;
            C[u] = bestCluster;
            ets_nodesMovedThisRound.local()++;
          }

          ts_deg.local() += G.degree(u);
          ts_runtime.local() += (tbb::tick_count::now() - t_node_move);
        };

      utils::Timer::instance().start_timer("local_moving_round", "Local Moving Round");

#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
      std::for_each(nodes.begin(), nodes.end(), moveNode);
#else
      tbb::parallel_for_each(nodes, moveNode);
#endif

      nodesMovedThisRound = ets_nodesMovedThisRound.combine(std::plus<size_t>());
      clusteringChanged |= nodesMovedThisRound > 0;

      if (debug) {
        std::stringstream os;
        os << "Thread Local Runtime : ";
        for (auto& t : ts_runtime) {
          os << t.seconds() << " ";
        }
        os << std::endl << "Thread Local Degree Sum: ";
        for (auto& d : ts_deg)
          os << d << " ";
        DBG << os.str();
      }

      utils::Timer::instance().stop_timer("local_moving_round");

      DBG << V(currentRound) << V(nodesMovedThisRound);
    }
    return clusteringChanged;
  }

  bool verifyGain(const Graph& G, ds::Clustering& C, const NodeID u, const PartitionID to, double gain, const IncidentClusterWeights& icw) {
#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
    const PartitionID from = C[u];

    // long double adjustedGain = adjustBasicModGain(gain);
    // long double adjustedGainRecomputed = adjustBasicModGain(modularityGain(icw.get(from), icw.get(to), clusterVolumes[from], clusterVolumes[to], G.nodeVolume(u)));
    long double adjustedGain = adjustAdvancedModGain(gain, icw.get(from), clusterVolumes[from], G.nodeVolume(u));
    const double volMultiplier = volMultiplierDivByNodeVol * G.nodeVolume(u);
    long double adjustedGainRecomputed = adjustAdvancedModGain(modularityGain(icw.get(to), clusterVolumes[to], volMultiplier), icw.get(from), clusterVolumes[from], G.nodeVolume(u));

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

    long double dTotalVolumeSquared = static_cast<long double>(G.totalVolume()) * static_cast<long double>(G.totalVolume());

    auto accBeforeMove = intraClusterWeights_And_SumOfSquaredClusterVolumes(G, C);
    long double coverageBeforeMove = static_cast<long double>(accBeforeMove.first) / G.totalVolume();
    long double expectedCoverageBeforeMove = accBeforeMove.second / dTotalVolumeSquared;
    long double modBeforeMove = coverageBeforeMove - expectedCoverageBeforeMove;

    // apply move
    C[u] = to;
    clusterVolumes[to] += G.nodeVolume(u);
    clusterVolumes[from] -= G.nodeVolume(u);

    auto accAfterMove = intraClusterWeights_And_SumOfSquaredClusterVolumes(G, C);
    long double coverageAfterMove = static_cast<long double>(accAfterMove.first) / G.totalVolume();
    long double expectedCoverageAfterMove = accAfterMove.second / dTotalVolumeSquared;
    long double modAfterMove = coverageAfterMove - expectedCoverageAfterMove;

    bool comp = eq(modBeforeMove + adjustedGain, modAfterMove);
    ASSERT(comp,
           V(modBeforeMove + adjustedGain) << V(modAfterMove) << V(gain) << V(adjustedGain)
           << V(coverageBeforeMove) << V(expectedCoverageBeforeMove) << V(modBeforeMove)
           << V(coverageAfterMove) << V(expectedCoverageAfterMove) << V(modAfterMove));

    // revert move
    C[u] = from;
    clusterVolumes[to] -= G.nodeVolume(u);
    clusterVolumes[from] += G.nodeVolume(u);

    return comp;
#else
    (void)G;
    (void)C;
    (void)u;
    (void)to;
    (void)gain;
    (void)icw;
    return true;
#endif
  }

  static std::pair<ArcWeight, ArcWeight> intraClusterWeights_And_SumOfSquaredClusterVolumes(const Graph& G, const ds::Clustering& C) {
    ArcWeight intraClusterWeights = 0;
    ArcWeight sumOfSquaredClusterVolumes = 0;
    std::vector<ArcWeight> clusterVolumes(G.numNodes(), 0);

    for (NodeID u : G.nodes()) {
      ArcWeight arcVol = 0;
      for (const Arc& arc : G.arcsOf(u)) {
        if (C[u] == C[arc.head])
          intraClusterWeights += arc.weight;
        arcVol += arc.weight;
      }

      ArcWeight selfLoopWeight = G.nodeVolume(u) - arcVol;          // already accounted for as twice!
      ASSERT(selfLoopWeight >= 0.0);
      intraClusterWeights += selfLoopWeight;
      clusterVolumes[C[u]] += G.nodeVolume(u);
    }

    for (NodeID cluster : G.nodes()) {      // unused cluster IDs have volume 0
      sumOfSquaredClusterVolumes += clusterVolumes[cluster] * clusterVolumes[cluster];
    }
    return std::make_pair(intraClusterWeights, sumOfSquaredClusterVolumes);
  }

  long double doubleMod(const Graph& G, std::pair<int64_t, int64_t>& icwAndSoscv) {
    long double coverage = static_cast<long double>(icwAndSoscv.first) / static_cast<long double>(G.totalVolume());
    long double expectedCoverage = static_cast<long double>(icwAndSoscv.second) / (static_cast<long double>(G.totalVolume()) * static_cast<long double>(G.totalVolume()));
    return coverage - expectedCoverage;
  }

 private:
/*
    inline int64_t integerModGain(const ArcWeight weightFrom, const ArcWeight weightTo, const ArcWeight volFrom, const ArcWeight volTo, const ArcWeight volNode) {
        return 2 * (totalVolume * (weightTo - weightFrom) + resolutionGamma * volNode * (volFrom - volNode - volTo));
    }
*/

  inline double modularityGain(const ArcWeight weightFrom, const ArcWeight weightTo, const ArcWeight volFrom, const ArcWeight volTo, const ArcWeight volNode) {
    return weightTo - weightFrom + volNode * (volFrom - volNode - volTo) / totalVolume;
  }

  // gain computed by the above function would be adjusted by gain = gain / totalVolume
  inline long double adjustBasicModGain(double gain) {
    return 2.0L * gain / totalVolume;
  }

  // ~factor 3 on human_gene
  // multiplier = reciprocalTotalVolume * resolutionGamma * volNode
  inline double modularityGain(const ArcWeight weightTo, const ArcWeight volTo, const double multiplier) {
    return weightTo - multiplier * volTo;
  }

  // gain of the above function is adjusted by gain = gain / totalVolume after applying gain = gain (-weightFrom + resolutionGamma * reciprocalTotalVolume * volNode * (volFrom - volNode))
  inline long double adjustAdvancedModGain(double gain, const ArcWeight weightFrom, const ArcWeight volFrom, const ArcWeight volNode) {
    return 2.0L * reciprocalTotalVolume * (gain - weightFrom + reciprocalTotalVolume * volNode * (volFrom - volNode));
  }

  const Context& _context;
  ArcWeight totalVolume = 0;
  double reciprocalTotalVolume = 0.0;
  double volMultiplierDivByNodeVol = 0.0;
  std::vector<AtomicArcWeight> clusterVolumes;
  tbb::enumerable_thread_specific<IncidentClusterWeights> ets_incidentClusterWeights;
};
}
