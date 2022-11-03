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

#include "local_moving_modularity.h"

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

#include "kahypar/utils/math.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <algorithm>

namespace mt_kahypar::metrics {
double modularity(const Graph& graph, const ds::Clustering& communities) {
  ASSERT(graph.canBeUsed());
  ASSERT(graph.numNodes() == communities.size());
  vec<NodeID> nodes(graph.numNodes());
  vec<double> cluster_mod(graph.numNodes(), 0.0);

  // make summation order deterministic!
  tbb::parallel_for(0UL, graph.numNodes(), [&](size_t pos) {
    nodes[pos] = pos;
  });
  tbb::parallel_sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs) {
    return std::tie(communities[lhs], lhs) < std::tie(communities[rhs], rhs);
  });

  // deterministic reduce doesn't have dynamic load balancing --> precompute the contributions and then sum them
  tbb::parallel_for(0UL, graph.numNodes(), [&](size_t pos) {
    NodeID x = nodes[pos];
    PartitionID comm = communities[x];
    double comm_vol = 0.0, internal = 0.0;
    if (pos == 0 || communities[nodes[pos - 1]] != comm) {
      for (size_t i = pos; i < nodes.size() && communities[nodes[i]] == comm; ++i) {
        NodeID u = nodes[i];
        comm_vol += graph.nodeVolume(u);
        for (const Arc& arc : graph.arcsOf(u)) {
          if (communities[arc.head] != comm) {
            internal -= arc.weight;
          }
        }
      }
      internal += comm_vol;
      cluster_mod[comm] = internal - (comm_vol * comm_vol) / graph.totalVolume();
    }
  });

  auto r = tbb::blocked_range<size_t>(0UL, graph.numNodes(), 1000);
  auto combine_range = [&](const tbb::blocked_range<size_t>& r, double partial) {
    return std::accumulate(cluster_mod.begin() + r.begin(), cluster_mod.begin() + r.end(), partial);
  };
  return tbb::parallel_deterministic_reduce(r, 0.0, combine_range, std::plus<>()) / graph.totalVolume();
}
}


namespace internal {
  template <typename T>
  struct Statistic {
    T min = 0;
    T q1 = 0;
    T med = 0;
    T q3 = 0;
    T max = 0;
    double avg = 0.0;
    double sd = 0.0;
  };

  template <typename T>
  Statistic<T> createStats(const std::vector<T>& vec, const double avg, const double stdev) {
    Statistic<T> stats;
    if (!vec.empty()) {
      const auto quartiles = kahypar::math::firstAndThirdQuartile(vec);
      stats.min = vec.front();
      stats.q1 = quartiles.first;
      stats.med = kahypar::math::median(vec);
      stats.q3 = quartiles.second;
      stats.max = vec.back();
      stats.avg = avg;
      stats.sd = stdev;
    }
    return stats;
  }
} // namespace internal

namespace mt_kahypar::community_detection {

bool ParallelLocalMovingModularity::localMoving(Graph& graph, ds::Clustering& communities, bool top_level) {
  ASSERT(graph.canBeUsed());
  _max_degree = graph.max_degree();
  _reciprocal_total_volume = 1.0 / graph.totalVolume();
  _vol_multiplier_div_by_node_vol = _reciprocal_total_volume;

  // init
  if (_context.partition.deterministic) {
    tbb::parallel_for(0UL, graph.numNodes(), [&](NodeID u) {
      communities[u] = u;
      _cluster_volumes[u].store(graph.nodeVolume(u), std::memory_order_relaxed);
      _cluster_weights[u].store(graph.nodeWeight(u), std::memory_order_relaxed);
    });
  } else {
    auto& nodes = permutation.permutation;
    nodes.resize(graph.numNodes());
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      nodes[u] = u;
      communities[u] = u;
      _cluster_volumes[u].store(graph.nodeVolume(u), std::memory_order_relaxed);
      _cluster_weights[u].store(graph.nodeWeight(u), std::memory_order_relaxed);
    });
  }

  DBG << "Louvain level" << V(graph.numNodes()) << V(graph.numArcs());

  // local moving
  bool clustering_changed = false;
  if ( graph.numArcs() > 0 ) {
    size_t number_of_nodes_moved = graph.numNodes();
    for (size_t round = 0;
        number_of_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * graph.numNodes()
        && round < _context.preprocessing.community_detection.max_pass_iterations; round++) {
      if (_context.partition.deterministic) {
        number_of_nodes_moved = synchronousParallelRound(graph, communities);
      } else {
        number_of_nodes_moved = parallelNonDeterministicRound(graph, communities);
      }
      clustering_changed |= number_of_nodes_moved > 0;
      DBG << "Louvain-Pass #" << round << " - num moves " << number_of_nodes_moved << " - Modularity:" << metrics::modularity(graph, communities);
    }
  }

  if (clustering_changed && _context.preprocessing.community_detection.use_isolated_nodes_treshold) {
    ASSERT(!_context.partition.deterministic);
    auto& nodes = permutation.permutation;
    std::vector<double> inv_gains;
    inv_gains.assign(nodes.size(), 0);
    tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
      double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
      inv_gains[i] = 1 / std::max(gain_to_iso, 1.0);
    });


    const double avg_inv_weigh_gain = utils::parallel_weighted_avg(inv_gains, graph.totalVolume(),
      [&](size_t i) {
        return graph.nodeVolume(nodes[i]);
      });
    const double stdev_inv_weigh_gain = utils::parallel_weighted_stdev(inv_gains, avg_inv_weigh_gain, graph.totalVolume(),
      [&](size_t i) {
        return graph.nodeVolume(nodes[i]);
      });

    const double stdev_factor = _context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor;
    HypernodeID first_separated = kInvalidHypernode;
    if (top_level && _context.preprocessing.community_detection.single_community_of_separated) {
      // we need to find one representative of the community of separated nodes
      for (size_t i = 0; i < nodes.size(); ++i) {
        const double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
        if (gain_to_iso <= 0 || gain_to_iso < 1 / (avg_inv_weigh_gain + stdev_factor * stdev_inv_weigh_gain)) {
          first_separated = nodes[i];
          break;
        }
      }
    }

    auto isolateNode = [&](const NodeID u) {
      const ArcWeight volU = graph.nodeVolume(u);
      const PartitionID from = communities[u];
      _cluster_volumes[from] -= volU;
      _cluster_weights[from] -= graph.nodeWeight(u);

      if (top_level && _context.preprocessing.community_detection.single_community_of_separated) {
        _cluster_volumes[first_separated] += volU;
        _cluster_weights[first_separated] += graph.nodeWeight(u);
        communities[u] = first_separated;
      } else {
        _cluster_volumes[u] += volU;
        _cluster_weights[u] += graph.nodeWeight(u);
        communities[u] = u;
      }
    };

    tbb::enumerable_thread_specific<int64_t> num_separated_local(0);

    tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
      const double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
      if (gain_to_iso <= 0 || gain_to_iso < 1 / (avg_inv_weigh_gain + stdev_factor * stdev_inv_weigh_gain)) {
        isolateNode(nodes[i]);
        graph.setIsolated(nodes[i]);
        num_separated_local.local()++;
      }
    });

    if (top_level) {
      tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
        if (!graph.isIsolated(nodes[i])) {
          const double comm_weight = _cluster_weights[communities[nodes[i]]].load(std::memory_order_relaxed);
          const double node_weight = graph.nodeWeight(nodes[i]);
          if (comm_weight == node_weight) {
            isolateNode(nodes[i]);
            graph.setIsolated(nodes[i]);
            num_separated_local.local()++;
          }
        }
      });
    }
    utils::Stats::instance().add_stat("num_separated_in_community_detection", num_separated_local.combine(std::plus<>()));


    if (top_level && _context.preprocessing.community_detection.collect_component_stats) {
      // some general stats
      ArcWeight adjacent_edges = 0;
      ArcWeight intern_edges = 0;
      ArcWeight adjacent_edges_non_iso = 0;
      HypernodeID truly_isolated_nodes = 0;
      for (const HypernodeID& node: nodes) {
        bool truly_isolated = true;
        ArcWeight current_adj = 0;
        if (graph.isIsolated(node)) {
          for (const Arc& arc : graph.arcsOf(node)) {
            truly_isolated &= !graph.isIsolated(arc.head);
            if (node < arc.head || !graph.isIsolated(arc.head)) {
              adjacent_edges += arc.weight;
              current_adj += arc.weight;
              if (graph.isIsolated(arc.head)) {
                intern_edges += arc.weight;
              }
            }
          }
          if (truly_isolated) {
            ++truly_isolated_nodes;
          } else {
            adjacent_edges_non_iso += current_adj;
          }
        }
      }
      utils::Stats::instance().add_stat("sep_total_edges", adjacent_edges);
      utils::Stats::instance().add_stat("sep_intern_edge_fraction", static_cast<double>(intern_edges) / static_cast<double>(adjacent_edges));
      utils::Stats::instance().add_stat<int32_t>("sep_single_nodes", truly_isolated_nodes);
      utils::Stats::instance().add_stat("sep_comp_intern_edge_fraction", static_cast<double>(intern_edges) / static_cast<double>(adjacent_edges_non_iso));

      // component analysis
      std::vector<bool> visited;
      visited.resize(graph.numNodes(), false);
      std::vector<HypernodeID> comp_sizes;
      std::vector<ArcWeight> comp_adjacent_edges;
      std::vector<double> comp_intern_edge_fracs;

      std::vector<HypernodeID> queue;
      for (const HypernodeID& node: nodes) {
        if (graph.isIsolated(node) && !visited.at(node)) {
          HypernodeID size = 0;
          ArcWeight adjacent = 0;
          ArcWeight intern = 0;
          queue.clear();
          queue.push_back(node);
          while (!queue.empty()) {
            const HypernodeID curr = queue.back();
            queue.pop_back();
            if (!visited[curr]) {
              visited[curr] = true;
              for (const Arc& arc : graph.arcsOf(curr)) {
                if (!visited[arc.head]) {
                  adjacent += arc.weight;
                  if (graph.isIsolated(arc.head)) {
                    intern += arc.weight;
                    queue.push_back(arc.head);
                  }
                }
              }
              ++size;
            }
          }

          ASSERT(size > 0);
          if (size > 1) {
            comp_sizes.push_back(size);
            comp_adjacent_edges.push_back(adjacent);
            comp_intern_edge_fracs.push_back(static_cast<double>(intern) / static_cast<double>(adjacent));
          }
        }
      }

      tbb::parallel_invoke([&] {
        std::sort(comp_sizes.begin(), comp_sizes.end());
      }, [&] {
        std::sort(comp_adjacent_edges.begin(), comp_adjacent_edges.end());
      }, [&] {
        std::sort(comp_intern_edge_fracs.begin(), comp_intern_edge_fracs.end());
      });
      double size_avg = utils::parallel_avg(comp_sizes, comp_sizes.size());
      auto size_stats = internal::createStats(comp_sizes, size_avg,
                                              utils::parallel_stdev(comp_sizes, size_avg, comp_sizes.size()));
      double adj_avg = utils::parallel_avg(comp_adjacent_edges, comp_adjacent_edges.size());
      auto adj_stats = internal::createStats(comp_adjacent_edges, adj_avg,
                                             utils::parallel_stdev(comp_adjacent_edges, adj_avg, comp_adjacent_edges.size()));
      double intern_avg = utils::parallel_avg(comp_intern_edge_fracs, comp_intern_edge_fracs.size());
      auto intern_stats = internal::createStats(comp_intern_edge_fracs, intern_avg,
                                                utils::parallel_stdev(comp_intern_edge_fracs, intern_avg, comp_intern_edge_fracs.size()));

      LOG << "name, total_num_node, sep_nodes, sep_single_nodes, sep_total_edges, sep_intern_edge_fraction, sep_comp_intern_edge_fraction"
             "num_components, comp_size_avg, comp_size_median, comp_size_sd, comp_size_max, "
             "adj_edges_avg, adj_edges_median, adj_edges_sd, adj_edges_max, "
             "intern_frac_avg, intern_frac_median, intern_frac_sd, intern_frac_max";
      const utils::Stats& stats = utils::Stats::instance();
      std::cout << _context.partition.graph_filename.substr(_context.partition.graph_filename.find_last_of('/') + 1) << ","
                << stats.get("total_num_nodes") << ","
                << stats.get("num_separated_in_community_detection") << ","
                << stats.get("sep_single_nodes") << ","
                << stats.get("sep_total_edges") << ","
                << stats.get("sep_intern_edge_fraction") << ","
                << stats.get("sep_comp_intern_edge_fraction") << ","
                << comp_sizes.size() << ","
                << size_stats.avg << "," << size_stats.med << "," << size_stats.sd << "," << size_stats.max << ","
                << adj_stats.avg << "," << adj_stats.med << "," << adj_stats.sd << "," << adj_stats.max << ","
                << intern_stats.avg << "," << intern_stats.med << "," << intern_stats.sd << "," << intern_stats.max
                << std::endl;
      std::exit(0);
    }

    if (top_level && !_context.preprocessing.community_detection.single_community_of_separated
        && _context.preprocessing.community_detection.separated_sub_communities) {
      size_t number_of_nodes_moved = graph.numNodes();
      for (size_t round = 0;
          number_of_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * graph.numNodes()
          && round < _context.preprocessing.community_detection.max_pass_iterations; round++) {
        number_of_nodes_moved = parallelNonDeterministicRound(graph, communities, true);
      }
    }
  }

  return clustering_changed;
}

size_t ParallelLocalMovingModularity::synchronousParallelRound(const Graph& graph, ds::Clustering& communities) {
  if (graph.numNodes() < 200) {
    return sequentialRound(graph, communities);
  }

  size_t seed = prng();
  permutation.random_grouping(graph.numNodes(), _context.shared_memory.static_balancing_work_packages, seed);
  size_t num_moved_nodes = 0;
  constexpr size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
  const size_t num_sub_rounds = _context.preprocessing.community_detection.num_sub_rounds_deterministic;
  size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);

  size_t max_round_size = 0;
  for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
    auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
    max_round_size = std::max(max_round_size,
                              size_t(permutation.bucket_bounds[last_bucket] - permutation.bucket_bounds[first_bucket]));
  }
  volume_updates_to.adapt_capacity(max_round_size);
  volume_updates_from.adapt_capacity(max_round_size);

  for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
    auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
    assert(first_bucket < last_bucket && last_bucket < permutation.bucket_bounds.size());
    size_t first = permutation.bucket_bounds[first_bucket];
    size_t last = permutation.bucket_bounds[last_bucket];

    tbb::enumerable_thread_specific<size_t> num_moved_local(0);
    tbb::parallel_for(first, last, [&](size_t pos) {
      HypernodeID u = permutation.at(pos);
      PartitionID best_cluster = computeMaxGainCluster(graph, communities, u, false);
      if (best_cluster != communities[u]) {
        volume_updates_from.push_back_buffered({ communities[u], u });
        volume_updates_to.push_back_buffered({best_cluster, u });
        num_moved_local.local() += 1;
      }
    });

    size_t num_moved_sub_round = num_moved_local.combine(std::plus<>());
    num_moved_nodes += num_moved_sub_round;

    // We can't do atomic adds of the volumes since they're not commutative and thus lead to non-deterministic decisions
    // Instead we sort the updates, and for each cluster let one thread sum up the updates.
    tbb::parallel_invoke([&] {
      volume_updates_to.finalize();
      tbb::parallel_sort(volume_updates_to.begin(), volume_updates_to.end());
    }, [&] {
      volume_updates_from.finalize();
      tbb::parallel_sort(volume_updates_from.begin(), volume_updates_from.end());
    });

    const size_t sz_to = volume_updates_to.size();
    tbb::parallel_for(0UL, sz_to, [&](size_t pos) {
      PartitionID c = volume_updates_to[pos].cluster;
      if (pos == 0 || volume_updates_to[pos - 1].cluster != c) {
        ArcWeight vol_delta = 0.0;
        for ( ; pos < sz_to && volume_updates_to[pos].cluster == c; ++pos) {
          vol_delta += graph.nodeVolume(volume_updates_to[pos].node);
          communities[volume_updates_to[pos].node] = c;
        }
        _cluster_volumes[c].store(_cluster_volumes[c].load(std::memory_order_relaxed) + vol_delta, std::memory_order_relaxed);
        // ...
      }
    });
    volume_updates_to.clear();

    const size_t sz_from = volume_updates_from.size();
    tbb::parallel_for(0UL, sz_from, [&](size_t pos) {
      PartitionID c = volume_updates_from[pos].cluster;
      if (pos == 0 || volume_updates_from[pos - 1].cluster != c) {
        ArcWeight vol_delta = 0.0;
        for ( ; pos < sz_from && volume_updates_from[pos].cluster == c; ++pos) {
          vol_delta -= graph.nodeVolume(volume_updates_from[pos].node);
        }
        _cluster_volumes[c].store(_cluster_volumes[c].load(std::memory_order_relaxed) + vol_delta, std::memory_order_relaxed);
        // ...
      }
    });
    volume_updates_from.clear();
  }

  return num_moved_nodes;
}

size_t ParallelLocalMovingModularity::sequentialRound(const Graph& graph, ds::Clustering& communities) {
  size_t seed = prng();
  permutation.sequential_fallback(graph.numNodes(), seed);
  size_t num_moved = 0;
  for (size_t i = 0; i < graph.numNodes(); ++i) {
    NodeID u = permutation.at(i);
    PartitionID best_cluster = computeMaxGainCluster(graph, communities, u, false);
    if (best_cluster != communities[u]) {
      _cluster_volumes[best_cluster] += graph.nodeVolume(u);
      _cluster_volumes[communities[u]] -= graph.nodeVolume(u);
      _cluster_weights[best_cluster] += graph.nodeWeight(u);
      _cluster_weights[communities[u]] -= graph.nodeWeight(u);
      communities[u] = best_cluster;
      num_moved++;
    }
  }
  return num_moved;
}

size_t ParallelLocalMovingModularity::parallelNonDeterministicRound(const Graph& graph, ds::Clustering& communities,
                                                                    const bool apply_to_isolated) {
  auto& nodes = permutation.permutation;
  if ( !_disable_randomization ) {
    utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
  }

  tbb::enumerable_thread_specific<size_t> local_number_of_nodes_moved(0);
  auto moveNode = [&](const NodeID u) {
    if (graph.isIsolated(u) == apply_to_isolated) {
      const ArcWeight volU = graph.nodeVolume(u);
      const PartitionID from = communities[u];
      PartitionID best_cluster = computeMaxGainCluster(graph, communities, u, apply_to_isolated);
      if (best_cluster != from) {
        _cluster_volumes[best_cluster] += volU;
        _cluster_volumes[from] -= volU;
        _cluster_weights[best_cluster] += graph.nodeWeight(u);
        _cluster_weights[from] -= graph.nodeWeight(u);
        communities[u] = best_cluster;
        ++local_number_of_nodes_moved.local();
      }
    }
  };

#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
  std::for_each(nodes.begin(), nodes.end(), moveNode);
#else
  tbb::parallel_for(0UL, nodes.size(), [&](size_t i) { moveNode(nodes[i]); });
#endif
  size_t number_of_nodes_moved = local_number_of_nodes_moved.combine(std::plus<>());
  return number_of_nodes_moved;
}


bool ParallelLocalMovingModularity::verifyGain(const Graph& graph, const ds::Clustering& communities, const NodeID u,
                                               const PartitionID to, double gain, double weight_from, double weight_to) {
  if (_context.partition.deterministic) {
    // the check is omitted, since changing the cluster volumes breaks determinism
    return true;
  }

  const PartitionID from = communities[u];

  long double adjustedGain = adjustAdvancedModGain(gain, weight_from, _cluster_volumes[from], graph.nodeVolume(u));
  const double volMultiplier = _vol_multiplier_div_by_node_vol * graph.nodeVolume(u);
  double modGain = modularityGain(weight_to, _cluster_volumes[to], volMultiplier);
  long double adjustedGainRecomputed = adjustAdvancedModGain(modGain, weight_from, _cluster_volumes[from], graph.nodeVolume(u));
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
  ds::Clustering communities_after_move = communities;
  communities_after_move[u] = to;
  _cluster_volumes[to] += graph.nodeVolume(u);
  _cluster_volumes[from] -= graph.nodeVolume(u);
  // ...

  auto accAfterMove = intraClusterWeightsAndSumOfSquaredClusterVolumes(graph, communities_after_move);
  long double coverageAfterMove = static_cast<long double>(accAfterMove.first) / graph.totalVolume();
  long double expectedCoverageAfterMove = accAfterMove.second / dTotalVolumeSquared;
  long double modAfterMove = coverageAfterMove - expectedCoverageAfterMove;

  const bool result = math::are_almost_equal_ld(modBeforeMove + adjustedGain, modAfterMove, 1e-8);
  ASSERT(result,
         V(modBeforeMove + adjustedGain) << V(modAfterMove) << V(gain) << V(adjustedGain)
                                         << V(coverageBeforeMove) << V(expectedCoverageBeforeMove) << V(modBeforeMove)
                                         << V(coverageAfterMove) << V(expectedCoverageAfterMove) << V(modAfterMove));

  _cluster_volumes[to] -= graph.nodeVolume(u);
  _cluster_volumes[from] += graph.nodeVolume(u);
  // ...

  return result;
}

std::pair<ArcWeight, ArcWeight> ParallelLocalMovingModularity::intraClusterWeightsAndSumOfSquaredClusterVolumes(
        const Graph& graph, const ds::Clustering& communities) {
  ArcWeight intraClusterWeights = 0;
  ArcWeight sumOfSquaredClusterVolumes = 0;
  vec<ArcWeight> cluster_volumes(graph.numNodes(), 0);

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
    cluster_volumes[communities[u]] += graph.nodeVolume(u);
  }

  for (NodeID cluster : graph.nodes()) {      // unused cluster IDs have volume 0
    sumOfSquaredClusterVolumes += cluster_volumes[cluster] * cluster_volumes[cluster];
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
/*
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
*/
}


}