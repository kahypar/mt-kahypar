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
using ds::Array;

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
    tbb::enumerable_thread_specific<int64_t> num_separated_local(0);
    HypernodeID first_separated = kInvalidHypernode;

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
      graph.setIsolated(u);
      num_separated_local.local()++;
    };

    const double stdev_factor = _context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor;
    const double stdev_factor_min = _context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor_min;
    if (!_context.preprocessing.community_detection.isolated_nodes_local_threshold) {
      std::vector<std::pair<double, HypernodeWeight>> inv_gains;
      inv_gains.assign(nodes.size(), std::make_pair(0, 0));
      tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
        double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
        inv_gains[i] = std::make_pair(1 / std::max(gain_to_iso, 1.0), graph.nodeWeight(nodes[i]));
      });
      auto map_to_first = std::function([](const std::pair<double, HypernodeWeight>& pair) { return pair.first; });

      const double avg_inv_weigh_gain = utils::parallel_weighted_avg(inv_gains, graph.totalVolume(),
        [&](size_t i) {
          return graph.nodeVolume(nodes[i]);
        }, map_to_first);
      const double stdev_inv_weigh_gain = utils::parallel_weighted_stdev(inv_gains, avg_inv_weigh_gain, graph.totalVolume(),
        [&](size_t i) {
          return graph.nodeVolume(nodes[i]);
        }, map_to_first);

      double min_gain = 1 / (avg_inv_weigh_gain + stdev_factor * stdev_inv_weigh_gain);
      const double min_gain_limit = 1 / (avg_inv_weigh_gain + stdev_factor_min * stdev_inv_weigh_gain);
      if (top_level && _context.preprocessing.community_detection.adjust_sd_factor) {
        LOG << "ADJUST SD FACTOR - " << V(min_gain);
        tbb::parallel_sort(inv_gains.begin(), inv_gains.end(), [](auto l, auto r) { return l.first > r.first; });
        HypernodeWeight sum = 0;
        for (auto pair: inv_gains) {
          sum += pair.second;
          if (sum >= (1 - _context.preprocessing.community_detection.sd_factor_core_size_target) * graph.totalWeight()) {
            const double max_inv_gain = pair.first;
            min_gain = std::max(min_gain, 1 / max_inv_gain);
            min_gain = std::min(min_gain, min_gain_limit);
            LOG << V(min_gain) << V(1 / max_inv_gain) << V(min_gain_limit);
            break;
          }
        }
      }

      if (top_level && _context.preprocessing.community_detection.single_community_of_separated) {
        // we need to find one representative of the community of separated nodes
        for (size_t i = 0; i < nodes.size(); ++i) {
          const double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
          if (gain_to_iso <= 0 || gain_to_iso <= min_gain) {
            first_separated = nodes[i];
            break;
          }
        }
      }

      tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
        const double gain_to_iso = computeGainComparedToIsolated(graph, communities, nodes[i]);
        if (gain_to_iso <= 0 || gain_to_iso <= min_gain) {
          isolateNode(nodes[i]);
        }
      });
    } else {
      // Compute mapped community ids with a parallel prefix sum
      Array<NodeID> mapping;
      mapping.assign(nodes.size(), 0);
      tbb::parallel_for(0UL, nodes.size(), [&](const HypernodeID& node) {
        if (!graph.isIsolated(node)) {
          mapping[communities[node]] = 1;
        }
      });
      parallel::TBBPrefixSum<NodeID, Array> mapping_prefix_sum(mapping);
      tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, mapping.size()), mapping_prefix_sum);
      const NodeID max_community_id = mapping_prefix_sum.total_sum();

      auto community_id = [&] (const NodeID node) {
        return mapping_prefix_sum[communities[node]];
      };

      Array<parallel::IntegralAtomicWrapper<NodeID>> first_node_index;
      first_node_index.assign(max_community_id + 1, parallel::IntegralAtomicWrapper<NodeID>(0)); // sentry
      tbb::parallel_for(0UL, nodes.size(), [&](const HypernodeID& node) {
        if (!graph.isIsolated(node)) {
          first_node_index[community_id(node) + 1].fetch_add(1);
        }
      });
      parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<NodeID>, Array> first_node_index_prefix_sum(first_node_index);
      tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, first_node_index.size()), first_node_index_prefix_sum);

      Array<NodeID> nodes_per_community;
      nodes_per_community.assign(nodes.size(), 0);
      tbb::parallel_for(0UL, nodes.size(), [&](const NodeID& node) {
        if (!graph.isIsolated(node)) {
          const NodeID pos = first_node_index[community_id(node)].fetch_add(1);
          nodes_per_community[pos] = node;
        }
      });

      std::vector<double> inv_gains;
      for (NodeID c = 0; c < nodes.size(); ++c) {
        if (mapping_prefix_sum[c + 1] > mapping_prefix_sum[c]) {
          const NodeID community = mapping_prefix_sum[c];
          inv_gains.clear();
          const NodeID start_pos = first_node_index_prefix_sum[community];
          const NodeID end_pos = first_node_index_prefix_sum[community + 1];

          inv_gains.assign(end_pos - start_pos, 0);
          tbb::parallel_for(start_pos, end_pos, [&](NodeID i) {
            const NodeID node = nodes_per_community[i];
            ASSERT(community_id(node) == community);
            double gain_to_iso = computeGainComparedToIsolated(graph, communities, node);
            inv_gains[i - start_pos] = 1 / std::max(gain_to_iso, 1.0);
          });
          // LOG << V(inv_gains.size()) << V(inv_gains[0]) << V(_cluster_volumes[c]) << V(graph.nodeVolume(nodes_per_community[start_pos]));

          const double avg_inv_weigh_gain = utils::parallel_weighted_avg(inv_gains, _cluster_volumes[c],
            [&](size_t i) {
              return graph.nodeVolume(nodes_per_community[i + start_pos]);
            });
          const double stdev_inv_weigh_gain = utils::parallel_weighted_stdev(inv_gains, avg_inv_weigh_gain, _cluster_volumes[c],
            [&](size_t i) {
              return graph.nodeVolume(nodes_per_community[i + start_pos]);
            });

          tbb::parallel_for(start_pos, end_pos, [&](NodeID i) {
            const NodeID node = nodes_per_community[i];
            const double gain_to_iso = computeGainComparedToIsolated(graph, communities, node);
            if (gain_to_iso <= 0 || gain_to_iso < 1 / (avg_inv_weigh_gain + stdev_factor * stdev_inv_weigh_gain)) {
              isolateNode(node);
            }
          });
          // LOG << V(community) << V(avg_inv_weigh_gain) << V(stdev_inv_weigh_gain);
        }
      }
    }

    if (top_level) {
      tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
        if (!graph.isIsolated(nodes[i])) {
          const double comm_weight = _cluster_weights[communities[nodes[i]]].load(std::memory_order_relaxed);
          const double node_weight = graph.nodeWeight(nodes[i]);
          if (comm_weight == node_weight) {
            isolateNode(nodes[i]);
          }
        }
      });
    }
    utils::Stats::instance().add_stat("num_separated_in_community_detection", num_separated_local.combine(std::plus<>()));
    LOG << V(num_separated_local.combine(std::plus<>()));


    if (top_level && _context.preprocessing.community_detection.collect_component_stats) {
      // some general stats
      tbb::enumerable_thread_specific<ArcWeight> adjacent_edges_l(0);
      tbb::enumerable_thread_specific<ArcWeight> intern_edges_l(0);
      tbb::enumerable_thread_specific<ArcWeight> core_edges_l(0);
      tbb::enumerable_thread_specific<HypernodeWeight> sep_weight_l(0);
      tbb::enumerable_thread_specific<HypernodeWeight> core_weight_l(0);
      tbb::parallel_for(0UL, nodes.size(), [&](size_t i) {
        const HypernodeID node = nodes[i];
        if (graph.isIsolated(node)) {
          for (const Arc& arc : graph.arcsOf(node)) {
            adjacent_edges_l.local() += arc.weight;
            if (graph.isIsolated(arc.head)) {
              intern_edges_l.local() += arc.weight;
            }
          }
          sep_weight_l.local() += graph.nodeWeight(node);
        } else {
          for (const Arc& arc : graph.arcsOf(node)) {
            if (!graph.isIsolated(arc.head)) {
              core_edges_l.local() += arc.weight;
            } else {
              adjacent_edges_l.local() += arc.weight;
            }
          }
          core_weight_l.local() += graph.nodeWeight(node);
        }
      });

      const ArcWeight adjacent_edges = adjacent_edges_l.combine(std::plus<>());
      const ArcWeight intern_edges = intern_edges_l.combine(std::plus<>());
      const ArcWeight core_edges = core_edges_l.combine(std::plus<>());
      const HypernodeWeight sep_weight = sep_weight_l.combine(std::plus<>());
      const HypernodeWeight core_weight = core_weight_l.combine(std::plus<>());
      utils::Stats::instance().update_stat("sep_total_edges", adjacent_edges);
      utils::Stats::instance().update_stat("sep_intern_edge_fraction", static_cast<double>(intern_edges) / static_cast<double>(adjacent_edges));
      utils::Stats::instance().update_stat("total_density", graph.totalVolume() / static_cast<double>(graph.totalWeight()));
      utils::Stats::instance().update_stat("core_density", core_edges / static_cast<double>(core_weight));
      utils::Stats::instance().update_stat("sep_density", intern_edges / static_cast<double>(sep_weight));
      utils::Stats::instance().update_stat("core_weight", core_weight / static_cast<double>(graph.totalWeight()));
      utils::Stats::instance().update_stat("stdev_factor", _context.preprocessing.community_detection.isolated_nodes_threshold_stdev_factor);
    }

    if (top_level && !_context.preprocessing.community_detection.single_community_of_separated
        && _context.preprocessing.community_detection.separated_sub_communities) {
      size_t number_of_nodes_moved = graph.numNodes();
      for (size_t round = 0;
          number_of_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * graph.numNodes()
          && round < _context.preprocessing.community_detection.max_pass_iterations; round++) {
        number_of_nodes_moved = parallelNonDeterministicRound(graph, communities, true);
        LOG << V(round) << V(number_of_nodes_moved);
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