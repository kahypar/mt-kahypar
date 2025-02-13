/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "deterministic_multilevel_coarsener.h"

#include <tbb/parallel_sort.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/hash.h"

namespace mt_kahypar {

template<typename TypeTraits>
bool DeterministicMultilevelCoarsener<TypeTraits>::coarseningPassImpl() {
  auto& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  const auto pass_start_time = std::chrono::high_resolution_clock::now();
  timer.start_timer("coarsening_pass", "Clustering");

  const Hypergraph& hg = Base::currentHypergraph();
  size_t num_nodes = Base::currentNumNodes();
  const double num_nodes_before_pass = num_nodes;
  vec<HypernodeID> clusters(num_nodes, kInvalidHypernode);
  tbb::parallel_for(UL(0), num_nodes, [&](HypernodeID u) {
    cluster_weight[u] = hg.nodeWeight(u);
    opportunistic_cluster_weight[u] = cluster_weight[u];
    propositions[u] = u;
    clusters[u] = u;
  });

  if (_context.coarsening.use_adaptive_edge_size) {
    hyperedge_size.resize(hg.initialNumEdges());
    tbb::parallel_for(HyperedgeID(0), hg.initialNumEdges(), [&](HyperedgeID e) {
      hyperedge_size[e] = hg.edgeSize(e);
    });
  }

  {
    const size_t contractable_nodes_per_subround = std::ceil(static_cast<double>(num_nodes - currentLevelContractionLimit()) / config.num_sub_rounds);
    permutation.random_grouping(num_nodes, _context.shared_memory.static_balancing_work_packages, config.prng());
    for (size_t sub_round = 0; sub_round < config.num_sub_rounds && num_nodes > currentLevelContractionLimit(); ++sub_round) {
      auto [first_bucket, last_bucket] = parallel::chunking::bounds(
        sub_round, config.num_buckets, config.num_buckets_per_sub_round);
      size_t first = permutation.bucket_bounds[first_bucket], last = permutation.bucket_bounds[last_bucket];
      // each vertex finds a cluster it wants to join
      tbb::parallel_for(first, last, [&](size_t pos) {
        const HypernodeID u = permutation.at(pos);
        if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {
          calculatePreferredTargetCluster(u, clusters);
        }
      });
      if (passed_nodes_from_previous_subround.size() > 0) {
        tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
          const HypernodeID u = passed_nodes_from_previous_subround[i];
          if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {
            calculatePreferredTargetCluster(u, clusters);
          }
        });
      }

      handleNodeSwaps(first, last, hg);

      if (_context.coarsening.split_contraction_limit_between_subrounds && _context.type == ContextType::main) {
        contractable_nodes.clear();
        for (size_t i = first; i < last; ++i) {
          HypernodeID u = permutation.at(i);
          HypernodeID target = propositions[u];
          if (target != u) {
            contractable_nodes.push_back(u);
          }
        }
        if (contractable_nodes.size() > contractable_nodes_per_subround) {
          std::shuffle(contractable_nodes.begin(), contractable_nodes.end(), std::mt19937(_context.partition.seed));
        }
        const size_t nodes_before_subround = num_nodes;
        size_t start = 0UL;
        size_t end = std::min(contractable_nodes.size(), contractable_nodes_per_subround);
        size_t actually_contracted_nodes = 0UL;
        while (end < contractable_nodes.size() && actually_contracted_nodes < contractable_nodes_per_subround) {
          tbb::enumerable_thread_specific<size_t> num_contracted_nodes{ 0 };
          if (start == 0UL) {
            tbb::parallel_for(end, contractable_nodes.size(), [&](const size_t i) {
              HypernodeID u = contractable_nodes[i];
              HypernodeID target = propositions[u];
              if (target != u) {
                __atomic_fetch_sub(&opportunistic_cluster_weight[target], hg.nodeWeight(u), __ATOMIC_RELAXED);
              }
            });
          } else {
            tbb::parallel_for(start, end, [&](const size_t i) {
              HypernodeID u = contractable_nodes[i];
              HypernodeID target = propositions[u];
              if (target != u) {
                __atomic_fetch_add(&opportunistic_cluster_weight[target], hg.nodeWeight(u), __ATOMIC_RELAXED);
              }
            });
          }

          tbb::parallel_for(start, end, [&](const size_t i) {
            HypernodeID u = contractable_nodes[i];
            HypernodeID target = propositions[u];
            if (target != u) {
              if (opportunistic_cluster_weight[target] <= maxAllowedNodeWeightInPass()) {
                // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
                if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
                  num_contracted_nodes.local() += 1;
                }
                clusters[u] = target;
                cluster_weight[target] = opportunistic_cluster_weight[target];
              } else {
                nodes_in_too_heavy_clusters.push_back_buffered(u);
              }
            }
          });
          num_nodes -= num_contracted_nodes.combine(std::plus<>());
          nodes_in_too_heavy_clusters.finalize();
          if (nodes_in_too_heavy_clusters.size() > 0) {
            handleNodesInTooHeavyClusters(num_nodes, clusters, hg);
            nodes_in_too_heavy_clusters.clear();
          }
          actually_contracted_nodes = nodes_before_subround - num_nodes;
          start = end;
          end = std::min(start + contractable_nodes_per_subround - actually_contracted_nodes, contractable_nodes.size());
          passed_nodes_from_previous_subround.clear();
        }

      } else {
        tbb::enumerable_thread_specific<size_t> num_contracted_nodes{ 0 };
        // already approve if we can grant all requests for proposed cluster
        // otherwise insert to shared vector so that we can group vertices by cluster
        tbb::parallel_for(first, last, [&](size_t pos) {
          HypernodeID u = permutation.at(pos);
          HypernodeID target = propositions[u];
          if (target != u) {
            if (opportunistic_cluster_weight[target] <= maxAllowedNodeWeightInPass()) {
              // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
              if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
                num_contracted_nodes.local() += 1;
              } else {
                cluster_weights_to_fix.push_back_buffered(u);
              }
              clusters[u] = target;
              cluster_weight[target] = opportunistic_cluster_weight[target];
            } else {
              if (opportunistic_cluster_weight[u] != hg.nodeWeight(u)) {
                // node u could still not move
                cluster_weights_to_fix.push_back_buffered(u);
              }
              nodes_in_too_heavy_clusters.push_back_buffered(u);
            }
          }
        });
        if (passed_nodes_from_previous_subround.size() > 0) {
          tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t pos) {
            const HypernodeID u = passed_nodes_from_previous_subround[pos];
            HypernodeID target = propositions[u];
            if (target != u) {
              if (opportunistic_cluster_weight[target] <= maxAllowedNodeWeightInPass()) {
                // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
                if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
                  num_contracted_nodes.local() += 1;
                } else {
                  cluster_weights_to_fix.push_back_buffered(u);
                }
                clusters[u] = target;
                cluster_weight[target] = opportunistic_cluster_weight[target];
              } else {
                if (opportunistic_cluster_weight[u] != hg.nodeWeight(u)) {
                  // node u could still not move
                  cluster_weights_to_fix.push_back_buffered(u);
                }
                nodes_in_too_heavy_clusters.push_back_buffered(u);
              }
            }
          });
        }
        const size_t contracted = num_contracted_nodes.combine(std::plus<>());
        DBG << "subround: " << sub_round << ", " << "contracted_nodes: " << contracted << "/" << contractable_nodes_per_subround;
        num_nodes -= contracted;
        nodes_in_too_heavy_clusters.finalize();
        if (nodes_in_too_heavy_clusters.size() > 0) {
          handleNodesInTooHeavyClusters(num_nodes, clusters, hg);
          nodes_in_too_heavy_clusters.clear();
        }

        // TODO: no need to fix cluster weights in the last subround
        cluster_weights_to_fix.finalize();
        if (cluster_weights_to_fix.size() > 0) {
          tbb::parallel_for(0UL, cluster_weights_to_fix.size(), [&](const size_t i) {
            const HypernodeID hn = cluster_weights_to_fix[i];
            const HypernodeID cluster = clusters[hn];
            if (cluster != hn) {
              cluster_weight[hn] -= hg.nodeWeight(hn);
              opportunistic_cluster_weight[hn] -= hg.nodeWeight(hn);
            }
          });
          cluster_weights_to_fix.clear();
        }

        if (_context.coarsening.use_adaptive_edge_size) {
          // update hyperedge sizes
          tbb::parallel_for(first, last, [&](size_t pos) {
            HypernodeID u = permutation.at(pos);
            if (u == propositions[u] || u == clusters[u]) {
              return;
            }

            // another idea to speed this up. this is slow if degree(clusters[u]) is unnecessarily large --> can mark smaller?
            // mark hg.incidentEdges(clusters[u]) in bitset
            // for each e in hg.incidentEdges(u)
            //     if e is marked --> reduce its size by 1
            // this is assuming that the vertex clusters[u] is still in that cluster. If it left in this subround, the number of edge size reductions is reduced by 1

            auto& ratings = default_rating_maps.local();
            for (HyperedgeID he : hg.incidentEdges(u)) {
              // this could be optimized to run once per affected hyperedge
              if (hg.edgeSize(he) >= _context.partition.ignore_hyperedge_size_threshold) continue;
              ratings.clear();
              for (HypernodeID v : hg.pins(he)) {
                ratings[clusters[v]] += 1;
              }
              hyperedge_size[he] = ratings.size();  // benign race
            }
          });
        }

        HEAVY_COARSENING_ASSERT([&] {
          vec<HypernodeWeight> cluster_weight_recalced(cluster_weight.size(), 0);
          for (const HypernodeID hn : hg.nodes()) {
            const HypernodeID cluster = clusters[hn];
            cluster_weight_recalced[cluster] += hg.nodeWeight(hn);
          }
          for (const HypernodeID c : hg.nodes()) {
            if (cluster_weight_recalced[c] > 0 && (cluster_weight_recalced[c] != cluster_weight[c] || cluster_weight_recalced[c] != opportunistic_cluster_weight[c])) {
              LOG << "ERROR wrong cluster_weight: " << V(sub_round) << ", " << V(c) << ", " << V(cluster_weight_recalced[c]) << ", " << V(cluster_weight[c]) << ", " << V(opportunistic_cluster_weight[c]) << std::endl;
            }
          }
        }(), "Clustering calculated wrong cluster-weights/opportunistic-cluster-weights");

        passed_nodes_from_previous_subround.clear();
      }
    }

  }
  timer.stop_timer("coarsening_pass");
  ++pass;
  if (num_nodes_before_pass / num_nodes <= _context.coarsening.minimum_shrink_factor) {
    return false;
  }
  _timer.start_timer("contraction", "Contraction");
  _uncoarseningData.performMultilevelContraction(std::move(clusters), true /* deterministic */, pass_start_time);
  _timer.stop_timer("contraction");
  passed_nodes_from_previous_subround.clear();
  return true;
}

template<typename TypeTraits>
void DeterministicMultilevelCoarsener<TypeTraits>::calculatePreferredTargetCluster(HypernodeID u, const vec<HypernodeID>& clusters) {
  const Hypergraph& hg = Base::currentHypergraph();
  auto& ratings = default_rating_maps.local();
  ratings.clear();

  // calculate ratings
  for (HyperedgeID he : hg.incidentEdges(u)) {
    HypernodeID he_size = hg.edgeSize(he);
    if (he_size < _context.partition.ignore_hyperedge_size_threshold) {
      // TODO should he_size filter use the original edge size or the adaptive one?
      he_size = _context.coarsening.use_adaptive_edge_size ? hyperedge_size[he] : he_size;
      double he_score = static_cast<double>(hg.edgeWeight(he)) / he_size;
      for (HypernodeID v : hg.pins(he)) {
        ratings[clusters[v]] += he_score;
      }
    }
  }

  // find highest rated, feasible cluster
  const PartitionID comm_u = hg.communityID(u);
  const HypernodeWeight weight_u = hg.nodeWeight(u);
  vec<HypernodeID>& best_targets = ties.local();
  double best_score = 0.0;

  for (const auto& entry : ratings) {
    HypernodeID target_cluster = entry.key;
    double target_score = entry.value;
    if (target_score >= best_score && target_cluster != u && hg.communityID(target_cluster) == comm_u
      && cluster_weight[target_cluster] + weight_u <= maxAllowedNodeWeightInPass()) {
      if (target_score > best_score) {
        best_targets.clear();
        best_score = target_score;
      }
      best_targets.push_back(target_cluster);
    }
  }

  HypernodeID best_target;
  if (best_targets.size() == 1) {
    best_target = best_targets[0];
  } else if (best_targets.empty()) {
    best_target = u;
  } else {
    // TODO: "rebase conflict"
    size_t pos = cluster_tie_breaker->select(best_targets.size() - 1, u);
    assert(pos < best_targets.size());
    best_target = best_targets[pos];
  }
  best_targets.clear();

  if (best_target != u) {
    propositions[u] = best_target;
    __atomic_fetch_add(&opportunistic_cluster_weight[best_target], hg.nodeWeight(u), __ATOMIC_RELAXED);
  }
}

template<typename TypeTraits>
size_t DeterministicMultilevelCoarsener<TypeTraits>::approveVerticesInTooHeavyClusters(vec<HypernodeID>& clusters) {
  const Hypergraph& hg = Base::currentHypergraph();
  tbb::enumerable_thread_specific<size_t> num_contracted_nodes { 0 };

  // group vertices by desired cluster, if their cluster is too heavy. approve the lower weight nodes first
  auto comp = [&](HypernodeID lhs, HypernodeID rhs) {
    HypernodeWeight wl = hg.nodeWeight(lhs), wr = hg.nodeWeight(rhs);
    return std::tie(propositions[lhs], wl, lhs) < std::tie(propositions[rhs], wr, rhs);
  };
  tbb::parallel_sort(nodes_in_too_heavy_clusters.begin(), nodes_in_too_heavy_clusters.end(), comp);

  tbb::parallel_for(UL(0), nodes_in_too_heavy_clusters.size(), [&](size_t pos) {
    HypernodeID target = propositions[nodes_in_too_heavy_clusters[pos]];
    // the first vertex for this cluster handles the approval
    size_t num_contracted_local = 0;
    if (pos == 0 || propositions[nodes_in_too_heavy_clusters[pos - 1]] != target) {
      HypernodeWeight target_weight = cluster_weight[target];
      size_t first_rejected = pos;
      // could be parallelized without extra memory but factor 2 work overhead and log(n) depth via binary search
      for (; ; ++first_rejected) {
        // we know that this cluster is too heavy, so the loop will terminate before
        assert(first_rejected < nodes_in_too_heavy_clusters.size());
        assert(propositions[nodes_in_too_heavy_clusters[first_rejected]] == target);
        HypernodeID v = nodes_in_too_heavy_clusters[first_rejected];
        if (target_weight + hg.nodeWeight(v) > maxAllowedNodeWeightInPass()) {
          break;
        }
        clusters[v] = target;
        target_weight += hg.nodeWeight(v);
        if (opportunistic_cluster_weight[v] == hg.nodeWeight(v)) {
          num_contracted_local += 1;
        }
      }
      cluster_weight[target] = target_weight;
      opportunistic_cluster_weight[target] = target_weight;
      num_contracted_nodes.local() += num_contracted_local;
    }
  });

  return num_contracted_nodes.combine(std::plus<>());
}

template<typename TypeTraits>
size_t DeterministicMultilevelCoarsener<TypeTraits>::recalculateForPassedOnHypernodes(vec<HypernodeID>& clusters) {
  size_t num_nodes = 0;
  const Hypergraph& hg = Base::currentHypergraph();
  tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
    const HypernodeID u = passed_nodes_from_previous_subround[i];
    if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {
      calculatePreferredTargetCluster(u, clusters);
    }
  });
  switch (_context.coarsening.swapStrategy) {
  case SwapResolutionStrategy::stay:
    tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
      const HypernodeID u = passed_nodes_from_previous_subround[i];
      const HypernodeID cluster_u = propositions[u];
      const HypernodeID cluster_v = propositions[cluster_u];
      if (u < cluster_u && u == cluster_v) {
        propositions[u] = u;
        propositions[cluster_u] = cluster_u;
        opportunistic_cluster_weight[cluster_u] -= hg.nodeWeight(u);
        opportunistic_cluster_weight[u] -= hg.nodeWeight(cluster_u);
      }
    });

    break;
  case SwapResolutionStrategy::to_smaller:
    tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
      const HypernodeID u = passed_nodes_from_previous_subround[i];
      const HypernodeID cluster_u = propositions[u];
      const HypernodeID cluster_v = propositions[cluster_u];
      if (u < cluster_u && u == cluster_v) {
        const HypernodeID target = opportunistic_cluster_weight[u] < opportunistic_cluster_weight[cluster_u] ? u : cluster_u;
        const HypernodeID source = target == u ? cluster_u : u;
        propositions[u] = target;
        propositions[cluster_u] = target;
        opportunistic_cluster_weight[source] -= hg.nodeWeight(target);
      }
    });

    break;
  case SwapResolutionStrategy::to_larger:
    tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
      const HypernodeID u = passed_nodes_from_previous_subround[i];
      const HypernodeID cluster_u = propositions[u];
      const HypernodeID cluster_v = propositions[cluster_u];
      if (u < cluster_u&& u == cluster_v) {
        const HypernodeID target = opportunistic_cluster_weight[u] > opportunistic_cluster_weight[cluster_u] ? u : cluster_u;
        const HypernodeID source = target == u ? cluster_u : u;
        propositions[u] = target;
        propositions[cluster_u] = target;
        opportunistic_cluster_weight[source] -= hg.nodeWeight(target);
      }
    });

    break;
  default:
    break;
  }

  tbb::enumerable_thread_specific<size_t> num_contracted_nodes{ 0 };

  tbb::parallel_for(0UL, passed_nodes_from_previous_subround.size(), [&](const size_t i) {
    //for (size_t i = 0; i < passed_nodes_from_previous_subround.size(); ++i) {
    const HypernodeID u = passed_nodes_from_previous_subround[i];
    HypernodeID target = propositions[u];
    if (target != u) {
      if (opportunistic_cluster_weight[target] <= maxAllowedNodeWeightInPass()) {
        // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
        if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
          num_contracted_nodes.local() += 1;
        }
        clusters[u] = target;
        cluster_weight[target] = opportunistic_cluster_weight[target];
      } else {
        nodes_in_too_heavy_clusters.push_back_buffered(u);
      }
    }
  });
  num_nodes += num_contracted_nodes.combine(std::plus<>());
  nodes_in_too_heavy_clusters.finalize();
  if (nodes_in_too_heavy_clusters.size() > 0) {
    num_nodes += approveVerticesInTooHeavyClusters(clusters);
  }
  return num_nodes;
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(DeterministicMultilevelCoarsener)

}
