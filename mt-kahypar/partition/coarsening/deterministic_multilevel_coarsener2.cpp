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

#include "deterministic_multilevel_coarsener2.h"

#include <tbb/parallel_sort.h>

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

template<typename TypeTraits>
bool DeterministicMultilevelCoarsener2<TypeTraits>::coarseningPassImpl() {
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
  const size_t num_edges_before = hg.initialNumEdges();
  const size_t num_pins_before = hg.initialNumPins();
  parallel::scalable_vector<HypernodeID> permutation(num_nodes);
  tbb::parallel_for(0UL, num_nodes, [&](size_t i) {
    permutation[i] = i;
  });
  std::shuffle(permutation.begin(), permutation.end(), std::mt19937(_context.partition.seed));
  if constexpr (prefixDoubling) {
    size_t first = 0;
    size_t last = 1;
    while (first < permutation.size()) {
      DBG << V(first) << ", " << V(last);
      // each vertex finds a cluster it wants to join
      tbb::parallel_for(first, last, [&](size_t pos) {
        const HypernodeID u = permutation.at(pos);
        if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {
          calculatePreferredTargetCluster(u, clusters);
        }
      });
      handleNodeSwaps(permutation, first, last, hg);

      tbb::enumerable_thread_specific<size_t> num_contracted_nodes{ 0 };

      // already approve if we can grant all requests for proposed cluster
      // otherwise insert to shared vector so that we can group vertices by cluster
      tbb::parallel_for(first, last, [&](size_t pos) {
        HypernodeID u = permutation.at(pos);
        HypernodeID target = propositions[u];
        if (target != u) {
          if (opportunistic_cluster_weight[target] <= _context.coarsening.max_allowed_node_weight) {
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
        num_nodes -= approveVerticesInTooHeavyClusters(clusters);
      }

      nodes_in_too_heavy_clusters.clear();
      first = last;
      last = std::min(permutation.size(), last * 2);
    }
  } else {
    // one node at a time
    for (const HypernodeID u : permutation) {
      if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {
        calculatePreferredTargetCluster(u, clusters);
      }
      const HypernodeID target = propositions[u];
        if (target != u) {
          if (opportunistic_cluster_weight[target] <= _context.coarsening.max_allowed_node_weight) {
            // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
            if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
              num_nodes--;
            }
            clusters[u] = target;
            cluster_weight[target] = opportunistic_cluster_weight[target];
          } else {
            std::cout << "THIS SHOULD NOT HAPPEN" << std::endl;
          }
        }
    }
  }

  timer.stop_timer("coarsening_pass");
  ++pass;
  if (num_nodes_before_pass / num_nodes <= _context.coarsening.minimum_shrink_factor) {
    return false;
  }
  utils::Measurements& measurements = utils::Utilities::instance().getMeasurements(_context.utility_id);
  if (_context.type == ContextType::main) {
    std::unordered_map<HypernodeID, HypernodeWeight> cluster_sizes;
    for (const HypernodeID& hn : hg.nodes()) {
      const auto cluster_id = clusters[hn];
      const auto weight = hg.nodeWeight(hn);
      auto it = cluster_sizes.find(cluster_id);
      if (it == cluster_sizes.end()) {
        cluster_sizes.insert({ cluster_id, weight });
      } else {
        it->second += weight;
      }
    }
    parallel::scalable_vector<HypernodeWeight> weights(cluster_sizes.size());
    size_t index = 0;
    size_t sum = 0;
    for (const auto& e : cluster_sizes) {
      weights[index++] = e.second;
      sum += e.second;
    }
    std::sort(weights.begin(), weights.end());
    const size_t min = weights[0];
    const size_t max = weights[weights.size() - 1];
    const size_t median = weights[weights.size() / 2];
    const double avg = static_cast<double>(sum) / weights.size();
    const size_t count = weights.size();
    measurements.min_cluster_size.push_back(min);
    measurements.max_cluster_size.push_back(max);
    measurements.median_cluster_size.push_back(median);
    measurements.avg_cluster_size.push_back(avg);
    measurements.cluster_count.push_back(count);
    size_t singletons = 0;
    while (weights[singletons] == 1) {
      singletons++;
    }
    measurements.num_singletons.push_back(singletons);
  }
  _timer.start_timer("contraction", "Contraction");
  _uncoarseningData.performMultilevelContraction(std::move(clusters), true /* deterministic */, pass_start_time);
  _timer.stop_timer("contraction");
  if (_context.type == ContextType::main) {
    Hypergraph& after = Base::currentHypergraph();
    const size_t eliminatedEdges = num_edges_before - after.initialNumEdges();
    const size_t eliminatedPins = num_pins_before - after.initialNumPins();
    measurements.eliminated_edges.push_back(eliminatedEdges);
    measurements.eliminated_pins.push_back(eliminatedPins);
    size_t score = 0;
    for (auto edge : after.edges()) {
      const HyperedgeWeight weight = after.edgeWeight(edge);
      const HypernodeWeight size = after.edgeSize(edge);
      score += weight * size;
    }
    measurements.score.push_back(score);
  }
  return true;
}

template<typename TypeTraits>
void DeterministicMultilevelCoarsener2<TypeTraits>::calculatePreferredTargetCluster(HypernodeID u, const vec<HypernodeID>& clusters) {
  const Hypergraph& hg = Base::currentHypergraph();
  auto& ratings = default_rating_maps.local();
  ratings.clear();

  // calculate ratings
  for (HyperedgeID he : hg.incidentEdges(u)) {
    HypernodeID he_size = hg.edgeSize(he);
    if (he_size < _context.partition.ignore_hyperedge_size_threshold) {
      double he_score = static_cast<double>(hg.edgeWeight(he)) / (he_size - 1);
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
      && cluster_weight[target_cluster] + weight_u <= _context.coarsening.max_allowed_node_weight) {
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
    hashing::SimpleIntHash<uint32_t> sih;
    hashing::HashRNG hash_prng(sih, u);
    size_t pos = std::uniform_int_distribution<uint32_t>(0, best_targets.size() - 1)(hash_prng);
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
size_t DeterministicMultilevelCoarsener2<TypeTraits>::approveVerticesInTooHeavyClusters(vec<HypernodeID>& clusters) {
  const Hypergraph& hg = Base::currentHypergraph();
  tbb::enumerable_thread_specific<size_t> num_contracted_nodes{ 0 };

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
        if (target_weight + hg.nodeWeight(v) > _context.coarsening.max_allowed_node_weight) {
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

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(DeterministicMultilevelCoarsener2)

}
