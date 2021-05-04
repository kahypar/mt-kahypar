/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "deterministic_multilevel_coarsener.h"
#include "mt-kahypar/utils/progress_bar.h"

#include <tbb/parallel_sort.h>


namespace mt_kahypar {

void DeterministicMultilevelCoarsener::coarsenImpl() {
  HypernodeID initial_num_nodes = currentNumNodes();
  utils::ProgressBar progress_bar(initial_num_nodes, 0,
                                  _context.partition.verbose_output && _context.partition.enable_progress_bar);

  std::mt19937 prng(_context.partition.seed);
  size_t pass = 0;
  while (currentNumNodes() > _context.coarsening.contraction_limit) {
    auto pass_start_time = std::chrono::high_resolution_clock::now();
    const Hypergraph& hg = currentHypergraph();


    size_t num_nodes = currentNumNodes();
    double num_nodes_before_pass = num_nodes;
    vec<HypernodeID> clusters(num_nodes, kInvalidHypernode);

    vec<HypernodeID> first_solution;

    size_t num_reps = 5;
    for (size_t i = 0; i < num_reps; ++i) {
      num_nodes = currentNumNodes();
      num_nodes_before_pass = num_nodes;
      clusters = vec<HypernodeID>(num_nodes, kInvalidHypernode);
      tbb::parallel_for(0UL, num_nodes, [&](HypernodeID u) {
        cluster_weight[u] = hg.nodeWeight(u);
        opportunistic_cluster_weight[u] = cluster_weight[u];
        propositions[u] = u;
        clusters[u] = u;
      });

      permutation.random_grouping(num_nodes, _context.shared_memory.num_threads, prng());
      size_t num_sub_rounds = 16;
      size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
      size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);
      for (size_t sub_round = 0; sub_round < num_sub_rounds && num_nodes > currentLevelContractionLimit(); ++sub_round) {
        auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
        assert(first_bucket < last_bucket && last_bucket < permutation.bucket_bounds.size());
        size_t first = permutation.bucket_bounds[first_bucket];
        size_t last = permutation.bucket_bounds[last_bucket];

        // each vertex finds a cluster it wants to join
        tbb::parallel_for(first, last, [&](size_t pos) {
          assert(pos < num_nodes_before_pass);
          HypernodeID u = permutation.at(pos);
          assert(u < num_nodes_before_pass);
          if (cluster_weight[u] == hg.nodeWeight(u) && hg.nodeIsEnabled(u)) {  // u is a singleton
            calculatePreferredTargetCluster(permutation.at(pos), clusters);
          }
        });

        /*
        auto in_round = [&](HypernodeID v) {
          size_t bucket_v = permutation.get_bucket(v);
          return bucket_v >= first_bucket && bucket_v < last_bucket;
        };
        */
        // --> don't need two separate arrays 'clusters' and 'propositions'
        // if in_round(v), the entry is the proposition, otherwise its cluster

        tbb::enumerable_thread_specific<size_t> num_contracted_nodes {0};

        // already approve if we can grant all requests for proposed cluster
        // otherwise insert to shared vector so that we can group vertices by cluster
        tbb::parallel_for(first, last, [&](size_t pos) {
          HypernodeID u = permutation.at(pos);
          HypernodeID target = propositions[u];
          assert(target < num_nodes_before_pass);
          if (target != u) {
            if (opportunistic_cluster_weight[target] <= _context.coarsening.max_allowed_node_weight) {
              // if other nodes joined cluster u but u itself leaves for a different cluster, it doesn't count
              if (opportunistic_cluster_weight[u] == hg.nodeWeight(u)) {
                num_contracted_nodes.local() += 1;
              }
              clusters[u] = target;
              cluster_weight[target] = opportunistic_cluster_weight[target];
              // could subtract node weight from cluster_weight[u] to encourage more vertices to join in the second round
              // if the leader abandons its cluster. however, this is not so likely for the problematic cases
              // and introduces addtional complexities
            } else {
              nodes_in_too_heavy_clusters.push_back_buffered(u);
              // nodes_in_too_heavy_clusters.push_back_atomic(u);
            }
          }
        });

        nodes_in_too_heavy_clusters.finalize();
        if (nodes_in_too_heavy_clusters.size() > 0) {
          // group vertices by desired cluster, if their cluster is too heavy. approve the lower weight nodes first

          // if this is too slow, check out IPS4O as sorting algorithm, or packing the data into a struct?
          auto comp = [&](HypernodeID lhs, HypernodeID rhs) {
            HypernodeWeight wl = hg.nodeWeight(lhs), wr = hg.nodeWeight(rhs);
            return std::tie(propositions[lhs], wl, lhs) < std::tie(propositions[rhs], wr, rhs);
          };
          tbb::parallel_sort(nodes_in_too_heavy_clusters.begin(), nodes_in_too_heavy_clusters.end(), comp);

          tbb::parallel_for(0UL, nodes_in_too_heavy_clusters.size(), [&](size_t pos) {
            HypernodeID target = propositions[nodes_in_too_heavy_clusters[pos]];
            // the first vertex for this cluster handles the approval
            size_t num_contracted_local = 0;
            if (pos == 0 || propositions[nodes_in_too_heavy_clusters[pos - 1]] != target) {
              HypernodeWeight target_weight = cluster_weight[target];
              size_t first_rejected = pos;
              for (; ; ++first_rejected) {    // could be parallelized without extra memory but factor 2 work overhead and log(n) depth via binary search
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
        }
        nodes_in_too_heavy_clusters.clear();
        num_nodes -= num_contracted_nodes.combine(std::plus<>());
      }

      if (i == 0) {
        first_solution = clusters;
      } else {
        if (first_solution != clusters) {
          LOG << V(i) << "non-determinism in sync coarsening clustering";
        }
      }

    }

    ++pass;
    DBG << V(pass) << V(num_nodes) << V(num_nodes_before_pass);
    if (num_nodes_before_pass / num_nodes <= _context.coarsening.minimum_shrink_factor) {
      break;
    }
    performMultilevelContraction(std::move(clusters), pass_start_time);
    assert(num_nodes == currentNumNodes());
  }

  progress_bar += (initial_num_nodes - progress_bar.count());
  progress_bar.disable();
  finalize();
}

void DeterministicMultilevelCoarsener::calculatePreferredTargetCluster(HypernodeID u, const vec<HypernodeID>& clusters) {
  const Hypergraph& hg = currentHypergraph();
  auto& ratings = default_rating_maps.local();
  ratings.clear();

  // calculate ratings
  for (HyperedgeID he : hg.incidentEdges(u)) {
    HypernodeID he_size = hg.edgeSize(he);
    if (he_size < _context.partition.ignore_hyperedge_size_threshold) {
      double he_score = static_cast<double>(hg.edgeWeight(he)) / he_size;
      for (HypernodeID v : hg.pins(he)) {
        // TODO test impact!
        // the original code only counts the first occurence of each cluster in the same hyperedge
        // i.e., if multiple pins have the same cluster, only one counts
        // PaToH accounts for each occurrence
        // for n-level this obviously doesn't matter (where the code was taken from)
        ratings[clusters[v]] += he_score;
        // r[v] += he_score;    // avoid cluster lookup --> later aggregate in the rating eval loop
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

   /* // version that aggregates cluster scores post net/pin iteration
  for (auto it = ratings.begin(); it != ratings.end(); ++it) {
    HypernodeID neighbor = it->key;
    double neighbor_score = it->value;
    HypernodeID target_cluster = clusters[neighbor];
    if (target_cluster != neighbor) {     // this doesn't work with unstable leaders
      ratings[target_cluster] += neighbor_score;
    }
    double target_score = ratings[target_cluster];

    if (target_cluster != u && target_score >= best_score && hg.communityID(target_cluster) == comm_u
        && cluster_weight[target_cluster] + weight_u <= _context.coarsening.max_allowed_node_weight) {
      if (target_score > best_score) {
        best_targets.clear();
        best_score = target_score;
      }
      best_targets.push_back(target_cluster);
    }
  }
  */

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

}
