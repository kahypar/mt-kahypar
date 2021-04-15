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

namespace mt_kahypar {

void DeterministicMultilevelCoarsener::coarsenImpl() {
  HypernodeID initial_num_nodes = currentNumNodes();
  utils::ProgressBar progress_bar(initial_num_nodes, 0,
                                  _context.partition.verbose_output && _context.partition.enable_progress_bar);

  std::mt19937 prng(_context.partition.seed);
  while (currentNumNodes() > _context.coarsening.contraction_limit) {
    auto pass_start_time = std::chrono::high_resolution_clock::now();
    const Hypergraph& hg = currentHypergraph();
    size_t num_nodes = currentNumNodes();
    double num_nodes_before_pass = num_nodes;
    vec<HypernodeID> clusters(num_nodes, kInvalidHypernode);
    tbb::parallel_for(0UL, num_nodes, [&](HypernodeID u) {
      cluster_weight[u] = hg.nodeWeight(u);
      clusters[u] = u;
    });

    permutation.random_grouping(num_nodes, _context.shared_memory.num_threads, prng());

    size_t num_sub_rounds = 16;
    size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
    size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);
    for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
      if (num_nodes < currentLevelContractionLimit()) {
        break;
      }

      auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
      size_t first = permutation.bucket_bounds[first_bucket];
      size_t last = permutation.bucket_bounds[last_bucket];

      // each vertex finds a cluster it wants to join
      tbb::parallel_for(first, last, [&](size_t pos) {
        assert(pos < permutation.permutation.size());
        HypernodeID u = permutation.at(pos);
        if (cluster_weight[u] == hg.nodeWeight(u)) {  // u is a singleton
          calculatePreferredTargetCluster(permutation.at(pos), clusters);
        }
      });

      // can ask whether vertex v is in current sub_round via
      // size_t bucket_v = permutation.get_bucket(v);
      // bool in_sub_round = bucket_v >= first_bucket && bucket_v < last_bucket;
      // --> wouldn't need two separate arrays 'clusters' and 'propositions'
      // if the vertex is in a current bucket, the entry is the proposition, otherwise its cluster

      // approve and execute or deny their requests
      tbb::parallel_for(first, last, [&](size_t pos) {
        HypernodeID u = permutation.at(pos);
      });

    }

    if (num_nodes_before_pass / num_nodes <= _context.coarsening.minimum_shrink_factor) {
      break;
    }
    performMultilevelContraction(std::move(clusters), pass_start_time);
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
    if (target_cluster != neighbor) {
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

  if (best_targets.size() == 1) {
    proposition[u] = best_targets[0];
  } else if (best_targets.size() == 0) {
    proposition[u] = u;
  } else {
    hashing::SimpleIntHash<uint32_t> sih;
    hashing::HashRNG hash_prng(sih, u);
    size_t pos = std::uniform_int_distribution<uint32_t>(0, best_targets.size() - 1)(hash_prng);
    proposition[u] = best_targets[pos];
  }
  best_targets.clear();
}

}
