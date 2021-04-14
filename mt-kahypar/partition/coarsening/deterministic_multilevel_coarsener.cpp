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
    size_t num_nodes = currentNumNodes();
    vec<HypernodeID> cluster(num_nodes, kInvalidHypernode);
    permutation.random_grouping(num_nodes, _context.shared_memory.num_threads, prng());

    size_t num_sub_rounds = 16;
    size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
    size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);
    for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
      auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
      size_t first = permutation.bucket_bounds[first_bucket];
      size_t last = permutation.bucket_bounds[last_bucket];
      tbb::parallel_for(first, last, [&](size_t pos) {
        assert(pos < permutation.permutation.size());
        HypernodeID u = permutation.at(pos);

        // TODO find cluster

        // can ask whether vertex v is in current sub_round via
        // size_t bucket_v = permutation.get_bucket(v);
        // bool in_sub_round = bucket_v >= first_bucket && bucket_v < last_bucket;

        cluster[u] = u;
      });

    }


  }

  progress_bar += (initial_num_nodes - progress_bar.count());
  progress_bar.disable();
}

}
