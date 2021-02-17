/*******************************************************************************
 * This file is part of MT-KaHyPar.
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

#include <tbb/parallel_sort.h>
#include "deterministic_label_propagation.h"


#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/chunking.h"
#include "mt-kahypar/parallel/parallel_counting_sort.h"

namespace mt_kahypar {

  bool DeterministicLabelPropagationRefiner::refineImpl(PartitionedHypergraph& phg,
                                                        const vec<HypernodeID>& ,
                                                        kahypar::Metrics& best_metrics,
                                                        const double)  {
    Gain overall_improvement = 0;
    size_t num_sub_rounds = context.refinement.deterministic_refinement.num_sub_rounds_sync_lp;
    size_t sub_round_size = parallel::chunking::idiv_ceil(phg.initialNumNodes(), num_sub_rounds);

    for (size_t iter = 0; iter < context.refinement.label_propagation.maximum_iterations; ++iter) {

      if (context.refinement.deterministic_refinement.use_coloring) {
        // get vertex permutation from coloring
        coloring(phg);
      } else if (context.refinement.deterministic_refinement.feistel_shuffling) {
        // get rid of this constant after initial tests, and use std::array in FeistelPermutation to store the keys
        constexpr size_t num_feistel_rounds = 4;
        feistel_permutation.create_permutation(num_feistel_rounds, phg.initialNumNodes(), prng);
      } else {
        // adapt number of nodes to active node set
        permutation.create_integer_permutation(phg.initialNumNodes(), context.shared_memory.num_threads, prng);
      }

      for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
        // calculate moves
        auto [first, last] = parallel::chunking::bounds(sub_round, phg.initialNumNodes(), sub_round_size);
        tbb::parallel_for(HypernodeID(first), HypernodeID(last), [&](const HypernodeID position) {
          calculateAndSaveBestMove(phg, permutation.at(position));
        });

        // sync, then apply moves
        applyAllMoves(phg);

        // for now apply all moves. temporary!
        // realistic approaches: sort move direction queues, find an improvement for k-way
        // unrealistic approaches: linear program (Ugander), negative cycle detection (Chris)
      }


    }


    best_metrics.km1 -= overall_improvement;
    best_metrics.imbalance = metrics::imbalance(phg, context);
    return overall_improvement > 0;
  }

  void DeterministicLabelPropagationRefiner::applyAllMoves(PartitionedHypergraph& phg) {
    std::atomic<Gain> gain(0);
    tbb::parallel_for(0UL, moves_back.load(std::memory_order_relaxed), [&](const size_t move_pos) {
      gain.fetch_add(performMoveWithAttributedGain(phg, moves[move_pos]), std::memory_order_relaxed);
    });
    moves_back.store(0, std::memory_order_relaxed);
  }

  void DeterministicLabelPropagationRefiner::applyMovesByMaximalPrefixesInBlockPairs(PartitionedHypergraph& phg) {

    PartitionID k = context.partition.k;
    PartitionID max_key = k*k;
    auto index = [&](PartitionID b1, PartitionID b2) {
      return b1 * k + b2;
    };
    auto get_key = [&](const Move& m) {
      return index(m.from, m.to);
    };

    // aggregate moves by direction. not in-place because of counting sort.
    // but it gives us the positions of the buckets right away
    auto positions = parallel::counting_sort(moves, sorted_moves, max_key, get_key, context.shared_memory.num_threads);

    vec<std::pair<PartitionID, PartitionID>> relevant_block_pairs;
    for (PartitionID i = 0; i < k; ++i) {
      for (PartitionID j = i + 1; j < k; ++j) {
        if (positions[index(i,j) + 1] != positions[index(i,j)]
            && positions[index(j,i) + 1] != positions[index(j,i)]) { // neither direction (i,j) nor (j,i) empty
          relevant_block_pairs.emplace_back(j, i);
        }
      }
    }


    tbb::parallel_for(0UL, relevant_block_pairs.size(), [&](size_t i) {
      auto [b1, b2] = relevant_block_pairs[i];
      auto comp = [](const Move& m1, const Move& m2) { return m1.gain > m2.gain; };
      const auto b = sorted_moves.begin();

      // sort both directions by gain (alternative: gain / weight?)
      std::sort(b + positions[index(b1, b2)],
                b + positions[index(b1, b2) + 1],
                comp);

      std::sort(b + positions[index(b2, b1)],
                b + positions[index(b2, b1) + 1],
                comp);

      // get best prefix with feasible weight

    });

  }

  vec<size_t> DeterministicLabelPropagationRefiner::aggregateDirectionBucketsInplace() {
    // this can be done more efficiently, i.e. in linear time with counting sort. leave as is for now, for simplicity
    tbb::parallel_sort(moves.begin(),
                       moves.begin() + moves_back.load(std::memory_order_relaxed),
                       [](const Move& m1, const Move& m2) { return std::tie(m1.from, m1.to) < std::tie(m2.from, m2.to); }
    );

    PartitionID k = context.partition.k;
    auto index = [&](PartitionID b1, PartitionID b2) {
      return b1 * k + b2;
    };

    vec<size_t> positions(k*k + 1, 0);
    size_t pos = 0;
    for (PartitionID i = 0; i < k; ++i) {
      for (PartitionID j = 0; j < k; ++j) {
        while (pos < moves.size() && moves[pos].from == i && moves[pos].to == j) {
          ++pos;
        }
        positions[index(i,j) + 1] = pos;
      }
    }
  }

  size_t DeterministicLabelPropagationRefiner::coloring(PartitionedHypergraph& phg) {
    // idea ignore large hyperedges to get good colorings
    // then vertices of the same color can be moved independently
    // correct the error from large hyperedges later with gain recalc for these specific hyperedges
    // the error is likely small, since large hyperedges rarely admit improvement

    LOG << "not yet implemented";
    return 0;
  }

}