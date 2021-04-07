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

#include "deterministic_label_propagation.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/chunking.h"
#include "mt-kahypar/parallel/parallel_counting_sort.h"

#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>

namespace mt_kahypar {

  bool DeterministicLabelPropagationRefiner::refineImpl(PartitionedHypergraph& phg,
                                                        const vec<HypernodeID>& ,
                                                        kahypar::Metrics& best_metrics,
                                                        const double)  {
    Gain overall_improvement = 0;
    size_t num_sub_rounds = context.refinement.deterministic_refinement.num_sub_rounds_sync_lp;

    for (size_t iter = 0; iter < context.refinement.label_propagation.maximum_iterations; ++iter) {
      size_t num_moves = 0;
      Gain round_improvement = 0;
      size_t n;
      if (context.refinement.deterministic_refinement.feistel_shuffling) {
        // get rid of this constant after initial tests, and use std::array in FeistelPermutation to store the keys
        constexpr size_t num_feistel_rounds = 4;
        feistel_permutation.create_permutation(num_feistel_rounds, phg.initialNumNodes(), prng);
        n = feistel_permutation.max_num_entries();
      } else {
        n = phg.initialNumNodes();
        // permutation.create_integer_permutation(n, context.shared_memory.num_threads, prng);
      }

      size_t sub_round_size = parallel::chunking::idiv_ceil(n, num_sub_rounds);
      for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
        // calculate moves
        moves_back.store(0, std::memory_order_relaxed);
        auto [first, last] = parallel::chunking::bounds(sub_round, n, sub_round_size);
        if (context.refinement.deterministic_refinement.feistel_shuffling) {
          tbb::parallel_for(HypernodeID(first), HypernodeID(last), [&](const HypernodeID cleartext) {
            const HypernodeID ciphertext = feistel_permutation.encrypt(cleartext);
            if (ciphertext < phg.initialNumNodes()) {
              calculateAndSaveBestMove(phg, ciphertext);
            }
          });
        } else {
          tbb::parallel_for(HypernodeID(first), HypernodeID(last), [&](const HypernodeID position) {
            // assert(position < permutation.permutation.size());
            // calculateAndSaveBestMove(phg, permutation.at(position));
            calculateAndSaveBestMove(phg, position);
          });
        }

        Gain sub_round_improvement = 0;
        size_t num_moves_in_sub_round = moves_back.load(std::memory_order_relaxed);
        if (num_moves_in_sub_round > 0) {
          // sync. then apply moves
          if (context.refinement.deterministic_refinement.apply_moves_by_maximal_prefix_in_block_pairs) {
            sub_round_improvement = applyMovesByMaximalPrefixesInBlockPairs(phg);
          } else {
            sub_round_improvement = applyMovesSortedByGainAndRevertUnbalanced(phg);
          }
        }
        round_improvement += sub_round_improvement;
        num_moves += num_moves_in_sub_round;
      }

      overall_improvement += round_improvement;
      DBG << V(iter) << V(num_moves) << V(round_improvement) << V(phg.initialNumNodes());

      if (num_moves == 0) {
        // no vertices with positive gain --> stop
        break;
      }
    }


    best_metrics.km1 -= overall_improvement;
    best_metrics.imbalance = metrics::imbalance(phg, context);

    return overall_improvement > 0;
  }

/*
 * for configs where we don't know exact gains --> have to trace the overall improvement with attributed gains
 * TODO active node set
 */
  Gain performMoveWithAttributedGain(PartitionedHypergraph& phg, const Move& m) {
    Gain attributed_gain = 0;
    auto objective_delta = [&](HyperedgeID he, HyperedgeWeight edge_weight, HypernodeID edge_size,
                               HypernodeID pin_count_in_from_part_after, HypernodeID pin_count_in_to_part_after) {
      attributed_gain -= km1Delta(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
    };
    phg.changeNodePart(m.node, m.from, m.to, objective_delta);
    return attributed_gain;
  }

  template<typename Predicate>
  Gain applyMovesIf(PartitionedHypergraph& phg, const vec<Move>& moves, size_t end, Predicate&& predicate) {
    auto range = tbb::blocked_range<size_t>(0UL, end);
    auto accum = [&](const tbb::blocked_range<size_t>& r, const Gain& init) -> Gain {
      Gain my_gain = init;
      for (size_t i = r.begin(); i < r.end(); ++i) {
        if (predicate(i)) {
          my_gain += performMoveWithAttributedGain(phg, moves[i]);
        }
      }
      return my_gain;
    };
    return tbb::parallel_reduce(range, 0, accum, std::plus<Gain>());
  }

  vec<HypernodeWeight> aggregatePartWeightDeltas(PartitionedHypergraph& phg, const vec<Move>& moves, size_t end) {
    // parallel reduce makes way too many vector copies
    tbb::enumerable_thread_specific<vec<HypernodeWeight>> ets_part_weight_diffs(phg.k(), 0);
    auto accum = [&](const tbb::blocked_range<size_t>& r) {
      auto& part_weights = ets_part_weight_diffs.local();
      for (size_t i = r.begin(); i < r.end(); ++i) {
        part_weights[moves[i].from] -= phg.nodeWeight(moves[i].node);
        part_weights[moves[i].to] += phg.nodeWeight(moves[i].node);
      }
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0UL, end), accum);
    vec<HypernodeWeight> res(phg.k(), 0);
    auto combine = [&](const vec<HypernodeWeight>& a) {
      for (size_t i = 0; i < res.size(); ++i) {
        res[i] += a[i];
      }
    };
    ets_part_weight_diffs.combine_each(combine);
    return res;
  }

  Gain DeterministicLabelPropagationRefiner::applyMovesSortedByGainAndRevertUnbalanced(PartitionedHypergraph& phg) {
    size_t num_moves = moves_back.load(std::memory_order_relaxed);
    tbb::parallel_sort(moves.begin(), moves.begin() + num_moves, [](const Move& m1, const Move& m2) {
      return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
    });

    size_t num_overloaded_blocks = 0;
    vec<HypernodeWeight> part_weights = aggregatePartWeightDeltas(phg, moves, num_moves);
    for (PartitionID i = 0; i < phg.k(); ++i) {
      part_weights[i] += phg.partWeight(i);
      if (part_weights[i] > context.partition.max_part_weights[i]) {
        num_overloaded_blocks++;
      }
    }

    // DBG << V(num_overloaded_blocks);

    size_t num_reverted_moves = 0;
    size_t j = num_moves;
    while (num_overloaded_blocks > 0 && j > 0) {
      Move& m = moves[--j];
      if (part_weights[m.to] > context.partition.max_part_weights[m.to]
          && part_weights[m.from] + phg.nodeWeight(m.node) <= context.partition.max_part_weights[m.from]) {
        part_weights[m.to] -= phg.nodeWeight(m.node);
        part_weights[m.from] += phg.nodeWeight(m.node);
        m.invalidate();
        num_reverted_moves++;
        if (part_weights[m.to] <= context.partition.max_part_weights[m.to]) {
          num_overloaded_blocks--;
        }
      }
    }

    if (num_overloaded_blocks > 0) {
      DBG << "still overloaded" << V(phg.initialNumNodes()) << V(phg.initialNumPins()) << V(num_overloaded_blocks) << V(num_moves);
    }

    // DBG << V(num_reverted_moves);

    // apply all moves that were not invalidated
    auto is_valid = [&](size_t pos) { return moves[pos].isValid(); };
    Gain gain = applyMovesIf(phg, moves, num_moves, is_valid);
    if (gain < 0) {
      DBG << "Kommando zurück" << V(gain) << V(num_moves) << V(num_reverted_moves);
      tbb::parallel_for(0UL, num_moves, [&](size_t i) {
        if (moves[i].isValid()) {
          std::swap(moves[i].from, moves[i].to);
        }
      });
      gain += applyMovesIf(phg, moves, num_moves, is_valid);
      assert(gain == 0);
    }
    return gain;
  }

  Gain DeterministicLabelPropagationRefiner::applyMovesByMaximalPrefixesInBlockPairs(PartitionedHypergraph& phg) {
    PartitionID k = phg.k();
    PartitionID max_key = k*k;
    auto index = [&](PartitionID b1, PartitionID b2) { return b1 * k + b2; };
    auto get_key = [&](const Move& m) { return index(m.from, m.to); };
    struct MovesWrapper {
      const Move& operator[](size_t i) const { return moves[i]; }
      size_t size() const { return sz; }
      const vec<Move>& moves;
      const size_t sz = 0;
    };

    MovesWrapper moves_wrapper { moves, moves_back.load(std::memory_order_relaxed) };

    // aggregate moves by direction. not in-place because of counting sort.
    // but it gives us the positions of the buckets right away
    auto positions = parallel::counting_sort(moves_wrapper, sorted_moves, max_key, get_key,
                                             context.shared_memory.num_threads);

    auto has_moves = [&](PartitionID p1, PartitionID p2) {
      size_t direction = index(p1, p2);
      return positions[direction + 1] != positions[direction];
    };

    vec<std::pair<PartitionID, PartitionID>> relevant_block_pairs;
    vec<size_t> involvements(k, 0);
    for (PartitionID p1 = 0; p1 < k; ++p1) {
      for (PartitionID p2 = p1 + 1; p2 < k; ++p2) {
        if (has_moves(p1,p2) && has_moves(p2,p1)) { // both directions have moves
          relevant_block_pairs.emplace_back(p1, p2);
          involvements[p1]++;
          involvements[p2]++;
        }
      }
    }

    // swap_prefix[index(p1,p2)] stores the first position of moves to revert out of the sequence of moves from p1 to p2
    vec<size_t> swap_prefix(max_key);
    vec<int64_t> part_weight_deltas(k, 0);

    tbb::parallel_for(0UL, relevant_block_pairs.size(), [&](size_t bp_index) {
      // sort both directions by gain (alternative: gain / weight?)
      auto [p1, p2] = relevant_block_pairs[bp_index];
      auto comp = [&](const Move& m1, const Move& m2) {
        return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
      };
      const auto b = sorted_moves.begin();
      size_t  i = positions[index(p1, p2)], i_last = positions[index(p1, p2) + 1],
              j = positions[index(p2, p1)], j_last = positions[index(p2, p1) + 1];
      std::sort(b + i, b + i_last, comp);
      std::sort(b + j, b + j_last, comp);


      // get balanced swap prefix
      HypernodeWeight budget_p1 = context.partition.max_part_weights[p1] - phg.partWeight(p1),
                      budget_p2 = context.partition.max_part_weights[p2] - phg.partWeight(p2);
      HypernodeWeight slack_p1 = budget_p1 / involvements[p1],
                      slack_p2 = budget_p2 / involvements[p2];

      int64_t balance = 0;
      std::tuple<size_t, size_t, int64_t> best {0,0,0};

      /*
       * this can be parallelized as follows.
       * 1. prefix sums of node weights over both move sequences
       * 2. pick middle of larger sequence, binary search for its prefix sum in the smaller sequence
       * 3. search for prefixes independently in both halves, and pick the better one
       *
       * in most cases we're expecting to take roughly as many moves as the size of the shorter sequence, from each of the sequences.
       * can we bias the search towards that?
       */

      // gain > 0 first. alternate depending on balance
      while (i < i_last && sorted_moves[i].gain > 0 && j < j_last && sorted_moves[j].gain > 0) {
        if (balance < 0 || (balance == 0 && sorted_moves[i].gain > sorted_moves[j].gain)) {
          // perform next move from p1 to p2
          balance += phg.nodeWeight(sorted_moves[i++].node);
        } else {
          // perform next move from p2 to p1
          balance -= phg.nodeWeight(sorted_moves[j++].node);
        }

        if (-balance <= slack_p1 && balance <= slack_p2) {
          best = {i,j,balance};
        }
      }

      // if one sequence is depleted or gain == 0. only do rebalancing in the other direction
      if (j == j_last || sorted_moves[j].gain == 0) {
        while (i < i_last && balance <= slack_p2 && (balance < 0 || sorted_moves[i].gain > 0)) {
          balance+= phg.nodeWeight(sorted_moves[i++].node);
          if (-balance <= slack_p1 && balance <= slack_p2) {
            best = {i,j,balance};
          }
        }
      } else if (i == i_last || sorted_moves[i].gain == 0) {
        while (j < j_last && -balance <= slack_p1 && (balance > 0 || sorted_moves[j].gain > 0)) {
          balance -= phg.nodeWeight(sorted_moves[j++].node);
          if (-balance <= slack_p1 && balance <= slack_p2) {
            best = {i,j,balance};
          }
        }
      }

      //swap_prefixes[bp_index] = best;
      swap_prefix[index(p1,p2)] = std::get<0>(best);
      swap_prefix[index(p2,p1)] = std::get<1>(best);
      int64_t best_balance = std::get<2>(best);

      // balance < 0 --> p1 got more weight, balance > 0 --> p2 got more weight
      __atomic_fetch_add(&part_weight_deltas[p1], best_balance, __ATOMIC_RELAXED);
      __atomic_fetch_sub(&part_weight_deltas[p2], best_balance, __ATOMIC_RELAXED);
    });

    auto remaining_moves = [&](PartitionID p1, PartitionID p2) {
      size_t direction = index(p1,p2);
      return positions[direction + 1] - swap_prefix[direction];
    };

    // simple check in case the slacks were too restrictive
    for (PartitionID p1 = 0; p1< k; ++p1) {
      for (PartitionID p2 = 0; p2 < k; ++p2) {
        if (remaining_moves(p2, p1) > 0) {
          // walking from the swap prefixes again might be too expensive
          // especially if they're still zero
          // we could walk from the end of the sequences we did in the pair-wise checks
        }
      }
    }

    auto in_prefix = [&](size_t pos) {
      return pos < swap_prefix[index(sorted_moves[pos].from, sorted_moves[pos].to)];
    };
    return applyMovesIf(phg, sorted_moves, sorted_moves.size(), in_prefix);
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
    return positions;
  }


}