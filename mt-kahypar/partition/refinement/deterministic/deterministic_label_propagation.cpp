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
    constexpr size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
    const size_t num_sub_rounds = context.refinement.deterministic_refinement.num_sub_rounds_sync_lp;
    const size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);

    const bool log = false && context.type == kahypar::ContextType::main;

    for (size_t iter = 0; iter < context.refinement.label_propagation.maximum_iterations; ++iter) {
      size_t num_moves = 0;
      Gain round_improvement = 0;
      size_t n = phg.initialNumNodes();
      permutation.random_grouping(n, context.shared_memory.static_balancing_work_packages, prng());
      for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
        auto [first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
        assert(first_bucket < last_bucket && last_bucket < permutation.bucket_bounds.size());
        size_t first = permutation.bucket_bounds[first_bucket], last = permutation.bucket_bounds[last_bucket];
        moves_back.store(0, std::memory_order_relaxed);

        // calculate moves
        auto t1 = tbb::tick_count::now();
        if (phg.k() == 2) {
          tbb::parallel_for(HypernodeID(first), HypernodeID(last), [&](const HypernodeID position) {
            assert(position < permutation.permutation.size());
            calculateAndSaveBestMoveTwoWay(phg, permutation.at(position));
          });
        } else {
          tbb::parallel_for(HypernodeID(first), HypernodeID(last), [&](const HypernodeID position) {
            assert(position < permutation.permutation.size());
            calculateAndSaveBestMove(phg, permutation.at(position));
          });
        }
        if (log) LOG << "calc moves time" << (tbb::tick_count::now() - t1).seconds();

        Gain sub_round_improvement = 0;
        size_t num_moves_in_sub_round = moves_back.load(std::memory_order_relaxed);
        DBG << V(num_moves_in_sub_round);
        if (num_moves_in_sub_round > 0) {
          auto t2 = tbb::tick_count::now();
          sub_round_improvement = applyMovesByMaximalPrefixesInBlockPairs(phg);
          auto t3 = tbb::tick_count::now();
          if (sub_round_improvement > 0 && moves_back.load(std::memory_order_relaxed) > 0) {
            sub_round_improvement += applyMovesSortedByGainAndRevertUnbalanced(phg);
          }
          auto t4 = tbb::tick_count::now();
          if (log) LOG << "apply by prefix" << (t3-t2).seconds();
          if (log) LOG << "apply by gain and seq revert" << (t4-t3).seconds();
        }
        if (log) LOG << V(sub_round_improvement) << V(num_moves_in_sub_round);
        round_improvement += sub_round_improvement;
        num_moves += num_moves_in_sub_round;
      }
      overall_improvement += round_improvement;

      if (num_moves == 0) {
        break; // no vertices with positive gain --> stop
      }
    }

    best_metrics.km1 -= overall_improvement;
    best_metrics.imbalance = metrics::imbalance(phg, context);
    if (context.type == kahypar::ContextType::main) {
      DBG << V(best_metrics.km1) << V(best_metrics.imbalance);
    }
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
    const size_t num_moves = moves_back.load(std::memory_order_relaxed);
    tbb::parallel_sort(moves.begin(), moves.begin() + num_moves, [](const Move& m1, const Move& m2) {
      return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
    });
    
    const auto& max_part_weights = context.partition.max_part_weights;
    size_t num_overloaded_blocks = 0, num_overloaded_before_round = 0;
    vec<HypernodeWeight> part_weights = aggregatePartWeightDeltas(phg, moves, num_moves);
    for (PartitionID i = 0; i < phg.k(); ++i) {
      part_weights[i] += phg.partWeight(i);
      if (part_weights[i] > max_part_weights[i]) {
        num_overloaded_blocks++;
      }
      if (phg.partWeight(i) > max_part_weights[i]) {
        num_overloaded_before_round++;
      }
    }

    size_t num_overloaded_before_first_pass = num_overloaded_blocks;
    size_t num_reverted_moves = 0;
    size_t j = num_moves;

    auto revert_move = [&](Move& m) {
      part_weights[m.to] -= phg.nodeWeight(m.node);
      part_weights[m.from] += phg.nodeWeight(m.node);
      m.invalidate();
      num_reverted_moves++;
      if (part_weights[m.to] <= max_part_weights[m.to]) {
        num_overloaded_blocks--;
      }
    };

    while (num_overloaded_blocks > 0 && j > 0) {
      Move& m = moves[--j];
      if (part_weights[m.to] > max_part_weights[m.to]
          && part_weights[m.from] + phg.nodeWeight(m.node) <= max_part_weights[m.from]) {
        revert_move(m);
      }
    }

    if (num_overloaded_blocks > 0) {
      DBG << "still overloaded" << num_overloaded_blocks << V(num_moves) << V(num_reverted_moves)
          << V(num_overloaded_before_round) << V(num_overloaded_before_first_pass) << "trigger second run";

      size_t num_extra_rounds = 1;
      j = num_moves;
      size_t last_valid_move = 0;
      while (num_overloaded_blocks > 0) {
        if (j == 0) {
          j = last_valid_move;
          last_valid_move = 0;
          num_extra_rounds++;
        }
        Move& m = moves[j-1];
        if (m.isValid() && part_weights[m.to] > max_part_weights[m.to]) {
          if (part_weights[m.from] + phg.nodeWeight(m.node) > max_part_weights[m.from]
              && part_weights[m.from] <= max_part_weights[m.from]) {
            num_overloaded_blocks++;
          }
          revert_move(m);
        }

        if (last_valid_move == 0 && m.isValid()) {
          last_valid_move = j;
        }
        --j;
      }

      DBG << V(num_reverted_moves) << V(num_extra_rounds);
    }

    // apply all moves that were not invalidated
    Gain gain = applyMovesIf(phg, moves, num_moves, [&](size_t pos) { return moves[pos].isValid(); });

    // if that decreased solution quality, revert it all
    if (gain < 0) {
      DBG << "Kommando zurück" << V(gain) << V(num_moves) << V(num_reverted_moves);
      gain += applyMovesIf(phg, moves, num_moves, [&](size_t pos) {
        if (moves[pos].isValid()) {
          std::swap(moves[pos].from, moves[pos].to);
          return true;
        } else {
          return false;
        }
      });
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

    const size_t num_moves = moves_back.load(std::memory_order_relaxed);

    MovesWrapper moves_wrapper { moves, num_moves };

    Gain estimated_gain = 0;
    for (size_t i = 0; i < moves_wrapper.size(); ++i) {
      estimated_gain += moves[i].gain;
    }

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
        if (has_moves(p1,p2) || has_moves(p2,p1)) {
          relevant_block_pairs.emplace_back(p1, p2);
        }
        // more involvements reduce slack --> only increment involvements if vertices are moved into that block
        if (has_moves(p1,p2)) {
          involvements[p2]++;
        }
        if (has_moves(p2,p1)) {
          involvements[p1]++;
        }
      }
    }

    // swap_prefix[index(p1,p2)] stores the first position of moves to revert out of the sequence of moves from p1 to p2
    vec<size_t> swap_prefix(max_key, 0);
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
      HypernodeWeight slack_p1 = budget_p1 / std::max(1UL, involvements[p1]),
                      slack_p2 = budget_p2 / std::max(1UL, involvements[p2]);

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
      }
      if (i == i_last || sorted_moves[i].gain == 0) {
        while (j < j_last && -balance <= slack_p1 && (balance > 0 || sorted_moves[j].gain > 0)) {
          balance -= phg.nodeWeight(sorted_moves[j++].node);
          if (-balance <= slack_p1 && balance <= slack_p2) {
            best = {i,j,balance};
          }
        }
      }

      swap_prefix[index(p1,p2)] = std::get<0>(best);
      swap_prefix[index(p2,p1)] = std::get<1>(best);
      int64_t best_balance = std::get<2>(best);

      // balance < 0 --> p1 got more weight, balance > 0 --> p2 got more weight
      __atomic_fetch_add(&part_weight_deltas[p1], best_balance, __ATOMIC_RELAXED);
      __atomic_fetch_sub(&part_weight_deltas[p2], best_balance, __ATOMIC_RELAXED);
    });

    moves_back.store(0, std::memory_order_relaxed);
    Gain actual_gain = applyMovesIf(phg, sorted_moves, num_moves, [&](size_t pos) {
      if (pos < swap_prefix[index(sorted_moves[pos].from, sorted_moves[pos].to)]) {
        return true;
      } else {
        size_t second_try_pos = moves_back.fetch_add(1, std::memory_order_relaxed);
        moves[second_try_pos] = sorted_moves[pos];
        return false;
      }
    });

    // revert everything if that decreased solution quality
    if (actual_gain < 0) {
      DBG << "Kommando zurück" << V(actual_gain);
      actual_gain += applyMovesIf(phg, sorted_moves, num_moves, [&](size_t pos) {
        if (pos < swap_prefix[index(sorted_moves[pos].from, sorted_moves[pos].to)]) {
          std::swap(sorted_moves[pos].from, sorted_moves[pos].to);
          return true;
        } else {
          return false;
        }
      });

      assert(actual_gain == 0);
    }

    DBG << V(num_moves) << V(estimated_gain) << V(actual_gain) << V(metrics::imbalance(phg, context));
    return actual_gain;
  }
}