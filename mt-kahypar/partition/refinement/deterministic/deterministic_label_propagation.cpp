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
                                                        const vec<HypernodeID>&,
                                                        kahypar::Metrics& best_metrics,
                                                        const double) {
    Gain overall_improvement = 0;
    constexpr size_t num_buckets = utils::ParallelPermutation<HypernodeID>::num_buckets;
    const size_t num_sub_rounds = context.refinement.deterministic_refinement.num_sub_rounds_sync_lp;
    const size_t num_buckets_per_sub_round = parallel::chunking::idiv_ceil(num_buckets, num_sub_rounds);

    const bool log = false && context.type == kahypar::ContextType::main;

    for (size_t iter = 0; iter < context.refinement.label_propagation.maximum_iterations; ++iter) {
      if (context.refinement.deterministic_refinement.use_active_node_set && ++round == 0) {
        std::fill(last_moved_in_round.begin(), last_moved_in_round.end(), CAtomic<uint32_t>(0));
      }

      // size == 0 means no node was moved last round, but there were positive gains --> try again with different permutation
      if (!context.refinement.deterministic_refinement.use_active_node_set || iter == 0 || active_nodes.size() == 0) {
        permutation.random_grouping(phg.initialNumNodes(), context.shared_memory.static_balancing_work_packages,prng());
      } else {
        tbb::parallel_sort(active_nodes.begin(), active_nodes.end());
        permutation.sample_buckets_and_group_by(active_nodes.range(),
                                                context.shared_memory.static_balancing_work_packages, prng());
      }
      active_nodes.clear();

      size_t num_moves = 0;
      Gain round_improvement = 0;
      for (size_t sub_round = 0; sub_round < num_sub_rounds; ++sub_round) {
        auto[first_bucket, last_bucket] = parallel::chunking::bounds(sub_round, num_buckets, num_buckets_per_sub_round);
        assert(first_bucket < last_bucket && last_bucket < permutation.bucket_bounds.size());
        size_t first = permutation.bucket_bounds[first_bucket], last = permutation.bucket_bounds[last_bucket];
        moves.clear();

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
        moves.finalize();
        if (log) LOG << "calc moves time" << (tbb::tick_count::now() - t1).seconds();

        Gain sub_round_improvement = 0;
        size_t num_moves_in_sub_round = moves.size();
        DBG << V(num_moves_in_sub_round);
        if (num_moves_in_sub_round > 0) {
          auto t2 = tbb::tick_count::now();
          sub_round_improvement = applyMovesByMaximalPrefixesInBlockPairs(phg);
          auto t3 = tbb::tick_count::now();
          if (sub_round_improvement > 0 && moves.size() > 0) {
            // sub_round_improvement += applyMovesSortedByGainAndRevertUnbalanced(phg);
            sub_round_improvement += applyMovesSortedByGainWithRecalculation(phg);
          }
          auto t4 = tbb::tick_count::now();
          if (log) LOG << "apply by prefix" << (t3 - t2).seconds();
          if (log) LOG << "apply by gain and seq revert" << (t4 - t3).seconds();
        }
        if (log) LOG << V(sub_round_improvement) << V(num_moves_in_sub_round);
        round_improvement += sub_round_improvement;
        num_moves += num_moves_in_sub_round;
      }
      overall_improvement += round_improvement;
      active_nodes.finalize();

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
 */
  Gain DeterministicLabelPropagationRefiner::performMoveWithAttributedGain(
          PartitionedHypergraph& phg, const Move& m, bool activate_neighbors) {
    Gain attributed_gain = 0;
    auto objective_delta = [&](HyperedgeID he, HyperedgeWeight edge_weight, HypernodeID edge_size,
                               HypernodeID pin_count_in_from_part_after, HypernodeID pin_count_in_to_part_after) {
      attributed_gain -= km1Delta(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
    };
    const bool was_moved = phg.changeNodePart(m.node, m.from, m.to, objective_delta);
    if (context.refinement.deterministic_refinement.use_active_node_set && activate_neighbors && was_moved) {
      // activate neighbors for next round
      const HypernodeID n = phg.initialNumNodes();
      for (HyperedgeID he : phg.incidentEdges(m.node)) {
        if (phg.edgeSize(he) <= context.refinement.label_propagation.hyperedge_size_activation_threshold) {
          if (last_moved_in_round[he + n].load(std::memory_order_relaxed) != round) {
            last_moved_in_round[he + n].store(round, std::memory_order_relaxed);   // no need for atomic semantics
            for (HypernodeID v : phg.pins(he)) {
              uint32_t lrv = last_moved_in_round[v].load(std::memory_order_relaxed);
              if (lrv != round &&
                  last_moved_in_round[v].compare_exchange_strong(lrv, round, std::memory_order_acq_rel)) {
                active_nodes.push_back_buffered(v);
              }
            }
          }
        }
      }
    }
    return attributed_gain;
  }

  template<typename Predicate>
  Gain DeterministicLabelPropagationRefiner::applyMovesIf(
          PartitionedHypergraph& phg, const vec<Move>& my_moves, size_t end, Predicate&& predicate) {
    auto range = tbb::blocked_range<size_t>(0UL, end);
    auto accum = [&](const tbb::blocked_range<size_t>& r, const Gain& init) -> Gain {
      Gain my_gain = init;
      for (size_t i = r.begin(); i < r.end(); ++i) {
        if (predicate(i)) {
          my_gain += performMoveWithAttributedGain(phg, my_moves[i], true);
        }
      }
      return my_gain;
    };
    return tbb::parallel_reduce(range, 0, accum, std::plus<>());
  }

  vec<HypernodeWeight> aggregatePartWeightDeltas(PartitionedHypergraph& phg, const vec<Move>& moves, size_t end) {
    // parallel reduce makes way too many vector copies
    tbb::enumerable_thread_specific<vec< HypernodeWeight>>
    ets_part_weight_diffs(phg.k(), 0);
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
    const size_t num_moves = moves.size();
    tbb::parallel_sort(moves.begin(), moves.begin() + num_moves, [](const Move& m1, const Move& m2) {
      return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
    });

    const auto& max_part_weights = context.partition.max_part_weights;
    size_t num_overloaded_blocks = 0, num_overloaded_before_round = 0;
    vec<HypernodeWeight> part_weights = aggregatePartWeightDeltas(phg, moves.getData(), num_moves);
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
        Move& m = moves[j - 1];
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
    Gain gain = applyMovesIf(phg, moves.getData(), num_moves, [&](size_t pos) { return moves[pos].isValid(); });

    // if that decreased solution quality, revert it all
    if (gain < 0) {
      DBG << "Kommando zurück" << V(gain) << V(num_moves) << V(num_reverted_moves);
      gain += applyMovesIf(phg, moves.getData(), num_moves, [&](size_t pos) {
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
    PartitionID max_key = k * k;
    auto index = [&](PartitionID b1, PartitionID b2) { return b1 * k + b2; };
    auto get_key = [&](const Move& m) { return index(m.from, m.to); };

    const size_t num_moves = moves.size();

    // aggregate moves by direction. not in-place because of counting sort.
    // but it gives us the positions of the buckets right away
    auto positions = parallel::counting_sort(moves, sorted_moves, max_key, get_key,
                                             context.shared_memory.num_threads);

    auto has_moves = [&](PartitionID p1, PartitionID p2) {
      size_t direction = index(p1, p2);
      return positions[direction + 1] != positions[direction];
    };

    vec<std::pair<PartitionID, PartitionID>> relevant_block_pairs;
    vec<size_t> involvements(k, 0);
    for (PartitionID p1 = 0; p1 < k; ++p1) {
      for (PartitionID p2 = p1 + 1; p2 < k; ++p2) {
        if (has_moves(p1, p2) || has_moves(p2, p1)) {
          relevant_block_pairs.emplace_back(p1, p2);
        }
        // more involvements reduce slack --> only increment involvements if vertices are moved into that block
        if (has_moves(p1, p2)) {
          involvements[p2]++;
        }
        if (has_moves(p2, p1)) {
          involvements[p1]++;
        }
      }
    }

    // swap_prefix[index(p1,p2)] stores the first position of moves to revert out of the sequence of moves from p1 to p2
    vec<size_t> swap_prefix(max_key, 0);
    vec<HypernodeWeight> part_weights(k, 0);
    for (PartitionID i = 0; i < phg.k(); ++i) { part_weights[i] = phg.partWeight(i); }

    tbb::parallel_for(0UL, relevant_block_pairs.size(), [&](size_t bp_index) {
      // sort both directions by gain (alternative: gain / weight?)
      auto sort_by_gain_and_prefix_sum_node_weights = [&](PartitionID p1, PartitionID p2) {
        size_t begin = positions[index(p1, p2)], end = positions[index(p1, p2) + 1];
        auto comp = [&](const Move& m1, const Move& m2) {
          return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
        };
        tbb::parallel_sort(sorted_moves.begin() + begin, sorted_moves.begin() + end, comp);
        tbb::parallel_for(begin, end, [&](size_t pos) {
          cumulative_node_weights[pos] = phg.nodeWeight(sorted_moves[pos].node);
        });
        parallel_prefix_sum(cumulative_node_weights.begin() + begin, cumulative_node_weights.begin() + end,
                            cumulative_node_weights.begin() + begin, std::plus<>(), 0);
      };

      const auto& bp = relevant_block_pairs[bp_index];
      tbb::parallel_invoke([&] {
        sort_by_gain_and_prefix_sum_node_weights(bp.first, bp.second);
      }, [&] {
        sort_by_gain_and_prefix_sum_node_weights(bp.second, bp.first);
      });
    });

    for (const auto& [p1,p2] : relevant_block_pairs) {
      HypernodeWeight lb_p1 = (part_weights[p1]) - context.partition.max_part_weights[p1],
                      ub_p2 = context.partition.max_part_weights[p2] - (part_weights[p2]);

      size_t p1_begin = positions[index(p1, p2)], p1_end = positions[index(p1, p2) + 1],
             p2_begin = positions[index(p2, p1)], p2_end = positions[index(p2, p1) + 1];

      auto balance = [&](size_t p1_ind, size_t p2_ind) {
        const auto a = (p1_ind == p1_begin - 1) ? 0 : cumulative_node_weights[p1_ind];
        const auto b = (p2_ind == p2_begin - 1) ? 0 : cumulative_node_weights[p2_ind];
        return a - b;
      };

      LOG << V(p1) << V(p2) << V(p1_begin) << V(p1_end) << V(p2_begin) << V(p2_begin) << (lb_p1) << V(ub_p2);

      auto best_prefix = findBestPrefixesRecursive(p1_begin, p1_end, p2_begin, p2_end,
                                                   p1_begin - 1, p2_begin - 1, lb_p1, ub_p2);


      swap_prefix[index(p1, p2)] = best_prefix.first;
      swap_prefix[index(p2, p1)] = best_prefix.second;
      HypernodeWeight best_balance = balance(best_prefix.first - 1, best_prefix.second - 1);
      LOG << V(best_prefix.first) << V(best_prefix.second) << V(best_balance);
      __atomic_fetch_sub(&part_weights[p1], best_balance, __ATOMIC_RELAXED);
      __atomic_fetch_add(&part_weights[p2], best_balance, __ATOMIC_RELAXED);
    }
    LOG << "block-pair handling done";

    moves.clear();
    Gain actual_gain = applyMovesIf(phg, sorted_moves, num_moves, [&](size_t pos) {
      if (pos < swap_prefix[index(sorted_moves[pos].from, sorted_moves[pos].to)]) {
        return true;
      } else {
        moves.push_back_buffered(sorted_moves[pos]);
        return false;
      }
    });
    moves.finalize();

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

    DBG << V(num_moves) << V(actual_gain) << V(metrics::imbalance(phg, context));
    return actual_gain;
  }

  std::pair<size_t, size_t> DeterministicLabelPropagationRefiner::findBestPrefixesRecursive(
          size_t p1_begin, size_t p1_end, size_t p2_begin, size_t p2_end,
          size_t p1_invalid, size_t p2_invalid,
          HypernodeWeight lb_p1, HypernodeWeight ub_p2)
  {
    auto balance = [&](size_t p1_ind, size_t p2_ind) {
      const auto a = (p1_ind == p1_invalid) ? 0 : cumulative_node_weights[p1_ind];
      const auto b = (p2_ind == p2_invalid) ? 0 : cumulative_node_weights[p2_ind];
      return a - b;
    };

    auto is_feasible = [&](size_t p1_ind, size_t p2_ind) {
      const HypernodeWeight bal = balance(p1_ind, p2_ind);
      return lb_p1 <= bal && bal <= ub_p2;
    };

    const size_t n_p1 = p1_end - p1_begin, n_p2 = p2_end - p2_begin;

    static constexpr size_t sequential_cutoff = 5000;
    if (n_p1 < sequential_cutoff && n_p2 < sequential_cutoff) {
      while (true) {
        if (is_feasible(p1_end - 1, p2_end - 1)) { return std::make_pair(p1_end, p2_end); }
        if (balance(p1_end - 1, p2_end - 1) < 0) {
          if (p2_end == p2_begin) { break; }
          p2_end--;
        } else {
          if (p1_end == p1_begin) { break; }
          p1_end--;
        }
      }
      return std::make_pair(invalid_pos, invalid_pos);
    }

    const auto c = cumulative_node_weights.begin();
    if (n_p1 > n_p2) {
      size_t p1_mid = p1_begin + n_p1 / 2;
      auto p2_match_it = std::lower_bound(c + p2_begin, c + p2_end, cumulative_node_weights[p1_mid]);
      size_t p2_match = std::distance(cumulative_node_weights.begin(), p2_match_it);

      if (p2_match != p2_end && p1_mid != p1_end && is_feasible(p1_mid, p2_match)) {
        // no need to search left range
        return findBestPrefixesRecursive(p1_mid + 1, p1_end, p2_match + 1, p2_end, p1_invalid, p2_invalid, lb_p1, ub_p2);
      }

      // TODO apply pruning
      // i.e. if p1_mid cannot be compensated, don't recurse on the tail range
      // could also binary search for the end of a range + slack and discard everything behind that

      std::pair<size_t, size_t> left, right;
      tbb::parallel_invoke([&] {
        left = findBestPrefixesRecursive(p1_begin, p1_mid, p2_begin, p2_match, p1_invalid, p2_invalid, lb_p1, ub_p2);
      }, [&] {
        right = findBestPrefixesRecursive(p1_mid, p1_end, p2_match, p2_end, p1_invalid, p2_invalid, lb_p1, ub_p2);
      });

      return right.first != invalid_pos ? right : left;
    } else {

    }

  }

  Gain DeterministicLabelPropagationRefiner::applyMovesSortedByGainWithRecalculation(PartitionedHypergraph& phg) {
    if (last_recalc_round.empty() || ++recalc_round == std::numeric_limits<uint32_t>::max()) {
      last_recalc_round.assign(max_num_edges, CAtomic<uint32_t>(0));
    }
    constexpr MoveID invalid_move_id = std::numeric_limits<MoveID>::max();
    if (move_pos_of_node.empty()) {
      move_pos_of_node.resize(max_num_nodes, invalid_move_id);
    }

    const size_t num_moves = moves.size();
    tbb::parallel_sort(moves.begin(), moves.begin() + num_moves, [](const Move& m1, const Move& m2) {
      return m1.gain > m2.gain || (m1.gain == m2.gain && m1.node < m2.node);
    });

    tbb::parallel_for(0UL, num_moves, [&](size_t pos) {
      move_pos_of_node[moves[pos].node] = pos + 1;    // pos + 1 to handle zero init of last_out
      moves[pos].gain = 0;
    });

    auto was_node_moved_in_this_round = [&](HypernodeID u) { return move_pos_of_node[u] != invalid_move_id; };

    // recalculate gains
    tbb::parallel_for(0UL, num_moves, [&](size_t pos) {
      auto& r = ets_recalc_data.local();

      HypernodeID u = moves[pos].node;
      for (HyperedgeID e : phg.incidentEdges(u)) {
        uint32_t expected = last_recalc_round[e].load(std::memory_order_relaxed);
        if (expected < recalc_round && last_recalc_round[e].exchange(recalc_round, std::memory_order_acq_rel) == expected) {
          for (HypernodeID v : phg.pins(e)) {
            if (was_node_moved_in_this_round(v)) {
              const MoveID m_id = move_pos_of_node[v];
              const Move& m = moves[m_id - 1];
              r[m.to].first_in = std::min(r[m.to].first_in, m_id);
              r[m.from].last_out = std::max(r[m.from].last_out, m_id);
            } else {
              r[phg.partID(v)].remaining_pins++;
            }
          }

          const HyperedgeWeight we = phg.edgeWeight(e);
          for (HypernodeID v : phg.pins(e)) {
            if (was_node_moved_in_this_round(v)) {
              const MoveID m_id = move_pos_of_node[v];
              Move& m = moves[m_id - 1];
              const bool benefit = r[m.from].last_out == m_id && r[m.from].first_in > m_id && r[m.from].remaining_pins == 0;
              const bool penalty = r[m.to].first_in == m_id && r[m.to].last_out < m_id && r[m.to].remaining_pins == 0;
              if (benefit && !penalty) {
                __atomic_fetch_add(&m.gain, we, __ATOMIC_RELAXED);
              }
              if (!benefit && penalty) {
                __atomic_fetch_sub(&m.gain, we, __ATOMIC_RELAXED);
              }
            }
          }

          if (phg.k() <= static_cast<int>(2 * phg.edgeSize(e))) {
            for (PartitionID i = 0; i < phg.k(); ++i) {
              r[i] = RecalculationData();
            }
          } else {
            for (HypernodeID v : phg.pins(e)) {
              if (was_node_moved_in_this_round(v)) {
                const Move& m = moves[move_pos_of_node[v] - 1];
                r[m.from] = RecalculationData();
                r[m.to] = RecalculationData();
              } else {
                r[phg.partID(v)] = RecalculationData();
              }
            }
          }
        }
      }
    });

#ifndef NDEBUG
    for (size_t pos = 0; pos < num_moves; ++pos) {
      const Move& m = moves[pos];
      Gain move_gain = performMoveWithAttributedGain(phg, m, false);
      ASSERT(move_gain == m.gain);
    }

    for (int64_t pos = num_moves - 1; pos >= 0; --pos) {
      Move reverse_move = moves[pos];
      std::swap(reverse_move.from, reverse_move.to);
      Gain move_gain = performMoveWithAttributedGain(phg, reverse_move, false);
      ASSERT(move_gain == -moves[pos].gain);
    }
#endif

    // remove markers again
    tbb::parallel_for(0UL, num_moves, [&](size_t pos) { move_pos_of_node[moves[pos].node] = invalid_move_id; });

    // calculate number of overloaded blocks
    size_t num_overloaded_blocks_before_pass = 0, num_overloaded_blocks = 0;
    const auto& max_part_weights = context.partition.max_part_weights;
    vec<HypernodeWeight> part_weights(phg.k());
    for (PartitionID i = 0; i < phg.k(); ++i) {
      part_weights[i] = phg.partWeight(i);
      if (part_weights[i] > max_part_weights[i]) {
        num_overloaded_blocks_before_pass++;
      }
    }
    num_overloaded_blocks = num_overloaded_blocks_before_pass;

    // prefix sum part weights and gains. (might incorporate parallel version if this takes too long)
    Gain best_gain = 0, gain_sum = 0;
    size_t best_index = 0;
    for (size_t pos = 0; pos < num_moves; ++pos) {
      const Move& m = moves[pos];
      num_overloaded_blocks -= (part_weights[m.from] > max_part_weights[m.from] &&
                                part_weights[m.from] - phg.nodeWeight(m.node) <= max_part_weights[m.from]);
      num_overloaded_blocks += (part_weights[m.to] <= max_part_weights[m.to] &&
                                part_weights[m.to] + phg.nodeWeight(m.node) > max_part_weights[m.to]);

      part_weights[m.from] -= phg.nodeWeight(m.node);
      part_weights[m.to] += phg.nodeWeight(m.node);
      gain_sum += m.gain;
      if (num_overloaded_blocks <= num_overloaded_blocks_before_pass && gain_sum >= best_gain) {
        best_index = pos + 1;
        best_gain = gain_sum;
      }
    }

    Gain attributed_gain = applyMovesIf(phg, moves.getData(), best_index, [&](size_t) { return true; });
    ASSERT(attributed_gain == best_gain); unused(attributed_gain);

    return best_gain;
  }

} // namespace mt_kahypar