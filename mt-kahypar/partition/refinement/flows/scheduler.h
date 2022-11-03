/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/quotient_graph.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {

class FlowRefinementScheduler final : public IRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  struct RefinementStats {
    RefinementStats(utils::Stats& stats) :
      _stats(stats),
      num_refinements(0),
      num_improvements(0),
      num_time_limits(0),
      correct_expected_improvement(0),
      zero_gain_improvement(0),
      failed_updates_due_to_conflicting_moves(0),
      failed_updates_due_to_conflicting_moves_without_rollback(0),
      failed_updates_due_to_balance_constraint(0),
      total_improvement(0) { }

    void reset() {
      num_refinements.store(0);
      num_improvements.store(0);
      num_time_limits.store(0);
      correct_expected_improvement.store(0);
      zero_gain_improvement.store(0);
      failed_updates_due_to_conflicting_moves.store(0);
      failed_updates_due_to_conflicting_moves_without_rollback.store(0);
      failed_updates_due_to_balance_constraint.store(0);
      total_improvement.store(0);
    }

    void update_global_stats();

    utils::Stats& _stats;
    CAtomic<int64_t> num_refinements;
    CAtomic<int64_t> num_improvements;
    CAtomic<int64_t> num_time_limits;
    CAtomic<int64_t> correct_expected_improvement;
    CAtomic<int64_t> zero_gain_improvement;
    CAtomic<int64_t> failed_updates_due_to_conflicting_moves;
    CAtomic<int64_t> failed_updates_due_to_conflicting_moves_without_rollback;
    CAtomic<int64_t> failed_updates_due_to_balance_constraint;
    CAtomic<HyperedgeWeight> total_improvement;
  };

  struct PartWeightUpdateResult {
    bool is_balanced = true;
    PartitionID overloaded_block = kInvalidPartition;
    HypernodeWeight overload_weight = 0;
  };

  friend std::ostream & operator<< (std::ostream& str, const RefinementStats& stats);

public:
  explicit FlowRefinementScheduler(const Hypergraph& hg,
                                   const Context& context) :
    _phg(nullptr),
    _context(context),
    _quotient_graph(hg, context),
    _refiner(hg, context),
    _constructor(hg, context),
    _was_moved(hg.initialNumNodes(), uint8_t(false)),
    _part_weights_lock(),
    _part_weights(context.partition.k, 0),
    _max_part_weights(context.partition.k, 0),
    _stats(utils::Utilities::instance().getStats(context.utility_id)),
    _apply_moves_lock() { }

  FlowRefinementScheduler(const FlowRefinementScheduler&) = delete;
  FlowRefinementScheduler(FlowRefinementScheduler&&) = delete;

  FlowRefinementScheduler & operator= (const FlowRefinementScheduler &) = delete;
  FlowRefinementScheduler & operator= (FlowRefinementScheduler &&) = delete;

  /**
   * Applies the sequence of vertex moves to the partitioned hypergraph.
   * The method ensures that the move sequence does not violate
   * the balance constaint and not worsen solution quality.
   * Returns, improvement in solution quality.
   */
  HyperedgeWeight applyMoves(const SearchID search_id,
                             MoveSequence& sequence);

  /**
   * Returns the current weight of each block.
   * Note, we do not want that the underlying refiner (ILP and Flows)
   * see partially updated part weight information. Thus, we perform
   * part weight updates for a move sequence as a transaction, which
   * we protect with a spin lock.
   */
  vec<HypernodeWeight> partWeights() {
    _part_weights_lock.lock();
    vec<HypernodeWeight> _copy_part_weights(_part_weights);
    _part_weights_lock.unlock();
    return _copy_part_weights;
  }

private:
  bool refineImpl(PartitionedHypergraph& phg,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& metrics,
                  double time_limit) final;

  void initializeImpl(PartitionedHypergraph& phg) final;

  PartWeightUpdateResult partWeightUpdate(const vec<HypernodeWeight>& part_weight_deltas,
                                          const bool rollback);

  std::string blocksOfSearch(const SearchID search_id) {
    const BlockPair blocks = _quotient_graph.getBlockPair(search_id);
    return "(" + std::to_string(blocks.i) + "," + std::to_string(blocks.j) + ")";
  }

  PartitionedHypergraph* _phg;
  const Context& _context;

  // ! Contains information of all cut hyperedges between the
  // ! blocks of the partition
  QuotientGraph _quotient_graph;

  // ! Maintains the flow refiner instances
  FlowRefinerAdapter _refiner;

  // ! Responsible for construction of an flow problems
  ProblemConstruction _constructor;

  // ! For each vertex it store wheather the corresponding vertex
  // ! was moved or not
  vec<uint8_t> _was_moved;

  // ! Maintains the part weights of each block
  SpinLock _part_weights_lock;
  vec<HypernodeWeight> _part_weights;
  vec<HypernodeWeight> _max_part_weights;

  // ! Contains refinement statistics
  RefinementStats _stats;

  SpinLock _apply_moves_lock;
};

}  // namespace kahypar
