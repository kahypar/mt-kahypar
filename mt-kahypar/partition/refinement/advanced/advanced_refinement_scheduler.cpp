/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/advanced/advanced_refinement_scheduler.h"

#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

void AdvancedRefinementScheduler::RefinementStats::update_global_stats() {
  utils::Stats& global_stats = utils::Stats::instance();
  global_stats.update_stat("num_advanced_refinements",
    num_refinements.load(std::memory_order_relaxed));
  global_stats.update_stat("num_advanced_improvement",
    num_improvements.load(std::memory_order_relaxed));
  global_stats.update_stat("correct_expected_improvement",
    correct_expected_improvement.load(std::memory_order_relaxed));
  global_stats.update_stat("failed_updates_due_to_conflicting_moves",
    failed_updates_due_to_conflicting_moves.load(std::memory_order_relaxed));
  global_stats.update_stat("failed_updates_due_to_conflicting_moves_without_rollback",
    failed_updates_due_to_conflicting_moves_without_rollback.load(std::memory_order_relaxed));
  global_stats.update_stat("failed_updates_due_to_balance_constraint",
    failed_updates_due_to_balance_constraint.load(std::memory_order_relaxed));
  global_stats.update_stat("total_advanced_refinement_improvement",
    total_improvement.load(std::memory_order_relaxed));
}

bool AdvancedRefinementScheduler::refineImpl(
                PartitionedHypergraph& phg,
                const parallel::scalable_vector<HypernodeID>&,
                kahypar::Metrics&,
                const double)  {
  unused(phg);
  ASSERT(_phg == &phg);

  _stats.update_global_stats();
  if ( _context.partition.paradigm == Paradigm::nlevel && phg.isGainCacheInitialized() ) {
    // Update Gain Cache
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( _was_moved[hn] ) {
        phg.recomputeMoveFromBenefit(hn);
        _was_moved[hn] = uint8_t(false);
      }
    });
  }
  HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
  _phg = nullptr;
  return false;
}

void AdvancedRefinementScheduler::initializeImpl(PartitionedHypergraph& phg)  {
  _phg = &phg;

  // Initialize Part Weights
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    _part_weights[i] = phg.partWeight(i);
  }

  utils::Timer::instance().start_timer("initialize_quotient_graph", "Initialize Quotient Graph");
  _quotient_graph.initialize(phg);
  utils::Timer::instance().stop_timer("initialize_quotient_graph");
}

namespace {

template<typename F>
bool changeNodePart(PartitionedHypergraph& phg,
                    const HypernodeID hn,
                    const PartitionID from,
                    const PartitionID to,
                    const F& objective_delta,
                    const bool is_nlevel) {
  bool success = false;
  if ( is_nlevel && phg.isGainCacheInitialized()) {
    success = phg.changeNodePartWithGainCacheUpdate(hn, from, to,
      std::numeric_limits<HypernodeWeight>::max(), [] { }, objective_delta);
  } else {
    success = phg.changeNodePart(hn, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, objective_delta);
  }
  ASSERT(success);
  return success;
}

template<typename F>
void applyMoveSequence(PartitionedHypergraph& phg,
                       const MoveSequence& sequence,
                       const F& objective_delta,
                       const bool is_nlevel,
                       vec<uint8_t>& was_moved) {
  for ( const Move& move : sequence.moves ) {
    changeNodePart(phg, move.node, move.from, move.to, objective_delta, is_nlevel);
    was_moved[move.node] = uint8_t(true);
  }
}

template<typename F>
void revertMoveSequence(PartitionedHypergraph& phg,
                        const MoveSequence& sequence,
                        const F& objective_delta,
                        const bool is_nlevel) {
  for ( const Move& move : sequence.moves ) {
    ASSERT(phg.partID(move.node) == move.to);
    changeNodePart(phg, move.node, move.to, move.from, objective_delta, is_nlevel);
  }
}

} // namespace

HyperedgeWeight AdvancedRefinementScheduler::applyMoves(MoveSequence& sequence) {
  ASSERT(_phg);
  ++_stats.num_refinements;

  // Compute Part Weight Deltas
  vec<HypernodeWeight> part_weight_deltas(_context.partition.k, 0);
  for ( const Move& move : sequence.moves ) {
    ASSERT(move.from != kInvalidPartition && move.from < _context.partition.k);
    ASSERT(move.to != kInvalidPartition && move.to < _context.partition.k);
    ASSERT(move.from != move.to);
    ASSERT(_phg->partID(move.node) == move.from);
    const HypernodeWeight node_weight = _phg->nodeWeight(move.node);
    part_weight_deltas[move.from] -= node_weight;
    part_weight_deltas[move.to] += node_weight;
  }

  HyperedgeWeight improvement = 0;
  auto delta_func = [&](const HyperedgeID he,
                        const HyperedgeWeight edge_weight,
                        const HypernodeID edge_size,
                        const HypernodeID pin_count_in_from_part_after,
                        const HypernodeID pin_count_in_to_part_after) {
    if ( _context.partition.objective == kahypar::Objective::km1 ) {
      improvement -= km1Delta(he, edge_weight, edge_size,
        pin_count_in_from_part_after, pin_count_in_to_part_after);
    } else if ( _context.partition.objective == kahypar::Objective::cut ) {
      improvement -= cutDelta(he, edge_weight, edge_size,
        pin_count_in_from_part_after, pin_count_in_to_part_after);
    }
  };


  // Update part weights atomically
  bool is_balanced = partWeightUpdate(part_weight_deltas, false);
  if ( is_balanced ) {
    // Apply move sequence to partition
    const bool is_nlevel = _context.partition.paradigm == Paradigm::nlevel;
    applyMoveSequence(*_phg, sequence, delta_func, is_nlevel, _was_moved);

    if ( improvement < 0 ) {
      is_balanced = partWeightUpdate(part_weight_deltas, true);
      if ( is_balanced ) {
        // Move sequence worsen solution quality => Rollback
        revertMoveSequence(*_phg, sequence, delta_func, is_nlevel);
        ++_stats.failed_updates_due_to_conflicting_moves;
        sequence.state = MoveSequenceState::WORSEN_SOLUTION_QUALITY;
      } else {
        // Rollback would violate balance constraint => Worst Case
        ++_stats.failed_updates_due_to_conflicting_moves_without_rollback;
        sequence.state = MoveSequenceState::WORSEN_SOLUTION_QUALITY_WITHOUT_ROLLBACK;
      }
    } else {
      ++_stats.num_improvements;
      _stats.correct_expected_improvement += (improvement == sequence.expected_improvement);
      sequence.state = MoveSequenceState::SUCCESS;
    }
  } else {
    ++_stats.failed_updates_due_to_balance_constraint;
    sequence.state = MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT;
  }

  _stats.total_improvement += improvement;
  return improvement;
}

bool AdvancedRefinementScheduler::partWeightUpdate(const vec<HypernodeWeight>& part_weight_deltas,
                                                   const bool rollback) {
  const HypernodeWeight multiplier = rollback ? -1 : 1;
  bool is_balanced = true;
  _part_weights_lock.lock();
  PartitionID i = 0;
  for ( ; i < _context.partition.k; ++i ) {
    if ( _part_weights[i] + multiplier * part_weight_deltas[i] > _context.partition.max_part_weights[i] ) {
      // Move Sequence Violates Balance Constraint => Rollback
      --i;
      for ( ; i >= 0; --i ) {
        _part_weights[i] -= multiplier * part_weight_deltas[i];
      }
      is_balanced = false;
      break;
    }
    _part_weights[i] += multiplier * part_weight_deltas[i];
  }
  _part_weights_lock.unlock();
  return is_balanced;
}

}