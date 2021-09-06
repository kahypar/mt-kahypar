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

#include "mt-kahypar/partition/refinement/advanced/scheduler.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/partitioning_output.h"
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
                kahypar::Metrics& best_metrics,
                const double)  {
  unused(phg);
  ASSERT(_phg == &phg);

  utils::Timer::instance().start_timer("advanced_refinement_scheduling", "Advanced Refinement Scheduling");
  std::atomic<HyperedgeWeight> overall_delta(0);
  tbb::parallel_for(0UL, _refiner.numAvailableRefiner(), [&](const size_t i) {
    while ( !_quotient_graph.terminate() && i < _quotient_graph.maximumRequiredRefiners() ) {
      SearchID search_id = _quotient_graph.requestNewSearch(_refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        DBG << "Start search" << search_id << "( Blocks =" << blocksOfSearch(search_id) << ", Thread =" << sched_getcpu() << ")";
        utils::Timer::instance().start_timer("construct_problem", "Construct Problem", true);
        const vec<HypernodeID> refinement_nodes =
          _constructor.construct(search_id, _quotient_graph, _refiner, phg);
        _quotient_graph.finalizeConstruction(search_id);
        utils::Timer::instance().stop_timer("construct_problem");

        bool improved_solution = false;
        if ( refinement_nodes.size() > 0 ) {
          utils::Timer::instance().start_timer("refine_problem", "Refine Problem", true);
          ++_stats.num_refinements;
          MoveSequence sequence = _refiner.refine(search_id, phg, refinement_nodes);
          utils::Timer::instance().stop_timer("refine_problem");

          if ( !sequence.moves.empty() ) {
            utils::Timer::instance().start_timer("apply_moves", "Apply Moves", true);
            HyperedgeWeight delta = applyMoves(search_id, sequence);
            overall_delta -= delta;
            improved_solution = sequence.state == MoveSequenceState::SUCCESS && delta > 0;
            utils::Timer::instance().stop_timer("apply_moves");
          }
        }

        _constructor.releaseNodes(search_id, refinement_nodes);
        _quotient_graph.finalizeSearch(search_id, improved_solution);
        _refiner.finalizeSearch(search_id);
      }
    }
    DBG << RED << "Refiner" << i << "terminates!" << END;
  });
  utils::Timer::instance().stop_timer("advanced_refinement_scheduling");

  DBG << "Num Refinements =" << _stats.num_refinements;
  DBG << "Num Improvements =" << _stats.num_improvements;
  DBG << "Failed due to Balance Constraint =" << _stats.failed_updates_due_to_balance_constraint;
  DBG << "Failed due to Conflicting Moves =" << _stats.failed_updates_due_to_conflicting_moves;
  DBG << "Total Improvement =" << _stats.total_improvement;
  DBG << "---------------------------------------------------------------";

  ASSERT([&]() {
    for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
      if ( _part_weights[i] != phg.partWeight(i) ) {
        LOG << V(_part_weights[i]) << V(phg.partWeight(i));
        return false;
      }
    }
    return true;
  }(), "Concurrent part weight updates failed!");

  // Update metrics statistics
  HyperedgeWeight current_metric = best_metrics.getMetric(
    kahypar::Mode::direct_kway, _context.partition.objective);
  HEAVY_REFINEMENT_ASSERT(current_metric + overall_delta ==
                          metrics::objective(phg, _context.partition.objective),
                          V(current_metric) << V(overall_delta) <<
                          V(metrics::objective(phg, _context.partition.objective)));
  best_metrics.updateMetric(current_metric + overall_delta,
    kahypar::Mode::direct_kway, _context.partition.objective);
  best_metrics.imbalance = metrics::imbalance(phg, _context);
  _stats.update_global_stats();

  // Update Gain Cache
  if ( _context.partition.paradigm == Paradigm::nlevel && phg.isGainCacheInitialized() ) {
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( _was_moved[hn] ) {
        phg.recomputeMoveFromBenefit(hn);
        _was_moved[hn] = uint8_t(false);
      }
    });
  }
  HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
  _phg = nullptr;
  return overall_delta.load(std::memory_order_relaxed) < 0;
}

void AdvancedRefinementScheduler::initializeImpl(PartitionedHypergraph& phg)  {
  _phg = &phg;

  // Initialize Part Weights
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    _part_weights[i] = phg.partWeight(i);
    _max_part_weights[i] = std::max(
      phg.partWeight(i), _context.partition.max_part_weights[i]);
  }

  _stats.reset();
  _refiner.reset();
  utils::Timer::instance().start_timer("initialize_quotient_graph", "Initialize Quotient Graph");
  _quotient_graph.initialize(phg);
  utils::Timer::instance().stop_timer("initialize_quotient_graph");
}

namespace {

struct NewCutHyperedge {
  HyperedgeID he;
  PartitionID block;
};

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
                       vec<uint8_t>& was_moved,
                       vec<NewCutHyperedge>& new_cut_hes) {
  for ( const Move& move : sequence.moves ) {
    ASSERT(move.from == phg.partID(move.node));
    if ( move.from != move.to ) {
      changeNodePart(phg, move.node, move.from, move.to, objective_delta, is_nlevel);
      was_moved[move.node] = uint8_t(true);
      // If move increases the pin count of some hyperedges in block 'move.to' to one 1
      // we set the corresponding block here.
      int i = new_cut_hes.size() - 1;
      while ( i >= 0 && new_cut_hes[i].block == kInvalidPartition ) {
        new_cut_hes[i].block = move.to;
        --i;
      }
    }
  }
}

template<typename F>
void revertMoveSequence(PartitionedHypergraph& phg,
                        const MoveSequence& sequence,
                        const F& objective_delta,
                        const bool is_nlevel) {
  for ( const Move& move : sequence.moves ) {
    if ( move.from != move.to ) {
      ASSERT(phg.partID(move.node) == move.to);
      changeNodePart(phg, move.node, move.to, move.from, objective_delta, is_nlevel);
    }
  }
}

void addCutHyperedgesToQuotientGraph(QuotientGraph& quotient_graph,
                                     const vec<NewCutHyperedge>& new_cut_hes) {
  for ( const NewCutHyperedge& new_cut_he : new_cut_hes ) {
    ASSERT(new_cut_he.block != kInvalidPartition);
    quotient_graph.addNewCutHyperedge(new_cut_he.he, new_cut_he.block);
  }
}

} // namespace

HyperedgeWeight AdvancedRefinementScheduler::applyMoves(const SearchID search_id,
                                                        MoveSequence& sequence,
                                                        const size_t num_retries) {
  unused(search_id);
  ASSERT(_phg);

  // TODO: currently we lock the applyMoves method, if overlapping search is enabled.
  // => find something smarter here
  if ( _context.refinement.advanced.use_overlapping_searches ) {
    _apply_moves_lock.lock();
  }

  // Compute Part Weight Deltas
  vec<HypernodeWeight> part_weight_deltas(_context.partition.k, 0);
  for ( Move& move : sequence.moves ) {
    move.from = _phg->partID(move.node);
    if ( move.from != move.to ) {
      const HypernodeWeight node_weight = _phg->nodeWeight(move.node);
      part_weight_deltas[move.from] -= node_weight;
      part_weight_deltas[move.to] += node_weight;
    }
  }

  HyperedgeWeight improvement = 0;
  vec<NewCutHyperedge> new_cut_hes;
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

    // Collect hyperedges with new blocks in its connectivity set
    if ( pin_count_in_to_part_after == 1 ) {
      // the corresponding block will be set in applyMoveSequence(...) function
      new_cut_hes.emplace_back(NewCutHyperedge { he, kInvalidPartition });
    }
  };

  // Update part weights atomically
  PartWeightUpdateResult update_res = partWeightUpdate(part_weight_deltas, false);
  if ( update_res.is_balanced ) {
    // Apply move sequence to partition
    const bool is_nlevel = _context.partition.paradigm == Paradigm::nlevel;
    applyMoveSequence(*_phg, sequence, delta_func, is_nlevel, _was_moved, new_cut_hes);

    if ( improvement < 0 ) {
      update_res = partWeightUpdate(part_weight_deltas, true);
      if ( update_res.is_balanced ) {
        // Move sequence worsen solution quality => Rollback
        DBG << RED << "Move sequence worsen solution quality ("
            << "Expected Improvement =" << sequence.expected_improvement
            << ", Real Improvement =" << improvement
            << ", Search ID =" << search_id << ")" << END;
        revertMoveSequence(*_phg, sequence, delta_func, is_nlevel);
        ++_stats.failed_updates_due_to_conflicting_moves;
        sequence.state = MoveSequenceState::WORSEN_SOLUTION_QUALITY;
      } else {
        // Rollback would violate balance constraint => Worst Case
        ++_stats.failed_updates_due_to_conflicting_moves_without_rollback;
        sequence.state = MoveSequenceState::WORSEN_SOLUTION_QUALITY_WITHOUT_ROLLBACK;
        DBG << RED << "Rollback of move sequence violated balance constraint ( Moved Nodes ="
            << sequence.moves.size()
            << ", Expected Improvement =" << sequence.expected_improvement
            << ", Real Improvement =" << improvement
            << ", Search ID =" << search_id << ")" << END;
      }
    } else {
      ++_stats.num_improvements;
      _stats.correct_expected_improvement += (improvement == sequence.expected_improvement);
      sequence.state = MoveSequenceState::SUCCESS;
      DBG << GREEN << "SUCCESS -"
          << "Moved Nodes =" << sequence.moves.size()
          << ", Expected Improvement =" << sequence.expected_improvement
          << ", Real Improvement =" << improvement
          << ", Search ID =" << search_id << END;
    }
  } else {
    ++_stats.failed_updates_due_to_balance_constraint;
    sequence.state = MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT;
    DBG << RED << "Move sequence violated balance constraint ( Moved Nodes ="
        << sequence.moves.size()
        << ", Expected Improvement =" << sequence.expected_improvement
        << ", Search ID =" << search_id << ")" << END;
  }

  if ( _context.refinement.advanced.use_overlapping_searches ) {
    _apply_moves_lock.unlock();
  }

  if ( sequence.state == MoveSequenceState::SUCCESS && improvement > 0 ) {
    addCutHyperedgesToQuotientGraph(_quotient_graph, new_cut_hes);
    _stats.total_improvement += improvement;
  } else if ( sequence.state == MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT &&
              sequence.expected_improvement > 0 && num_retries > 0 &&
              tryToFixMoveSequenceAfterBalanceViolation(sequence, update_res) ) {
    --_stats.failed_updates_due_to_balance_constraint;
    improvement = applyMoves(search_id, sequence, num_retries - 1);
  }

  return improvement;
}

AdvancedRefinementScheduler::PartWeightUpdateResult AdvancedRefinementScheduler::partWeightUpdate(
  const vec<HypernodeWeight>& part_weight_deltas, const bool rollback) {
  const HypernodeWeight multiplier = rollback ? -1 : 1;
  PartWeightUpdateResult res;
  _part_weights_lock.lock();
  PartitionID i = 0;
  for ( ; i < _context.partition.k; ++i ) {
    if ( _part_weights[i] + multiplier * part_weight_deltas[i] > _max_part_weights[i] ) {
      DBG << "Move sequence violated balance constraint of block" << i
          << "(Max =" << _max_part_weights[i]
          << ", Actual =" << (_part_weights[i] + multiplier * part_weight_deltas[i]) << ")";
      res.is_balanced = false;
      res.overloaded_block = i;
      res.overload_weight = ( _part_weights[i] + multiplier *
        part_weight_deltas[i] ) - _max_part_weights[i];
      // Move Sequence Violates Balance Constraint => Rollback
      --i;
      for ( ; i >= 0; --i ) {
        _part_weights[i] -= multiplier * part_weight_deltas[i];
      }
      break;
    }
    _part_weights[i] += multiplier * part_weight_deltas[i];
  }
  _part_weights_lock.unlock();
  return res;
}

bool AdvancedRefinementScheduler::tryToFixMoveSequenceAfterBalanceViolation(
  MoveSequence& sequence, const PartWeightUpdateResult& update_res) {
  // Sort nodes in decreasing order of their node degrees and remove
  // nodes from sequence with low degree. Assumption is that low degree
  // vertices have small gains
  std::sort(sequence.moves.begin(), sequence.moves.end(),
    [&](const Move& lhs, const Move& rhs) {
      const bool lhs_move_to_overloaded_block = lhs.to == update_res.overloaded_block;
      const bool rhs_move_to_overloaded_block = rhs.to == update_res.overloaded_block;
      return  lhs_move_to_overloaded_block < rhs_move_to_overloaded_block ||
        ( lhs_move_to_overloaded_block == rhs_move_to_overloaded_block &&
          _phg->nodeDegree(lhs.node) > _phg->nodeDegree(rhs.node) );
    });

  HypernodeWeight overload_weight = update_res.overload_weight;
  while ( !sequence.moves.empty() && overload_weight > 0 ) {
    const Move& move = sequence.moves.back();
    if ( move.to == update_res.overloaded_block ) {
      if ( _phg->partID(move.node) != update_res.overloaded_block ) {
        overload_weight -= _phg->nodeWeight(move.node);
      }
      sequence.moves.pop_back();
    } else {
      break;
    }
  }
  return overload_weight <= 0;
}

}