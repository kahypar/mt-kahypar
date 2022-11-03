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

#include "mt-kahypar/partition/refinement/flows/scheduler.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {

void FlowRefinementScheduler::RefinementStats::update_global_stats() {
  _stats.update_stat("num_flow_refinements",
    num_refinements.load(std::memory_order_relaxed));
  _stats.update_stat("num_flow_improvement",
    num_improvements.load(std::memory_order_relaxed));
  _stats.update_stat("num_time_limits",
    num_time_limits.load(std::memory_order_relaxed));
  _stats.update_stat("correct_expected_improvement",
    correct_expected_improvement.load(std::memory_order_relaxed));
  _stats.update_stat("zero_gain_improvement",
    zero_gain_improvement.load(std::memory_order_relaxed));
  _stats.update_stat("failed_updates_due_to_conflicting_moves",
    failed_updates_due_to_conflicting_moves.load(std::memory_order_relaxed));
  _stats.update_stat("failed_updates_due_to_conflicting_moves_without_rollback",
    failed_updates_due_to_conflicting_moves_without_rollback.load(std::memory_order_relaxed));
  _stats.update_stat("failed_updates_due_to_balance_constraint",
    failed_updates_due_to_balance_constraint.load(std::memory_order_relaxed));
  _stats.update_stat("total_flow_refinement_improvement",
    total_improvement.load(std::memory_order_relaxed));
}

namespace {

  static constexpr size_t PROGRESS_BAR_SIZE = 50;

  template<typename F>
  std::string progress_bar(const size_t value, const size_t max, const F& f) {
    const double percentage = static_cast<double>(value) / std::max(max,1UL);
    const size_t ticks = PROGRESS_BAR_SIZE * percentage;
    std::stringstream pbar_str;
    pbar_str << "|"
             << f(percentage) << std::string(ticks, '|') << END
             << std::string(PROGRESS_BAR_SIZE - ticks, ' ')
             << "| " << std::setprecision(2) << (100.0 * percentage) << "% (" << value << ")";
    return pbar_str.str();
  }
}

inline std::ostream & operator<< (std::ostream& str, const FlowRefinementScheduler::RefinementStats& stats) {
  str << "\n";
  str << "Total Improvement                   = " << stats.total_improvement << "\n";
  str << "Number of Flow-Based Refinements    = " << stats.num_refinements << "\n";
  str << "+ No Improvements                   = "
      << progress_bar(stats.num_refinements - stats.num_improvements, stats.num_refinements,
          [&](const double percentage) { return percentage > 0.9 ? RED : percentage > 0.75 ? YELLOW : GREEN; }) << "\n";
  str << "+ Number of Improvements            = "
      << progress_bar(stats.num_improvements, stats.num_refinements,
          [&](const double percentage) { return percentage < 0.05 ? RED : percentage < 0.15 ? YELLOW : GREEN; }) << "\n";
  str << "  + Correct Expected Improvements   = "
      << progress_bar(stats.correct_expected_improvement, stats.num_improvements,
          [&](const double percentage) { return percentage > 0.9 ? GREEN : percentage > 0.75 ? YELLOW : RED; }) << "\n";
  str << "  + Incorrect Expected Improvements = "
      << progress_bar(stats.num_improvements - stats.correct_expected_improvement, stats.num_improvements,
          [&](const double percentage) { return percentage < 0.1 ? GREEN : percentage < 0.25 ? YELLOW : RED; }) << "\n";
  str << "  + Zero-Gain Improvements          = "
      << progress_bar(stats.zero_gain_improvement, stats.num_improvements,
          [&](const double) { return WHITE; }) << "\n";
  str << "+ Failed due to Balance Constraint  = "
      << progress_bar(stats.failed_updates_due_to_balance_constraint, stats.num_refinements,
          [&](const double percentage) { return percentage < 0.01 ? GREEN : percentage < 0.05 ? YELLOW : RED; }) << "\n";
  str << "+ Failed due to Conflicting Moves   = "
      << progress_bar(stats.failed_updates_due_to_conflicting_moves, stats.num_refinements,
          [&](const double percentage) { return percentage < 0.01 ? GREEN : percentage < 0.05 ? YELLOW : RED; }) << "\n";
  str << "+ Time Limits                       = "
      << progress_bar(stats.num_time_limits, stats.num_refinements,
          [&](const double percentage) { return percentage < 0.0025 ? GREEN : percentage < 0.01 ? YELLOW : RED; }) << "\n";
  str << "---------------------------------------------------------------";
  return str;
}


bool FlowRefinementScheduler::refineImpl(
                PartitionedHypergraph& phg,
                const parallel::scalable_vector<HypernodeID>&,
                Metrics& best_metrics,
                const double)  {
  phg.resetMoveState();
  ASSERT(_phg == &phg);
  _quotient_graph.setObjective(best_metrics.getMetric(
    Mode::direct, _context.partition.objective));

  std::atomic<HyperedgeWeight> overall_delta(0);
  utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  tbb::parallel_for(0UL, _refiner.numAvailableRefiner(), [&](const size_t i) {
    while ( i < std::max(1UL, static_cast<size_t>(
        std::ceil(_context.refinement.flows.parallel_searches_multiplier *
            _quotient_graph.numActiveBlockPairs()))) ) {
      SearchID search_id = _quotient_graph.requestNewSearch(_refiner);
      if ( search_id != QuotientGraph::INVALID_SEARCH_ID ) {
        DBG << "Start search" << search_id
            << "( Blocks =" << blocksOfSearch(search_id)
            << ", Refiner =" << i << ")";
        timer.start_timer("region_growing", "Grow Region", true);
        const Subhypergraph sub_hg =
          _constructor.construct(search_id, _quotient_graph, phg);
        _quotient_graph.finalizeConstruction(search_id);
        timer.stop_timer("region_growing");

        HyperedgeWeight delta = 0;
        bool improved_solution = false;
        if ( sub_hg.numNodes() > 0 ) {
          ++_stats.num_refinements;
          MoveSequence sequence = _refiner.refine(search_id, phg, sub_hg);

          if ( !sequence.moves.empty() ) {
            timer.start_timer("apply_moves", "Apply Moves", true);
            delta = applyMoves(search_id, sequence);
            overall_delta -= delta;
            improved_solution = sequence.state == MoveSequenceState::SUCCESS && delta > 0;
            timer.stop_timer("apply_moves");
          } else if ( sequence.state == MoveSequenceState::TIME_LIMIT ) {
            ++_stats.num_time_limits;
            DBG << RED << "Search" << search_id << "reaches the time limit ( Time Limit ="
                << _refiner.timeLimit() << "s )" << END;
          }
        }
        _quotient_graph.finalizeSearch(search_id, improved_solution ? delta : 0);
        _refiner.finalizeSearch(search_id);
        DBG << "End search" << search_id
            << "( Blocks =" << blocksOfSearch(search_id)
            << ", Refiner =" << i
            << ", Running Time =" << _refiner.runningTime(search_id) << ")";
      } else {
        break;
      }
    }
    _refiner.terminateRefiner();
    DBG << RED << "Refiner" << i << "terminates!" << END;
  });

  DBG << _stats;

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
    Mode::direct, _context.partition.objective);
  HEAVY_REFINEMENT_ASSERT(current_metric + overall_delta ==
                          metrics::objective(phg, _context.partition.objective),
                          V(current_metric) << V(overall_delta) <<
                          V(metrics::objective(phg, _context.partition.objective)));
  best_metrics.updateMetric(current_metric + overall_delta,
    Mode::direct, _context.partition.objective);
  best_metrics.imbalance = metrics::imbalance(phg, _context);
  _stats.update_global_stats();

  // Update Gain Cache
  if ( ( _context.partition.paradigm == Paradigm::nlevel ||
         _context.refinement.refine_until_no_improvement ) &&
         phg.isGainCacheInitialized() ) {
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( _was_moved[hn] ) {
        phg.recomputeMoveFromPenalty(hn);
        _was_moved[hn] = uint8_t(false);
      }
    });
  }

  HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
  _phg = nullptr;
  return overall_delta.load(std::memory_order_relaxed) < 0;
}

void FlowRefinementScheduler::initializeImpl(PartitionedHypergraph& phg)  {
  _phg = &phg;

  // Initialize Part Weights
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    _part_weights[i] = phg.partWeight(i);
    _max_part_weights[i] = std::max(
      phg.partWeight(i), _context.partition.max_part_weights[i]);
  }

  _stats.reset();
  utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  timer.start_timer("initialize_quotient_graph", "Initialize Quotient Graph");
  _quotient_graph.initialize(phg);
  timer.stop_timer("initialize_quotient_graph");

  const size_t max_parallism = _context.refinement.flows.num_parallel_searches;
  DBG << "Initial Active Block Pairs =" << _quotient_graph.numActiveBlockPairs()
      << ", Initial Num Threads =" << max_parallism;
  _refiner.initialize(max_parallism);
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
                    const bool gain_cache_update) {
  bool success = false;
  if ( gain_cache_update && phg.isGainCacheInitialized()) {
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
                       const bool gain_cache_update,
                       vec<uint8_t>& was_moved,
                       vec<NewCutHyperedge>& new_cut_hes) {
  for ( const Move& move : sequence.moves ) {
    ASSERT(move.from == phg.partID(move.node));
    if ( move.from != move.to ) {
      changeNodePart(phg, move.node, move.from, move.to, objective_delta, gain_cache_update);
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
                        const bool gain_cache_update) {
  for ( const Move& move : sequence.moves ) {
    if ( move.from != move.to ) {
      ASSERT(phg.partID(move.node) == move.to);
      changeNodePart(phg, move.node, move.to, move.from, objective_delta, gain_cache_update);
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

HyperedgeWeight FlowRefinementScheduler::applyMoves(const SearchID search_id,
                                                        MoveSequence& sequence) {
  unused(search_id);
  ASSERT(_phg);

  // TODO: currently we lock the applyMoves method
  // => find something smarter here
  _apply_moves_lock.lock();

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
    if ( _context.partition.objective == Objective::km1 ) {
      improvement -= km1Delta(he, edge_weight, edge_size,
        pin_count_in_from_part_after, pin_count_in_to_part_after);
    } else if ( _context.partition.objective == Objective::cut ) {
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
    const bool gain_cache_update =
      _context.partition.paradigm == Paradigm::nlevel ||
      _context.refinement.refine_until_no_improvement;
    applyMoveSequence(*_phg, sequence, delta_func, gain_cache_update, _was_moved, new_cut_hes);

    if ( improvement < 0 ) {
      update_res = partWeightUpdate(part_weight_deltas, true);
      if ( update_res.is_balanced ) {
        // Move sequence worsen solution quality => Rollback
        DBG << RED << "Move sequence worsen solution quality ("
            << "Expected Improvement =" << sequence.expected_improvement
            << ", Real Improvement =" << improvement
            << ", Search ID =" << search_id << ")" << END;
        revertMoveSequence(*_phg, sequence, delta_func, gain_cache_update);
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
      _stats.zero_gain_improvement += (improvement == 0);
      sequence.state = MoveSequenceState::SUCCESS;
      DBG << ( improvement > 0 ? GREEN : "" ) << "SUCCESS -"
          << "Moved Nodes =" << sequence.moves.size()
          << ", Expected Improvement =" << sequence.expected_improvement
          << ", Real Improvement =" << improvement
          << ", Search ID =" << search_id << ( improvement > 0 ? END : "" );
    }
  } else {
    ++_stats.failed_updates_due_to_balance_constraint;
    sequence.state = MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT;
    DBG << RED << "Move sequence violated balance constraint ( Moved Nodes ="
        << sequence.moves.size()
        << ", Expected Improvement =" << sequence.expected_improvement
        << ", Search ID =" << search_id << ")" << END;
  }

  _apply_moves_lock.unlock();

  if ( sequence.state == MoveSequenceState::SUCCESS && improvement > 0 ) {
    addCutHyperedgesToQuotientGraph(_quotient_graph, new_cut_hes);
    _stats.total_improvement += improvement;
  }

  return improvement;
}

FlowRefinementScheduler::PartWeightUpdateResult FlowRefinementScheduler::partWeightUpdate(
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

}