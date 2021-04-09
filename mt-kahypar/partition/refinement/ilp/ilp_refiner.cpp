/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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


#include "mt-kahypar/partition/refinement/ilp/ilp_refiner.h"

#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {

MoveSequence ILPRefiner::refineImpl(const PartitionedHypergraph& phg,
                                    const AdvancedProblem& problem) {

  MoveSequence sequence { {}, 0 };
  if ( estimatedNumberOfNonZeros(problem.stats) >
       _context.refinement.advanced.ilp.min_non_zeros ) {
    // Construct ILP Model
    utils::Timer::instance().start_timer("construct_ilp", "Construct ILP", true);
    GRBModel ilp = _model.construct(phg, problem);
    utils::Timer::instance().stop_timer("construct_ilp");

    try {
      ilp.set(GRB_DoubleParam_TimeLimit, _context.refinement.advanced.ilp.time_limit);
      ilp.set(GRB_IntParam_LogToConsole,
        _context.partition.verbose_output && debug &&
        _context.shared_memory.num_threads == _context.refinement.advanced.num_threads_per_search);
      ilp.set(GRB_IntParam_Threads, _num_threads);
      ilp.set("MIPFocus", "1");
      ilp.set("Symmetry", "2");

      utils::Timer::instance().start_timer("solve_ilp", "Solve ILP", true);
      const HyperedgeWeight objective_before = _model.getInitialObjective();
      ilp.optimize();
      const HyperedgeWeight objective_after = ilp.get(GRB_DoubleAttr_ObjVal);
      utils::Timer::instance().stop_timer("solve_ilp");

      utils::Timer::instance().start_timer("retrieve_solution", "Retrieve Solution", true);
      int status = ilp.get(GRB_IntAttr_Status);
      if ( status == GRB_OPTIMAL || status == GRB_TIME_LIMIT ) {
        sequence.expected_improvement = objective_before - objective_after;
        if ( _context.refinement.advanced.ilp.apply_zero_gain_moves ||
             sequence.expected_improvement > 0 ) {
          for ( const HypernodeID& hn : problem.nodes ) {
            const PartitionID from = phg.partID(hn);
            const PartitionID to = _model.partID(hn, problem);
            ASSERT(to != kInvalidPartition);
            if ( from != to ) {
              sequence.moves.emplace_back(Move { from, to, hn, kInvalidGain });
            }
          }
        }
      }
      utils::Timer::instance().stop_timer("retrieve_solution");

      if ( status == GRB_OPTIMAL ) {
        utils::Stats::instance().update_stat("ilp_solved_optimal", 1);
      } else if ( status == GRB_TIME_LIMIT ) {
        utils::Stats::instance().update_stat("ilp_time_limit", 1);
      }
    } catch(GRBException e) {
      DBG << RED << "Error code = " << e.getErrorCode() << "Message =" << e.getMessage() << END;
    } catch(...) {
      ERROR("Exception during optimization");
    }
  }

  if ( sequence.moves.size() > 0 ) {
    DBG << V(problem.stats.numNodes())
        << V(problem.stats.numEdges())
        << V(problem.stats.numPins())
        << V(problem.stats.numContainedBlocks())
        << V(sequence.moves.size())
        << V(sequence.expected_improvement);
  }

  return sequence;
}

} // namespace mt_kahypar