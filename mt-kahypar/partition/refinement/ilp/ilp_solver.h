/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include "gurobi_c++.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_hypergraph.h"
#include "mt-kahypar/partition/refinement/ilp/ilp_model.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

class ILPSolver {

  static constexpr bool debug = false;

  struct RollbackElement {
    HypernodeID hn;
    PartitionID from;
    PartitionID to;
  };

 public:
  explicit ILPSolver(PartitionedHypergraph& phg,
                     const Context& context,
                     const GRBEnv& env) :
    _phg(phg),
    _context(context),
    _ilp_hg(phg),
    _model(context, env),
    _rollback_cache() { }

  /**
   * Solves an ILP on the subhypergraph defined by nodes and
   * returns the improvement with respect to the objective function
   * found by the ILP.
   */
  HyperedgeWeight solve(const vec<HypernodeID>& nodes) {
    _ilp_hg.reset();
    _model.reset();
    _rollback_cache.clear();

    // Setup ILP
    utils::Timer::instance().start_timer("construct_ilp", "Construct ILP", true);
    utils::Timer::instance().start_timer("construct_ilp_hypergraph", "Construct ILP Hypergraph", true);
    _ilp_hg.initialize(nodes);
    utils::Timer::instance().stop_timer("construct_ilp_hypergraph");
    utils::Timer::instance().start_timer("construct_ilp_problem", "Construct ILP Problem", true);
    _model.construct(_ilp_hg);
    utils::Timer::instance().stop_timer("construct_ilp_problem");
    utils::Timer::instance().stop_timer("construct_ilp");

    // Solve ILP
    utils::Timer::instance().start_timer("solve_ilp", "Solve ILP", true);
    int status = _model.solve();
    utils::Timer::instance().stop_timer("solve_ilp");

    HyperedgeWeight delta = 0;
    if ( status == GRB_OPTIMAL || status == GRB_TIME_LIMIT ) {
      DBG << "ILP solver improved objective function from"
          << _model.getInitialObjective()
          << "to" << _model.getOptimizedObjective()
          << "( Delta =" << (_model.getInitialObjective() - _model.getOptimizedObjective()) << ")";
      // Apply solution
      utils::Timer::instance().start_timer("move_vertices", "Move Vertices", true);
      delta = applyMoves(nodes);
      if ( delta < 0 ) {
        // Applying the moves found by the ILP worsen solution quality.
        // This can happen if several ILP's apply their moves concurrently.
        // => revert moves
        DBG << "Rollback: Move decreased objective function by" << delta;
        delta += rollback();
      } else {
        DBG << "Applying moves improved objective function by" << delta;
      }
      utils::Timer::instance().stop_timer("move_vertices");
    }

    return delta;
  }

 private:
  HyperedgeWeight applyMoves(const vec<HypernodeID>& nodes) {
    HyperedgeWeight delta = 0;
    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID edge_size,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      if ( _context.partition.objective == kahypar::Objective::km1 ) {
        delta += km1Delta(he, edge_weight, edge_size,
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      } else if ( _context.partition.objective == kahypar::Objective::cut ) {
        delta += cutDelta(he, edge_weight, edge_size,
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      }
    };
    for ( const HypernodeID& hn : nodes ) {
      const PartitionID from = _phg.partID(hn);
      const PartitionID to = _model.partID(hn);
      if ( from != to && _phg.changeNodePart(hn, from, to, delta_func) ) {
        DBG << "Move vertex" << hn << "from block" << from << "to" << to;
        _rollback_cache.push_back(RollbackElement { hn, from, to });
      }
    }
    return -delta;
  }

  HyperedgeWeight rollback() {
    HyperedgeWeight delta = 0;
    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID edge_size,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      if ( _context.partition.objective == kahypar::Objective::km1 ) {
        delta += km1Delta(he, edge_weight, edge_size,
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      } else if ( _context.partition.objective == kahypar::Objective::cut ) {
        delta += cutDelta(he, edge_weight, edge_size,
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      }
    };
    for ( const RollbackElement& element : _rollback_cache ) {
      ASSERT(_phg.partID(element.hn) == element.to);
      _phg.changeNodePart(element.hn, element.to, element.from, delta_func);
    }
    return -delta;
  }

  PartitionedHypergraph& _phg;
  const Context& _context;
  ILPHypergraph _ilp_hg;
  ILPModel _model;
  vec<RollbackElement> _rollback_cache;
};

} // namespace mt_kahypar