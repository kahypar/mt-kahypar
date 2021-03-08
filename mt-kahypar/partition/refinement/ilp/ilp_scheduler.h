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
#include "mt-kahypar/partition/refinement/ilp/ilp_solver.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

class ILPScheduler {

  static constexpr bool debug = false;

  struct VertexGain {
    HypernodeID hn;
    Gain gain;
    bool is_border_vertex;
  };

 public:
  explicit ILPScheduler(PartitionedHypergraph& phg,
                        const Context& context) :
    _phg(phg),
    _context(context),
    _env(),
    _solver(phg, context, _env),
    _num_vertices(0),
    _num_hyperedges(0),
    _num_pins(0),
    _k(0),
    _gains(),
    _visited_hns(),
    _visited_hes() { }

  bool refine();

 private:
  void computeGains() {
    if ( _context.partition.objective == kahypar::Objective::km1 ) {
      computeGainsForConnectivityMetric();
    } else if ( _context.partition.objective == kahypar::Objective::cut ) {
      computeGainsForCutMetric();
    }
  }

  void computeGainsForConnectivityMetric();

  void computeGainsForCutMetric();

  size_t estimatedNumberOfNonZeros() {
    if ( _context.partition.objective == kahypar::Objective::km1 ) {
      return _k * ( 2 * _num_vertices + 2 * _num_pins );
    } else if ( _context.partition.objective == kahypar::Objective::cut ) {
      return _k * ( 2 * _num_vertices + 3 * ( _num_pins - _num_hyperedges ) );
    }
    return std::numeric_limits<size_t>::max();
  }

  // ! Returns the number of vertices initially inserted into the bfs queue
  size_t bfs(vec<HypernodeID>& nodes,
             const size_t gains_start_idx,
             const size_t gains_end_idx,
             const Gain min_gain,
             const int max_distance);

  void bfs(vec<HypernodeID>& nodes,
           const HypernodeID start_hn,
           const int max_distance);

  PartitionedHypergraph& _phg;
  const Context& _context;
  GRBEnv _env;
  ILPSolver _solver;

  size_t _num_vertices;
  size_t _num_hyperedges;
  size_t _num_pins;
  PartitionID _k;
  vec<VertexGain> _gains;

  vec<bool> _visited_hns;
  vec<bool> _marked_hns;
  kahypar::ds::FastResetFlagArray<> _visited_hes;
};

} // namespace mt_kahypar