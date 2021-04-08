/*******************************************************************************
 * This file is part of MT-KaHyPar.
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

#pragma once

#include "gurobi_c++.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/advanced/problem_stats.h"

namespace mt_kahypar {

class ILPModel {

  static constexpr bool debug = false;

 public:
  explicit ILPModel(const Context& context) :
    _context(context),
    _env(),
    _num_nodes(0),
    _num_edges(0),
    _k(0),
    _initial_objective(0),
    _variables(),
    _contains_he_var(),
    _nodes_to_ilp(),
    _edges_to_ilp() {
    // Sanity check
    if ( _context.partition.objective != kahypar::Objective::km1 &&
         _context.partition.objective != kahypar::Objective::cut ) {
      ERROR("ILP Model is not able to optimize" << _context.partition.objective << "metric");
    }
  }

  void initialize(const AdvancedProblem& problem);

  GRBModel construct(const PartitionedHypergraph& phg,
                     const AdvancedProblem& problem);

  HyperedgeWeight getInitialObjective() const {
    return _initial_objective;
  }

  PartitionID partID(const HypernodeID hn,
                     const AdvancedProblem& problem);

 private:
  // ####################### Variable #######################

  void addVariablesToModel(GRBModel& model,
                           const PartitionedHypergraph& phg,
                           const AdvancedProblem& problem);

  // ####################### Objective Function #######################

  GRBLinExpr getConnectivityMetric(const PartitionedHypergraph& phg,
                                   const AdvancedProblem& problem);

  void addObjectiveFunction(GRBModel& model,
                            const PartitionedHypergraph& phg,
                            const AdvancedProblem& problem);

  // ####################### Constraints #######################

  void restrictVerticesToOneBlock(GRBModel& model,
                                  const AdvancedProblem& problem);

  void balanceConstraint(GRBModel& model,
                         const PartitionedHypergraph& phg,
                         const AdvancedProblem& problem);

  void modelHyperedgeConnectivity(GRBModel& model,
                                  const PartitionedHypergraph& phg,
                                  const AdvancedProblem& problem);

  // ####################### Helper Functions #######################

  std::string vertex_var_desc(const HypernodeID hn, const PartitionID k) {
    return "x_{" + std::to_string(hn) + "," + std::to_string(k) + "}";
  }

  std::string hyperedge_var_desc(const HyperedgeID he, const PartitionID k) {
    return "y_{" + std::to_string(he) + "," + std::to_string(k) + "}";
  }

  size_t vertex_offset(const HypernodeID hn, const PartitionID k) {
    return hn * _k + k;
  }

  size_t hyperedge_offset(const HyperedgeID he, const PartitionID k) {
    return _num_nodes * _k + he * _k + k;
  }

  const Context& _context;
  GRBEnv _env;

  HypernodeID _num_nodes;
  HyperedgeID _num_edges;
  PartitionID _k;
  HyperedgeWeight _initial_objective;
  vec<GRBVar> _variables;

  vec<bool> _contains_he_var;
  ds::DynamicSparseMap<HypernodeID, HypernodeID> _nodes_to_ilp;
  ds::DynamicSparseMap<HyperedgeID, HyperedgeID> _edges_to_ilp;
};

} // namespace mt_kahypar