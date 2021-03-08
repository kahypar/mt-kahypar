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

namespace mt_kahypar {

class ILPModel {

  static constexpr bool debug = false;

 public:
  explicit ILPModel(const Context& context,
                    const GRBEnv& env) :
    _hg(nullptr),
    _context(context),
    _is_constructed(false),
    _is_solved(false),
    _model(env),
    _initial_objective(0),
    _optimized_objective(0),
    _variables(),
    _contains_variable(),
    _unremovable_block() {
    // Sanity check
    if ( _context.partition.objective != kahypar::Objective::km1 &&
         _context.partition.objective != kahypar::Objective::cut ) {
      ERROR("ILP Model is not able to optimize" << _context.partition.objective << "metric");
    }
    _model.set(GRB_DoubleParam_TimeLimit, _context.refinement.ilp.time_limit);
  }

  void construct(ILPHypergraph& hg);

  int solve(const bool supress_output = false) {
    _model.set(GRB_IntParam_LogToConsole,
      _context.partition.verbose_output && !supress_output);
    ASSERT(_is_constructed);
    int status = -1;
    try {
      _model.optimize();
      status = _model.get(GRB_IntAttr_Status);
      if ( status == GRB_OPTIMAL ) {
        _optimized_objective = _model.get(GRB_DoubleAttr_ObjVal);
        _is_solved = true;
      }
    } catch(GRBException e) {
      ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
    } catch(...) {
      ERROR("Exception during optimization");
    }
    return status;
  }

  // ! Returns the objective given by the input solution
  HyperedgeWeight getInitialObjective() const {
    ASSERT(_is_constructed);
    return _initial_objective;
  }

  // ! Returns the objective given by the input solution
  HyperedgeWeight getOptimizedObjective() const {
    ASSERT(_is_solved);
    return _optimized_objective;
  }

  // ! Returns the block of vertex hn to which the ILP Optimizer
  // ! assigns it to
  PartitionID partID(const HypernodeID hn) {
    ASSERT(_hg);
    ASSERT(_is_solved);
    try {
      const HypernodeID u = _hg->mapToILPHypergraph(hn);
      for ( PartitionID i = 0; i < _hg->k(); ++i ) {
        if ( _variables[vertex_offset(u,i)].get(GRB_DoubleAttr_X) > 0 ) {
          return _hg->toOriginalBlock(i);
        }
      }
    } catch(GRBException e) {
      ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
    } catch(...) {
      ERROR("Exception during optimization");
    }
    ERROR("Invalid solution: Vertex" << hn << "must be assigned to one block");
    return kInvalidPartition;
  }

  void reset() {
    _hg = nullptr;
    _is_constructed = false;
    _is_solved = false;
    _model.reset();
    _initial_objective = 0;
    _optimized_objective = 0;
    _variables.clear();
    _contains_variable.clear();
    _unremovable_block.clear();
  }

 private:

  // ####################### Variable #######################

  void addVariablesToModel(ILPHypergraph& hg);

  PartitionID determineNumberOfUnremovableBlocks(ILPHypergraph& hg, const HyperedgeID he);

  void addHyperedgeVariablesToModelForConnectivityMetric(ILPHypergraph& hg);

  void addHyperedgeVariablesToModelForCutMetric(ILPHypergraph& hg);

  // ####################### Objective Function #######################

  void addObjectiveFunction(ILPHypergraph& hg);

  void addConnectivityMetric(ILPHypergraph& hg);

  void addCutMetric(ILPHypergraph& hg);

  // ####################### Constraints #######################

  void restrictVerticesToOneBlock(ILPHypergraph& hg);

  void balanceConstraint(ILPHypergraph& hg);

  void modelHyperedgeConstraint(ILPHypergraph& hg);

  void modelHyperedgeConnectivity(ILPHypergraph& hg);

  void modelHyperedgeCut(ILPHypergraph& hg);

  // ####################### Helper Functions #######################

  std::string vertex_var_desc(const HypernodeID hn, const PartitionID k) {
    ASSERT(_hg);
    return "x_{" + std::to_string(hn) + "," + std::to_string(k) + "}";
  }

  std::string hyperedge_var_desc(const HyperedgeID he, const PartitionID k) {
    ASSERT(_hg);
    return "y_{" + std::to_string(he) + "," + std::to_string(k) + "}";
  }

  size_t vertex_offset(const HypernodeID hn, const PartitionID k) {
    ASSERT(_hg);
    return hn * _hg->k() + k;
  }

  size_t hyperedge_offset(const HyperedgeID he, const PartitionID k) {
    ASSERT(_hg);
    return _hg->numNodes() * _hg->k() + he * _hg->k() + k;
  }

  ILPHypergraph* _hg;
  const Context& _context;

  // ! True, if ILP is constructed
  bool _is_constructed;
  // ! True, if ILP is solved
  bool _is_solved;
  // ! Gurobi model
  GRBModel _model;
  // ! Objective given by the input solution
  HyperedgeWeight _initial_objective;
  // ! Objective after ILP problem was solved
  HyperedgeWeight _optimized_objective;
  // ! Gurobi variables
  vec<GRBVar> _variables;
  // ! _contains_variable[i] == true, if variable i is contained in the Gurobi model
  vec<bool> _contains_variable;
  // ! Only for cut metric, contains the block of the hyperedge which we can not remove
  // ! with our ILP problem
  vec<PartitionID> _unremovable_block;
};

} // namespace mt_kahypar