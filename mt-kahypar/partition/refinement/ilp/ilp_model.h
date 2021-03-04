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
  explicit ILPModel(ILPHypergraph& hg,
                    const Context& context,
                    const GRBEnv& env) :
    _hg(hg),
    _context(context),
    _is_constructed(false),
    _is_solved(false),
    _model(env),
    _objective(0),
    _variables(hg.numNodes() * hg.k() + hg.numEdges() * hg.k()),
    _contains_variable(hg.numNodes() * hg.k() + hg.numEdges() * hg.k(), true),
    _num_unremovable_blocks(hg.numEdges(), 0) {
    _model.set(GRB_IntParam_LogToConsole, debug);
  }

  void construct();

  void solve() {
    ASSERT(_is_constructed);
    try {
      _model.optimize();
      _objective = _model.get(GRB_DoubleAttr_ObjVal);
      _is_solved = true;
    } catch(GRBException e) {
      ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
    } catch(...) {
      ERROR("Exception during optimization");
    }
  }

  // ! Returns the connectivity metric
  HyperedgeWeight getObjective() const {
    ASSERT(_is_constructed);
    return _objective;
  }

  // ! Returns the block of vertex hn to which the ILP Optimizer
  // ! assigns it to
  PartitionID partID(const HypernodeID hn) {
    ASSERT(_is_solved);
    try {
      const HypernodeID u = _hg.mapToILPHypergraph(hn);
      for ( PartitionID i = 0; i < _hg.k(); ++i ) {
        if ( _variables[vertex_offset(u,i)].get(GRB_DoubleAttr_X) > 0 ) {
          return _hg.toOriginalBlock(i);
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

 private:

  void addVariablesToModel();

  void addObjectiveFunction();

  void restrictVerticesToOneBlock();

  void balanceConstraint();

  void modelHyperedgeConnectivity();

  std::string vertex_var_desc(const HypernodeID hn, const PartitionID k) {
    return "x_{" + std::to_string(hn) + "," + std::to_string(k) + "}";
  }

  std::string hyperedge_var_desc(const HyperedgeID he, const PartitionID k) {
    return "y_{" + std::to_string(he) + "," + std::to_string(k) + "}";
  }

  size_t vertex_offset(const HypernodeID hn, const PartitionID k) {
    return hn * _hg.k() + k;
  }

  size_t hyperedge_offset(const HyperedgeID he, const PartitionID k) {
    return _hg.numNodes() * _hg.k() + he * _hg.k() + k;
  }

  ILPHypergraph& _hg;
  const Context& _context;

  // ! True, if ILP is constructed
  bool _is_constructed;
  // ! True, if ILP is solved
  bool _is_solved;
  // ! Gurobi model
  GRBModel _model;
  // ! Connectivity metric
  HyperedgeWeight _objective;
  // ! Gurobi variables
  vec<GRBVar> _variables;
  // ! _contains_variable[i] == true, if variable i is contained in the Gurobi model
  vec<bool> _contains_variable;
  // ! Indicates for each hyperedge the number of blocks which we can not
  // ! remove with our ILP formulation, because they contain pins which are
  // ! not contained in our subhypergraph.
  vec<PartitionID> _num_unremovable_blocks;
};

} // namespace mt_kahypar