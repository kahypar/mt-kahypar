/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019, 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"

#include "ilp_model.h"

namespace mt_kahypar {

/**
 * Constructs the following ILP Problem:
 *
 * Binary Decision Variables:
 * x_{v,i} -> If 1, then vertex v is assigned to block i
 * y_{e,i} -> If 1, then hyperedge e contains at least one pin in block i
 *
 * Constraints:
 *
 * Each block vertex v is assigned exactly to one block:
 * \forall v \in V: x_{v,1} + ... + x_{v,k} = 1  (1)
 *
 * Each block satisfies the balance constraint:
 * \forall k: \sum_{v \in V} x_{v,k} * c(v) <= L_max (2)
 *
 * y_{e,i} is equal to 1, if at least one pin is assigned to block i
 * \forall e \in E: \forall v \in e: \forall k: y_{e,k} >= x_{e,k} (3)
 *
 * Objective Function:
 * minimize \sum_{e \in E} (\lambda(e) - 1) * \omega(e) with \lambda(e) = \sum_{i = 1}^k y_{e,i}
 *
 * Optimizations (TODOs):
 *  - remove super vertices that represents vertices of a block not contained in the ILP
 *    o For a supervertex v representing block i, we adapt the balance constraint to ... <= L_max - c(v)
 *    o If a supervertex v representing block i is contained in a hyperedge e, we can remove
 *      variable y_{e,i} because each solution will have y_{e,i} = 1 (adapt objective function accordingly).
 *  - add balance constraint lazily, e.g. if optimal solution is found but did not satisfy the balance constraint
 */
void ILPModel::construct() {
  try {
    addVariablesToModel();
    addObjectiveFunction();
    restrictVerticesToOneBlock();
    balanceConstraint();
    modelHyperedgeConnectivity();
    _is_constructed = true;
  } catch(GRBException e) {
    ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
  } catch(...) {
    ERROR("Exception during optimization");
  }
}

void ILPModel::addVariablesToModel() {
  for ( const HypernodeID& hn : _hg.nodes() ) {
    for ( PartitionID i = 0; i < _hg.k(); ++i ) {
      _variables[vertex_offset(hn, i)] =
        _model.addVar(0.0, 1.0, 0.0, GRB_BINARY, vertex_var_desc(hn, i));
      _variables[vertex_offset(hn, i)].set(
        GRB_DoubleAttr_Start, (_hg.partID(hn) == i));
    }
  }

  for ( const HypernodeID& he : _hg.edges() ) {
    for ( PartitionID i = 0; i < _hg.k(); ++i ) {
      _variables[hyperedge_offset(he, i)] =
        _model.addVar(0.0, 1.0, 0.0, GRB_BINARY, hyperedge_var_desc(he, i));
      _variables[hyperedge_offset(he, i)].set(
        GRB_DoubleAttr_Start, (_hg.containsPinInPart(he, i)));
    }
  }
}

// ! Models the connectivity metric
void ILPModel::addObjectiveFunction() {
  GRBLinExpr connectivity_metric = 0;
  for ( const HypernodeID& he : _hg.edges() ) {
    GRBLinExpr connectivity = 0;
    for ( PartitionID i = 0; i < _hg.k(); ++i ) {
      connectivity += _variables[hyperedge_offset(he, i)];
    }
    connectivity_metric += (connectivity - 1) * _hg.edgeWeight(he);
  }
  _model.setObjective(connectivity_metric, GRB_MINIMIZE);
}

// ! Models constraints (1)
void ILPModel::restrictVerticesToOneBlock() {
  for ( const HypernodeID& hn : _hg.nodes() ) {
    GRBLinExpr only_one_block = 0;
    for ( PartitionID i = 0; i < _hg.k(); ++i ) {
      only_one_block += _variables[vertex_offset(hn, i)];
    }
    _model.addConstr(only_one_block == 1, "only_one_block_" + std::to_string(hn));
  }
}

// ! Models constraints (2)
void ILPModel::balanceConstraint() {
  for ( PartitionID i = 0; i < _hg.k(); ++i ) {
    const PartitionID original_block = _hg.toOriginalBlock(i);
    GRBLinExpr constraint = 0;
    for ( const HypernodeID& hn : _hg.nodes() ) {
      constraint += _variables[vertex_offset(hn, i)] * _hg.nodeWeight(hn);
    }
    _model.addConstr(constraint <= _context.partition.max_part_weights[original_block],
      "balance_constraint_" + std::to_string(i));
  }
}

// ! Models constraints (3)
void ILPModel::modelHyperedgeConnectivity() {
  for ( const HyperedgeID& he : _hg.edges() ) {
    for ( const HypernodeID& pin : _hg.pins(he) ) {
      for ( PartitionID i = 0; i < _hg.k(); ++i ) {
        _model.addConstr(_variables[hyperedge_offset(he,i)] >= _variables[vertex_offset(pin,i)],
          "connectivity_values_" + std::to_string(he) + "_" +
          std::to_string(pin) + "_" + std::to_string(i));
      }
    }
  }
}

} // namespace mt_kahypar

#pragma GCC diagnostic pop