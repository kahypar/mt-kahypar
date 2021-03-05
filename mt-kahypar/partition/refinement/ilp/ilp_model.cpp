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

#include "kahypar/datastructure/fast_reset_flag_array.h"

namespace mt_kahypar {

/**
 * Constructs the following ILP Problem:
 *
 * Binary Decision Variables:
 * x_{v,i} -> If 1, then vertex v is assigned to block i
 * For connectivity metric:
 * y_{e,i} -> If 1, then hyperedge e contains at least one pin in block i
 * For cut metric:
 * y_e -> If 1, then hyperedge e is a cut hyperedge
 *
 * Constraints:
 *
 * Each block vertex v is assigned exactly to one block:
 * \forall v \in V: x_{v,1} + ... + x_{v,k} = 1  (1)
 *
 * Each block satisfies the balance constraint:
 * \forall k: \sum_{v \in V} x_{v,k} * c(v) <= L_max (2)
 *
 * Only for connectivity metric:
 * y_{e,i} is equal to 1, if at least one pin is assigned to block i
 * \forall e \in E: \forall v \in e: \forall k: y_{e,k} >= x_{e,k} (3)
 *
 * Only for cut metric
 * y_e is equal to 1, if two pins are assigned to different blocks.
 * We denote with v_{e,i} the i-th pin of hyperedge e
 * \forall e \in E: \forall i \in [2,|e|]: \forall k: y_e >= x_{v_{e,1},k} - x_{v_{e,i},k}
 *
 * Connecivity Metric:
 * minimize \sum_{e \in E} (\lambda(e) - 1) * \omega(e) with \lambda(e) = \sum_{i = 1}^k y_{e,i}
 *
 * Cut Metric:
 * minimize \sum_{e \in E} y_e * \omega(e)
 *
 * Optimizations:
 *  - remove super vertices that represents vertices of a block not contained in the ILP
 *    o For a supervertex v representing block i, we adapt the balance constraint to ... <= L_max - c(v)
 *    o If a supervertex v representing block i is contained in a hyperedge e, we can remove
 *      variable y_{e,i} because each solution will have y_{e,i} = 1 (adapt objective function accordingly).
 *  - If we can assign all vertices of the subhypergraph to a block without violating the balance constraint,
 *    we do not add the corresponding constraint to the ILP problem.
 *  - add balance constraint lazily, e.g. if optimal solution is found but did not satisfy the balance constraint (TODO)
 */
void ILPModel::construct(ILPHypergraph& hg) {
  try {
    _hg = &hg;
    addVariablesToModel(hg);
    addObjectiveFunction(hg);         // Connectivity Metric
    restrictVerticesToOneBlock(hg);   // Constraint (1)
    balanceConstraint(hg);            // Constraint (2)
    modelHyperedgeConstraint(hg);     // Constraint (3)
    _is_constructed = true;
  } catch(GRBException e) {
    ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
  } catch(...) {
    ERROR("Exception during optimization");
  }
}

void ILPModel::addVariablesToModel(ILPHypergraph& hg) {
  // Initialize variables
  _variables.resize(hg.numNodes() * hg.k() + hg.numEdges() * hg.k());
  _contains_variable.assign(hg.numNodes() * hg.k() + hg.numEdges() * hg.k(), true);
  _unremovable_block.assign(hg.numEdges(), kInvalidPartition);

  for ( const HypernodeID& hn : hg.nodes() ) {
    for ( PartitionID i = 0; i < hg.k(); ++i ) {
      _variables[vertex_offset(hn, i)] =
        _model.addVar(0.0, 1.0, 0.0, GRB_BINARY, vertex_var_desc(hn, i));
      _variables[vertex_offset(hn, i)].set(
        GRB_DoubleAttr_Start, (hg.partID(hn) == i));
    }
  }

  if ( _context.partition.objective == kahypar::Objective::km1 ) {
    addHyperedgeVariablesToModelForConnectivityMetric(hg);
  } else if ( _context.partition.objective == kahypar::Objective::cut ) {
    addHyperedgeVariablesToModelForCutMetric(hg);
  }
}

// ! Counts the number of super vertices contained in hyperedge he
// ! Super vertices are not part of the ILP and consequently cannot be moved
PartitionID ILPModel::determineNumberOfUnremovableBlocks(ILPHypergraph& hg, const HyperedgeID he) {
  PartitionID num_unremovable_blocks = 0;
  for ( const HypernodeID& pin : hg.pins(he) ) {
    // Super vertices can not change their block
    if ( hg.isSuperVertex(pin) ) {
      const PartitionID block = hg.superVertexBlock(pin);
      size_t he_offset = hyperedge_offset(he, block);
      if ( _contains_variable[he_offset] ) {
        ++num_unremovable_blocks;
        _contains_variable[he_offset] = false;
        _unremovable_block[he] = block;
      }
    }
  }
  return num_unremovable_blocks;
}

void ILPModel::addHyperedgeVariablesToModelForConnectivityMetric(ILPHypergraph& hg) {
  for ( const HypernodeID& he : hg.edges() ) {
    PartitionID connectivity = determineNumberOfUnremovableBlocks(hg, he);
    for ( PartitionID i = 0; i < hg.k(); ++i ) {
      // Only add variable y_{e,i} to ILP, if we can potentially
      // remove block i from hyperedge e
      const size_t he_offset = hyperedge_offset(he, i);
      if ( _contains_variable[he_offset] ) {
        _variables[he_offset] =
          _model.addVar(0.0, 1.0, 0.0, GRB_BINARY, hyperedge_var_desc(he, i));
        const bool contains_pin_in_block = hg.containsPinInPart(he, i);
        _variables[he_offset].set(
          GRB_DoubleAttr_Start, contains_pin_in_block);
        connectivity += contains_pin_in_block;
      }
    }
    _initial_objective += (connectivity - 1) * hg.edgeWeight(he);
  }
}

void ILPModel::addHyperedgeVariablesToModelForCutMetric(ILPHypergraph& hg) {
  for ( const HypernodeID& he : hg.edges() ) {
    // Only add variable y_e to ILP, if we can potentially
    // remove hyperedge e from cut
    const bool is_cut = hg.isCut(he);
    const PartitionID num_unremovable_blocks = determineNumberOfUnremovableBlocks(hg, he);
    if ( num_unremovable_blocks < 2 ) {
      _variables[hyperedge_offset(he, 0)] =
        _model.addVar(0.0, 1.0, 0.0, GRB_BINARY, hyperedge_var_desc(he, 0));
      _variables[hyperedge_offset(he, 0)].set(GRB_DoubleAttr_Start, is_cut);
      _contains_variable[hyperedge_offset(he, 0)] = true;
    } else {
      _contains_variable[hyperedge_offset(he, 0)] = false;
    }
    _initial_objective += is_cut * hg.edgeWeight(he);
  }
}

// ! Models the connectivity metric
void ILPModel::addObjectiveFunction(ILPHypergraph& hg) {
  if ( _context.partition.objective == kahypar::Objective::km1 ) {
    addConnectivityMetric(hg);
  } else if ( _context.partition.objective == kahypar::Objective::cut ) {
    addCutMetric(hg);
  }
}

void ILPModel::addConnectivityMetric(ILPHypergraph& hg) {
  GRBLinExpr connectivity_metric = 0;
  for ( const HypernodeID& he : hg.edges() ) {
    GRBLinExpr connectivity = 0;
    for ( PartitionID i = 0; i < hg.k(); ++i ) {
      const size_t he_offset = hyperedge_offset(he, i);
      if ( _contains_variable[he_offset] ) {
        connectivity += _variables[he_offset];
      } else {
        connectivity += 1;
      }
    }
    connectivity_metric += (connectivity - 1) * hg.edgeWeight(he);
  }
  _model.setObjective(connectivity_metric, GRB_MINIMIZE);
}

void ILPModel::addCutMetric(ILPHypergraph& hg) {
  GRBLinExpr cut_metric = 0;
  for ( const HypernodeID& he : hg.edges() ) {
    const size_t he_offset = hyperedge_offset(he, 0);
    cut_metric += _contains_variable[he_offset] ?
      _variables[he_offset] * hg.edgeWeight(he) : hg.edgeWeight(he);
  }
  _model.setObjective(cut_metric, GRB_MINIMIZE);
}

// ! Models constraint (1)
void ILPModel::restrictVerticesToOneBlock(ILPHypergraph& hg) {
  for ( const HypernodeID& hn : hg.nodes() ) {
    GRBLinExpr only_one_block = 0;
    for ( PartitionID i = 0; i < hg.k(); ++i ) {
      only_one_block += _variables[vertex_offset(hn, i)];
    }
    _model.addConstr(only_one_block == 1, "only_one_block_" + std::to_string(hn));
  }
}

// ! Models constraint (2)
void ILPModel::balanceConstraint(ILPHypergraph& hg) {
  for ( PartitionID i = 0; i < hg.k(); ++i ) {
    const PartitionID original_block = hg.toOriginalBlock(i);
    const HypernodeWeight l_max =
      _context.partition.max_part_weights[original_block] - hg.superVertexWeight(i);
    // Only add balance constraint for block i, if we would violate the balance constraint
    // if all vertices are assigned to that block
    if ( hg.subhypergraphWeight() > l_max ) {
      GRBLinExpr constraint = 0;
      for ( const HypernodeID& hn : hg.nodes() ) {
        constraint += _variables[vertex_offset(hn, i)] * hg.nodeWeight(hn);
      }
      _model.addConstr(constraint <= l_max, "balance_constraint_" + std::to_string(i));
    }
  }
}

// ! Models constraint (3)
void ILPModel::modelHyperedgeConstraint(ILPHypergraph& hg) {
  if ( _context.partition.objective == kahypar::Objective::km1 ) {
    modelHyperedgeConnectivity(hg);
  } else if ( _context.partition.objective == kahypar::Objective::cut ) {
    modelHyperedgeCut(hg);
  }
}

// ! Models constraint for connectivity metric (3)
void ILPModel::modelHyperedgeConnectivity(ILPHypergraph& hg) {
  for ( const HyperedgeID& he : hg.edges() ) {
    for ( const HypernodeID& pin : hg.pins(he) ) {
      if ( !hg.isSuperVertex(pin) ) {
        for ( PartitionID i = 0; i < hg.k(); ++i ) {
          const size_t he_offset = hyperedge_offset(he,i);
          if ( _contains_variable[he_offset] ) {
            _model.addConstr(_variables[he_offset] >= _variables[vertex_offset(pin,i)],
              "connectivity_values_" + std::to_string(he) + "_" +
              std::to_string(pin) + "_" + std::to_string(i));
          }
        }
      }
    }
  }
}

// ! Models constraint for cut metric (3)
void ILPModel::modelHyperedgeCut(ILPHypergraph& hg) {
  for ( const HyperedgeID& he : hg.edges() ) {
    const size_t he_offset = hyperedge_offset(he,0);
    if ( _contains_variable[he_offset] ) {
      const bool contains_unremovable_block = _unremovable_block[he] != kInvalidPartition;
      if ( contains_unremovable_block ) {
        // If hyperedge he contains an unremovable block, we must assign all pins
        // to corresponding block to remove he from the cut
        const PartitionID i = _unremovable_block[he];
        for ( const HypernodeID& pin : hg.pins(he) ) {
          if ( !hg.isSuperVertex(pin) ) {
            _model.addConstr(_variables[he_offset] >= 1 - _variables[vertex_offset(pin,i)],
              "cut_values_" + std::to_string(he) + "_" +
              std::to_string(pin) + "_" + std::to_string(i));
          }
        }
      } else {
        // Otherwise, the hyperedge is cut if two pins of the hyperedge are
        // assigned to different blocks
        HypernodeID first_vertex = kInvalidHypernode;
        for ( const HypernodeID& pin : hg.pins(he) ) {
          if ( !hg.isSuperVertex(pin) && first_vertex != kInvalidHypernode && first_vertex != pin ) {
            for ( PartitionID i = 0; i < hg.k(); ++i ) {
              _model.addConstr(_variables[he_offset] >=
                _variables[vertex_offset(first_vertex,i)] - _variables[vertex_offset(pin,i)],
                "cut_values_" + std::to_string(he) + "_" +
                std::to_string(pin) + "_" + std::to_string(i));
            }
          } else if ( !hg.isSuperVertex(pin) && first_vertex == kInvalidHypernode ) {
            first_vertex = pin;
          }
        }
      }
    }
  }
}

} // namespace mt_kahypar

#pragma GCC diagnostic pop