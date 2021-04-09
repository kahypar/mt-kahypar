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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"

#include "mt-kahypar/partition/refinement/ilp/ilp_model.h"

namespace mt_kahypar {

void ILPModel::initialize(const AdvancedProblem& problem) {
  _num_nodes = problem.stats.numNodes();
  _num_edges = problem.stats.numEdges();
  _k = problem.stats.numContainedBlocks();
  _initial_objective = 0;
  _variables.resize(_k * (_num_nodes + _num_edges));
  _contains_he_var.assign(_k * _num_edges, true);
  _nodes_to_ilp.clear();
  _edges_to_ilp.clear();
}

GRBModel ILPModel::construct(const PartitionedHypergraph& phg,
                             const AdvancedProblem& problem) {
  GRBModel model(_env);

  try {
    initialize(problem);
    addVariablesToModel(model, phg, problem);
    addObjectiveFunction(model, phg, problem);
    restrictVerticesToOneBlock(model, problem);
    balanceConstraint(model, phg, problem);
    modelHyperedgeConnectivity(model, phg, problem);
  } catch(GRBException e) {
    ERROR("Error code = " << e.getErrorCode() << "Message =" << e.getMessage());
  } catch(...) {
    ERROR("Exception during optimization");
  }

  return model;
}

PartitionID ILPModel::partID(const HypernodeID hn,
                             const AdvancedProblem& problem) {
  ASSERT(_nodes_to_ilp.contains(hn));
  const HypernodeID ilp_hn = _nodes_to_ilp[hn];
  for ( PartitionID i = 0; i < _k; ++i ) {
    if ( _variables[vertex_offset(ilp_hn, i)].get(GRB_DoubleAttr_X) > 0 ) {
      return problem.stats.indexToBlock(i);
    }
  }
  return kInvalidPartition;
}

void ILPModel::addVariablesToModel(GRBModel& model,
                                   const PartitionedHypergraph& phg,
                                   const AdvancedProblem& problem) {
  // Add hypernode variables
  for ( HypernodeID i = 0; i < problem.nodes.size(); ++i ) {
    const HypernodeID hn = problem.nodes[i];
    _nodes_to_ilp[hn] = i;
    ASSERT(problem.stats.isBlockContained(phg.partID(hn)));
    const PartitionID block = problem.stats.blockIndex(phg.partID(hn));
    for ( PartitionID p = 0; p < _k; ++p ) {
      const size_t offset = vertex_offset(i, p);
      ASSERT(offset < _variables.size());
      _variables[offset] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, vertex_var_desc(i, p));
      _variables[offset].set(GRB_DoubleAttr_Start, block == p);
    }
  }

  const vec<HyperedgeID>& edges = problem.stats.containedEdges();

  for ( HyperedgeID i = 0; i < edges.size(); ++i ) {
    const HyperedgeID he = edges[i];
    _edges_to_ilp[he] = i;
    // We only add hyperedge variables to our ILP problem, if all pins of
    // the corresponding block are completly contained in our ILP problem.
    // Otherwise, we can not remove the corresponding block from the hyperedge.
    for ( const HypernodeID& pin : phg.pins(he) ) {
      const PartitionID block = phg.partID(pin);
      if ( problem.stats.isBlockContained(block) &&
           !_nodes_to_ilp.contains(pin) ) {
        const PartitionID p = problem.stats.blockIndex(block);
        _contains_he_var[hyperedge_offset(i, p) - _k * _num_nodes] = false;
      }
    }

    PartitionID connectivity = 0;
    for ( PartitionID p = 0; p < problem.stats.numContainedBlocks(); ++p ) {
      const size_t offset = hyperedge_offset(i, p);
      ASSERT(offset < _variables.size());
      if ( _contains_he_var[offset - _k * _num_nodes] ) {
        _variables[offset] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, hyperedge_var_desc(i, p));
        const HypernodeID pin_count_in_part = phg.pinCountInPart(he, problem.stats.indexToBlock(p));
        _variables[offset].set(GRB_DoubleAttr_Start, pin_count_in_part > 0);
        connectivity += pin_count_in_part > 0;
      } else {
        ++connectivity;
      }
    }
    _initial_objective += (connectivity - 1) * phg.edgeWeight(he);
  }
}

GRBLinExpr ILPModel::getConnectivityMetric(const PartitionedHypergraph& phg,
                                           const AdvancedProblem& problem) {
  GRBLinExpr connectivity_metric = 0;
  const vec<HyperedgeID>& edges = problem.stats.containedEdges();
  for ( HyperedgeID he = 0; he < problem.stats.numEdges(); ++he ) {
    GRBLinExpr connectivity = 0;
    for ( PartitionID i = 0; i < _k; ++i ) {
      const size_t offset = hyperedge_offset(he, i);
      if ( _contains_he_var[offset - _k * _num_nodes] ) {
        connectivity += _variables[offset];
      } else {
        connectivity += 1;
      }
    }
    connectivity_metric += (connectivity - 1) * phg.edgeWeight(edges[he]);
  }
  return connectivity_metric;
}

void ILPModel::addObjectiveFunction(GRBModel& model,
                                    const PartitionedHypergraph& phg,
                                    const AdvancedProblem& problem) {
  GRBLinExpr objective_function = 0;
  if ( _context.partition.objective == kahypar::Objective::km1 ) {
    objective_function = getConnectivityMetric(phg, problem);
  } else {
    ERROR("Objective function not supported!");
  }
  model.setObjective(objective_function, GRB_MINIMIZE);
}

void ILPModel::restrictVerticesToOneBlock(GRBModel& model,
                                          const AdvancedProblem& problem) {
  for ( HypernodeID hn = 0; hn < problem.stats.numNodes(); ++hn ) {
    GRBLinExpr only_one_block = 0;
    for ( PartitionID i = 0; i < _k; ++i ) {
      only_one_block += _variables[vertex_offset(hn, i)];
    }
    model.addConstr(only_one_block == 1, "only_one_block_" + std::to_string(hn));
  }
}

void ILPModel::balanceConstraint(GRBModel& model,
                                const PartitionedHypergraph& phg,
                                const AdvancedProblem& problem) {
  for ( PartitionID i = 0; i < _k; ++i ) {
    const PartitionID original_block = problem.stats.indexToBlock(i);
    const HypernodeWeight l_max =
      problem.stats.maxPartWeight(original_block) -
      ( phg.partWeight(original_block) - problem.stats.nodeWeightOfBlock(original_block) );
    // Only add balance constraint for block i, if we would violate the balance constraint
    // if all vertices are assigned to that block
    if ( problem.stats.totalWeight() > l_max ) {
      GRBLinExpr constraint = 0;
      for ( HypernodeID hn = 0; hn < problem.stats.numNodes(); ++hn ) {
        constraint += _variables[vertex_offset(hn, i)] * phg.nodeWeight(problem.nodes[hn]);
      }
      model.addConstr(constraint <= l_max, "balance_constraint_" + std::to_string(i));
    }
  }
}

void ILPModel::modelHyperedgeConnectivity(GRBModel& model,
                                          const PartitionedHypergraph& phg,
                                          const AdvancedProblem& problem) {
  const vec<HyperedgeID>& edges = problem.stats.containedEdges();
  for ( const HyperedgeID& he : edges ) {
    const HyperedgeID ilp_he = _edges_to_ilp[he];
    for ( const HypernodeID& pin : phg.pins(he) ) {
      if ( _nodes_to_ilp.contains(pin) ) {
        const HypernodeID hn = _nodes_to_ilp[pin];
        for ( PartitionID i = 0; i < _k; ++i ) {
          const size_t he_offset = hyperedge_offset(ilp_he, i);
          if ( _contains_he_var[he_offset - _k * _num_nodes] ) {
            model.addConstr(_variables[he_offset] >= _variables[vertex_offset(hn,i)],
              "connectivity_values_" + std::to_string(he) + "_" +
              std::to_string(pin) + "_" + std::to_string(i));
          }
        }
      }
    }
  }
}

} // namespace mt_kahypar

#pragma GCC diagnostic pop