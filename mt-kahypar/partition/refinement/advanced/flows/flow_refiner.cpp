/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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


#include "mt-kahypar/partition/refinement/advanced/flows/flow_refiner.h"

namespace mt_kahypar {

MoveSequence FlowRefiner::refineImpl(const PartitionedHypergraph& phg,
                                     const vec<HypernodeID>& refinement_nodes) {
  constructFlowHypergraph(phg, refinement_nodes);
  return MoveSequence { { }, 0 };
}

FlowRefiner::FlowProblem FlowRefiner::constructFlowHypergraph(const PartitionedHypergraph& phg,
                                                              const vec<HypernodeID>& refinement_nodes) {
  ASSERT(_block_0 != kInvalidPartition && _block_1 != kInvalidPartition);
  FlowProblem flow_problem;
  flow_problem.total_cut = 0;
  flow_problem.non_removable_cut = 0;
  // Add refinement nodes to flow network
  whfc::Node flow_hn(0);
  HypernodeWeight weight_block_0 = 0;
  HypernodeWeight weight_block_1 = 0;
  for ( const HypernodeID& u : refinement_nodes) {
    const HypernodeWeight u_weight = phg.nodeWeight(u);
    _node_to_whfc[u] = flow_hn++;
    _flow_hg.addNode(whfc::NodeWeight(u_weight));
    if ( phg.partID(u) == _block_0 ) {
      weight_block_0 += u_weight;
    } else {
      ASSERT(phg.partID(u) == _block_1);
      weight_block_1 += u_weight;
    }
    for ( const HyperedgeID& he : phg.incidentEdges(u) ) {
      _visited_hes[he] = ds::EmptyStruct { };
    }
  }
  // TODO(heuer): As mentioned by Lars, this could be a potentially source for bad quality.
  // In KaHyPar, the source has label 0 and the sink has label greater than all nodes in block 0.
  // Some weird implementation details might cause here some quality issues.
  flow_problem.source = flow_hn++;
  flow_problem.sink = flow_hn++;
  _flow_hg.addNode(whfc::NodeWeight(std::max(0, phg.partWeight(_block_0) - weight_block_0)));
  _flow_hg.addNode(whfc::NodeWeight(std::max(0, phg.partWeight(_block_1) - weight_block_1)));

  // Add hyperedge to flow network and configure source and sink
  for ( const auto& entry : _visited_hes ) {
    const HyperedgeID he = entry.key;
    if ( !canHyperedgeBeDropped(phg, he) ) {
      const HyperedgeWeight he_weight = phg.edgeWeight(he);
      _flow_hg.startHyperedge(whfc::Flow(he_weight));
      bool connectToSource = false;
      bool connectToSink = false;
      if ( phg.pinCountInPart(he, _block_0) > 0 && phg.pinCountInPart(he, _block_1) > 0 ) {
        flow_problem.total_cut += he_weight;
      }
      for ( const HypernodeID& pin : phg.pins(he) ) {
        if ( _node_to_whfc.contains(pin) ) {
          _flow_hg.addPin(_node_to_whfc[pin]);
        } else {
          connectToSource |= phg.partID(pin) == _block_0;
          connectToSink |= phg.partID(pin) == _block_1;
        }
      }

      if ( connectToSource && connectToSink ) {
        // Hyperedge is connected to source and sink which means we can not remove it
        // from the cut with the current flow problem => remove he from flow problem
        _flow_hg.removeCurrentHyperedge();
        flow_problem.non_removable_cut += he_weight;
      } else if ( connectToSource ) {
        _flow_hg.addPin(flow_problem.source);
      } else if ( connectToSink ) {
        _flow_hg.addPin(flow_problem.sink);
      }
    }
  }

  if ( _flow_hg.nodeWeight(flow_problem.source) == 0 ||
       _flow_hg.nodeWeight(flow_problem.sink) == 0 ) {
    // Source or sink not connected to vertices in the flow problem
    flow_problem.non_removable_cut = 0;
    flow_problem.total_cut = 0;
  } else {
  _flow_hg.finalize();
  }
  return flow_problem;
}

bool FlowRefiner::isMaximumProblemSizeReachedImpl(ProblemStats& stats) const {
  ASSERT(_phg);
  ASSERT(stats.numContainedBlocks() == 2);
  _block_0 = stats.block(0);
  _block_1 = stats.block(1);
  const HypernodeWeight max_weight_0 =
    _scaling * _context.partition.perfect_balance_part_weights[_block_1] - _phg->partWeight(_block_1);
  const HypernodeWeight max_weight_1 =
    _scaling * _context.partition.perfect_balance_part_weights[_block_0] - _phg->partWeight(_block_0);
  if ( stats.nodeWeightOfBlock(_block_0) >= max_weight_0 ) {
    stats.lockBlock(_block_0);
  }
  if ( stats.nodeWeightOfBlock(_block_1) >= max_weight_1 ) {
    stats.lockBlock(_block_1);
  }
  return stats.isLocked(_block_0) && stats.isLocked(_block_1);
}

} // namespace mt_kahypar