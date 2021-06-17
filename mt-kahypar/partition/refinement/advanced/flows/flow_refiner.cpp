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
  MoveSequence sequence { { }, 0 };
  // Construct flow network that contains all vertices given in refinement nodes
  FlowProblem flow_problem = constructFlowHypergraph(phg, refinement_nodes);
  if ( flow_problem.total_cut - flow_problem.non_removable_cut > 0 ) {
    // Set maximum allowed block weights for block 0 and 1
    _hfc.cs.setMaxBlockWeight(0, std::max(
      flow_problem.weight_of_block_0, _context.partition.max_part_weights[_block_0]));
    _hfc.cs.setMaxBlockWeight(1, std::max(
      flow_problem.weight_of_block_1, _context.partition.max_part_weights[_block_1]));

    _hfc.reset();
    _hfc.upperFlowBound = flow_problem.total_cut - flow_problem.non_removable_cut;
    // Solve max-flow min-cut problem
    bool flowcutter_succeeded =
      _hfc.runUntilBalancedOrFlowBoundExceeded(flow_problem.source, flow_problem.sink);
    if ( flowcutter_succeeded ) {
      // We apply the solution if it either improves the cut or the balance of
      // the bipartition induced by the two blocks
      HyperedgeWeight new_cut = flow_problem.non_removable_cut + _hfc.cs.flowValue;
      const bool improved_solution = new_cut < flow_problem.total_cut ||
        ( new_cut == flow_problem.total_cut &&
          static_cast<HypernodeWeight>(std::max(_hfc.cs.n.sourceWeight, _hfc.cs.n.targetWeight)) <
          std::max(flow_problem.weight_of_block_0, flow_problem.weight_of_block_1));

      // Extract move sequence
      if ( improved_solution ) {
        sequence.expected_improvement = flow_problem.total_cut - new_cut;
        for ( const HypernodeID& hn : refinement_nodes ) {
          const PartitionID from = phg.partID(hn);
          const PartitionID to = _hfc.cs.n.isSource(_node_to_whfc[hn]) ? _block_0 : _block_1;
          if ( from != to ) {
            sequence.moves.push_back(Move { from, to, hn, kInvalidGain });
          }
        }
      }
    }
  }
  return sequence;
}

FlowRefiner::FlowProblem FlowRefiner::constructFlowHypergraph(const PartitionedHypergraph& phg,
                                                              const vec<HypernodeID>& refinement_nodes) {
  ASSERT(_block_0 != kInvalidPartition && _block_1 != kInvalidPartition);
  FlowProblem flow_problem;
  flow_problem.total_cut = 0;
  flow_problem.non_removable_cut = 0;

  if ( _context.refinement.advanced.flows.determine_distance_from_cut ) {
    _cut_hes.clear();
  }

  // Add refinement nodes to flow network
  whfc::Node flow_hn(0);
  HypernodeWeight weight_block_0 = 0;
  HypernodeWeight weight_block_1 = 0;
  auto add_nodes = [&](const PartitionID block, HypernodeWeight& weight_of_block) {
    for ( const HypernodeID& u : refinement_nodes) {
      ASSERT(phg.partID(u) == _block_0 || phg.partID(u) == _block_1);
      if ( phg.partID(u) == block ) {
        const HypernodeWeight u_weight = phg.nodeWeight(u);
        _node_to_whfc[u] = flow_hn++;
        _flow_hg.addNode(whfc::NodeWeight(u_weight));
        weight_of_block += u_weight;
        for ( const HyperedgeID& he : phg.incidentEdges(u) ) {
          if ( _context.refinement.advanced.flows.determine_distance_from_cut &&
               phg.pinCountInPart(he, _block_0) > 0 && phg.pinCountInPart(he, _block_1) > 0 &&
               !_visited_hes.contains(he) ) {
            _cut_hes.push_back(he);
          }
          _visited_hes[he] = ds::EmptyStruct { };
        }
      }
    }
  };
  // Add source nodes
  flow_problem.source = flow_hn++;
  _flow_hg.addNode(whfc::NodeWeight(0));
  add_nodes(_block_0, weight_block_0);
  _flow_hg.nodeWeight(flow_problem.source) = whfc::NodeWeight(std::max(0, phg.partWeight(_block_0) - weight_block_0));
  // Add sink nodes
  flow_problem.sink = flow_hn++;
  _flow_hg.addNode(whfc::NodeWeight(0));
  add_nodes(_block_1, weight_block_1);
  _flow_hg.nodeWeight(flow_problem.sink) = whfc::NodeWeight(std::max(0, phg.partWeight(_block_1) - weight_block_1));
  flow_problem.weight_of_block_0 = _flow_hg.nodeWeight(flow_problem.source) + weight_block_0;
  flow_problem.weight_of_block_1 = _flow_hg.nodeWeight(flow_problem.sink) + weight_block_1;

  // Add hyperedge to flow network and configure source and sink
  for ( const auto& entry : _visited_hes ) {
    const HyperedgeID he = entry.key;
    if ( !canHyperedgeBeDropped(phg, he) ) {
      _tmp_pins.clear();
      const HyperedgeWeight he_weight = phg.edgeWeight(he);
      _flow_hg.startHyperedge(whfc::Flow(he_weight));
      bool connectToSource = false;
      bool connectToSink = false;
      if ( phg.pinCountInPart(he, _block_0) > 0 && phg.pinCountInPart(he, _block_1) > 0 ) {
        flow_problem.total_cut += he_weight;
      }
      for ( const HypernodeID& pin : phg.pins(he) ) {
        if ( _node_to_whfc.contains(pin) ) {
          _tmp_pins.push_back(pin);
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
      } else {
        // According to Lars: Adding to source or sink to the start of
        // each pin list improves running time
        if ( connectToSource ) {
          _flow_hg.addPin(flow_problem.source);
        } else if ( connectToSink ) {
          _flow_hg.addPin(flow_problem.sink);
        }
        for ( const HypernodeID& pin : _tmp_pins ) {
          _flow_hg.addPin(_node_to_whfc[pin]);
        }
      }
    }
  }

  if ( _context.refinement.advanced.flows.determine_distance_from_cut ) {
    // Determine the distance of each node contained in the flow network from the cut.
    // This technique improves piercing decision within the WHFC framework.
    determineDistanceFromCut(phg, flow_problem.source, flow_problem.sink);
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

void FlowRefiner::determineDistanceFromCut(const PartitionedHypergraph& phg,
                                           const whfc::Node source,
                                           const whfc::Node sink) {
  _hfc.cs.borderNodes.distance.distance.assign(_flow_hg.numNodes(), whfc::HopDistance(0));
  _visited_hes.clear();
  _visited_hns.clear();

  // Initialize bfs queue with vertices contained in cut hyperedges
  parallel::scalable_queue<HypernodeID> q;
  parallel::scalable_queue<HypernodeID> next_q;
  for ( const HyperedgeID& he : _cut_hes ) {
    for ( const HypernodeID& pin : phg.pins(he) ) {
      if ( _node_to_whfc.contains(pin) && !_visited_hns.contains(pin) ) {
        q.push(pin);
        _visited_hns[pin] = ds::EmptyStruct { };
      }
    }
    _visited_hes[he] = ds::EmptyStruct { };
  }

  // Perform BFS to determine distance of each vertex from cut
  whfc::HopDistance dist(1);
  whfc::HopDistance max_dist_source(0);
  whfc::HopDistance max_dist_sink(0);
  while ( !q.empty() ) {
    const HypernodeID u = q.front();
    q.pop();

    if ( phg.partID(u) == _block_0 ) {
      _hfc.cs.borderNodes.distance[_node_to_whfc[u]] = -dist;
      max_dist_source = std::max(max_dist_source, dist);
    } else {
      ASSERT(phg.partID(u) == _block_1);
      _hfc.cs.borderNodes.distance[_node_to_whfc[u]] = dist;
      max_dist_sink = std::max(max_dist_sink, dist);
    }

    for ( const HyperedgeID& he : phg.incidentEdges(u) ) {
      if ( !_visited_hes.contains(he) ) {
        for ( const HypernodeID& pin : phg.pins(he) ) {
          if ( _node_to_whfc.contains(pin) && !_visited_hns.contains(pin) ) {
            next_q.push(pin);
            _visited_hns[pin] = ds::EmptyStruct { };
          }
        }
        _visited_hes[he] = ds::EmptyStruct { };
      }
    }

    if ( q.empty() ) {
      std::swap(q, next_q);
      ++dist;
    }
  }
  _hfc.cs.borderNodes.distance[source] = -(max_dist_source + 1);
  _hfc.cs.borderNodes.distance[sink] = max_dist_sink + 1;
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