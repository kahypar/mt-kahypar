/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/flows/flow_refiner.h"

#include "tbb/concurrent_queue.h"

namespace mt_kahypar {

MoveSequence FlowRefiner::refineImpl(const PartitionedHypergraph& phg,
                                     const Subhypergraph& sub_hg,
                                     const HighResClockTimepoint& start) {
  MoveSequence sequence { { }, 0 };
  // Construct flow network that contains all vertices given in refinement nodes
  utils::Timer::instance().start_timer("construct_flow_network", "Construct Flow Network", true);
  FlowProblem flow_problem = constructFlowHypergraph(phg, sub_hg);
  utils::Timer::instance().stop_timer("construct_flow_network");
  if ( flow_problem.total_cut - flow_problem.non_removable_cut > 0 ) {

    // Solve max-flow min-cut problem
    bool time_limit_reached = false;
    utils::Timer::instance().start_timer("hyper_flow_cutter", "HyperFlowCutter", true);
    bool flowcutter_succeeded = runFlowCutter(flow_problem, start, time_limit_reached);
    utils::Timer::instance().stop_timer("hyper_flow_cutter");
    if ( flowcutter_succeeded ) {
      // We apply the solution if it either improves the cut or the balance of
      // the bipartition induced by the two blocks

      HyperedgeWeight new_cut = flow_problem.non_removable_cut;
      HypernodeWeight max_part_weight;
      const bool sequential = _context.shared_memory.num_threads == _context.refinement.flows.num_parallel_searches;
      if (sequential) {
        new_cut += _sequential_hfc.cs.flow_algo.flow_value;
        max_part_weight = std::max(_sequential_hfc.cs.source_weight, _sequential_hfc.cs.target_weight);
      } else {
        new_cut += _parallel_hfc.cs.flow_algo.flow_value;
        max_part_weight = std::max(_parallel_hfc.cs.source_weight, _parallel_hfc.cs.target_weight);
      }

      const bool improved_solution = new_cut < flow_problem.total_cut ||
        (new_cut == flow_problem.total_cut && max_part_weight < std::max(flow_problem.weight_of_block_0, flow_problem.weight_of_block_1));

      // Extract move sequence
      if ( improved_solution ) {
        sequence.expected_improvement = flow_problem.total_cut - new_cut;
        for ( const whfc::Node& u : _flow_hg.nodeIDs() ) {
          const HypernodeID hn = _whfc_to_node[u];
          if ( hn != kInvalidHypernode ) {
            const PartitionID from = phg.partID(hn);
            PartitionID to;
            if (sequential) {
              to = _sequential_hfc.cs.flow_algo.isSource(u) ? _block_0 : _block_1;
            } else {
              to = _parallel_hfc.cs.flow_algo.isSource(u) ? _block_0 : _block_1;
            }

            if ( from != to ) {
              sequence.moves.push_back(Move { from, to, hn, kInvalidGain });
            }
          }
        }
      }
    } else if ( time_limit_reached ) {
      sequence.state = MoveSequenceState::TIME_LIMIT;
    }
  }
  return sequence;
}

#define NOW std::chrono::high_resolution_clock::now()
#define RUNNING_TIME(X) std::chrono::duration<double>(NOW - X).count();

bool FlowRefiner::runFlowCutter(const FlowProblem& flow_problem,
                                const HighResClockTimepoint& start,
                                bool& time_limit_reached) {
  whfc::Node s = flow_problem.source;
  whfc::Node t = flow_problem.sink;
  bool result = false;

  size_t iteration = 0;
  auto on_cut = [&] {
    if (++iteration == 25) {
      iteration = 0;
      double elapsed = RUNNING_TIME(start);
      if (elapsed > _time_limit) {
        time_limit_reached = true;
        return false;
      }
    }
    return true;
  };


  const bool sequential = _context.shared_memory.num_threads == _context.refinement.flows.num_parallel_searches;
  if (sequential) {
    _sequential_hfc.cs.setMaxBlockWeight(0, std::max(
            flow_problem.weight_of_block_0, _context.partition.max_part_weights[_block_0]));
    _sequential_hfc.cs.setMaxBlockWeight(1, std::max(
            flow_problem.weight_of_block_1, _context.partition.max_part_weights[_block_1]));

    _sequential_hfc.reset();
    _sequential_hfc.setFlowBound(flow_problem.total_cut - flow_problem.non_removable_cut);
    result = _sequential_hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t, on_cut);
  } else {
    _parallel_hfc.cs.setMaxBlockWeight(0, std::max(
            flow_problem.weight_of_block_0, _context.partition.max_part_weights[_block_0]));
    _parallel_hfc.cs.setMaxBlockWeight(1, std::max(
            flow_problem.weight_of_block_1, _context.partition.max_part_weights[_block_1]));

    _parallel_hfc.reset();
    _parallel_hfc.setFlowBound(flow_problem.total_cut - flow_problem.non_removable_cut);
    result = _parallel_hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t, on_cut);
  }
  return result;
}

FlowProblem FlowRefiner::constructFlowHypergraph(const PartitionedHypergraph& phg,
                                                 const Subhypergraph& sub_hg) {
  _block_0 = sub_hg.block_0;
  _block_1 = sub_hg.block_1;
  ASSERT(_block_0 != kInvalidPartition && _block_1 != kInvalidPartition);
  FlowProblem flow_problem;


  const bool sequential = _context.shared_memory.num_threads == _context.refinement.flows.num_parallel_searches;
  if ( sequential ) {
    flow_problem = _sequential_construction.constructFlowHypergraph(
      phg, sub_hg, _block_0, _block_1, _whfc_to_node);
  } else {
    flow_problem = _parallel_construction.constructFlowHypergraph(
      phg, sub_hg, _block_0, _block_1, _whfc_to_node);
  }

  DBG << "Flow Hypergraph [ Nodes =" << _flow_hg.numNodes()
      << ", Edges =" << _flow_hg.numHyperedges()
      << ", Pins =" << _flow_hg.numPins()
      << ", Blocks = (" << _block_0 << "," << _block_1 << ") ]";

  return flow_problem;
}
} // namespace mt_kahypar
