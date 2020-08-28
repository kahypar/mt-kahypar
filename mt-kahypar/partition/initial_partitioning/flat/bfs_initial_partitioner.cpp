/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"

#include "mt-kahypar/partition/initial_partitioning/flat/policies/pseudo_peripheral_start_nodes.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

tbb::task* BFSInitialPartitioner::execute() {
    if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::bfs) ) {
      HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
      PartitionedHypergraph& hypergraph = _ip_data.local_partitioned_hypergraph();
      kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue =
        _ip_data.local_hypernode_fast_reset_flag_array();
      kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
        _ip_data.local_hyperedge_fast_reset_flag_array();

      _ip_data.reset_unassigned_hypernodes();
      parallel::scalable_vector<HypernodeID> start_nodes =
        PseudoPeripheralStartNodes::computeStartNodes(_ip_data, _context, kInvalidPartition, _rng);

      // Insert each start node for each block into its corresponding queue
      hypernodes_in_queue.reset();
      hyperedges_in_queue.reset();
      parallel::scalable_vector<Queue> queues(_context.partition.k);

      for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
        queues[block].push(start_nodes[block]);
        markHypernodeAsInQueue(hypergraph, hypernodes_in_queue, start_nodes[block], block);
      }

      HypernodeID num_assigned_hypernodes = 0;
      // We grow the k blocks of the partition starting from each start node in
      // a BFS-fashion. The BFS queues for each block are visited in round-robin-fashion.
      // Once a block is on turn, it pops it first hypernode and pushes
      // all adjacent vertices into its queue.
      const HypernodeID current_num_nodes =
        hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
      while ( num_assigned_hypernodes < current_num_nodes ) {
        for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
          HypernodeID hn = kInvalidHypernode;

          while ( !queues[block].empty() ) {
            const HypernodeID next_hn = queues[block].front();
            queues[block].pop();

            if ( hypergraph.partID(next_hn) == kInvalidPartition ) {
              // Hypernode is assigned to the current block, if it is not
              // assigned to an other block and if the assignment does not
              // violate the balanced constraint.
              // In case, there is no hypernode that fits into the current block,
              // we take the last unassigned hypernode popped from the queue.
              // Note, in that case the balanced constraint will be violated.
              hn = next_hn;
              if ( fitsIntoBlock(hypergraph, hn, block) ) {
                break;
              }
            }
          }

          if ( hn == kInvalidHypernode ) {
            // Special case, in case all hypernodes in the queue are already
            // assigned to an other block or the hypergraph is unconnected, we
            // choose an new unassigned hypernode (if one exists)
            hn = _ip_data.get_unassigned_hypernode();
          }

          if ( hn != kInvalidHypernode ) {
            ASSERT(hypergraph.partID(hn) == kInvalidPartition, V(block) << V(hypergraph.partID(hn)));
            hypergraph.setNodePart(hn, block);
            ++num_assigned_hypernodes;
            pushIncidentHypernodesIntoQueue(hypergraph, _context, queues[block],
              hypernodes_in_queue, hyperedges_in_queue, hn, block);
          } else {
            ASSERT(queues[block].empty());
          }
        }
      }

      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration<double>(end - start).count();
      _ip_data.commit(InitialPartitioningAlgorithm::bfs, time);
    }
    return nullptr;
}

// ! Pushes all adjacent hypernodes (not visited before) of hypernode hn
// ! into the BFS queue of the corresponding block.
inline void BFSInitialPartitioner::pushIncidentHypernodesIntoQueue(const PartitionedHypergraph& hypergraph,
                                                                   const Context& context,
                                                                   Queue& queue,
                                                                   kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                                                   kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                                                   const HypernodeID hn,
                                                                   const PartitionID block) {
  ASSERT(hn != kInvalidHypernode && block != kInvalidPartition);
  for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
    if ( !hyperedges_in_queue[block * hypergraph.initialNumEdges() + he] ) {
      if ( hypergraph.edgeSize(he) <= context.partition.ignore_hyperedge_size_threshold ) {
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( !hypernodes_in_queue[block * hypergraph.initialNumNodes() + pin] &&
                hypergraph.partID(pin) == kInvalidPartition ) {
            queue.push(pin);
            markHypernodeAsInQueue(hypergraph, hypernodes_in_queue, pin, block);
          }
        }
      }
      markHyperedgeAsInQueue(hypergraph, hyperedges_in_queue, he, block);
    }
  }
}

} // namespace mt_kahypar
