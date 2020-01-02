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

#pragma once

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/pseudo_peripheral_start_nodes.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
template<typename TypeTraits>
class BFSInitialPartitionerT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;
  using Queue = parallel::scalable_queue<HypernodeID>;


  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

 public:
  BFSInitialPartitionerT(InitialPartitioningDataContainer& ip_data,
                         const Context& context) :
    _ip_data(ip_data),
    _context(context) { }

  tbb::task* execute() override {
    HyperGraph& hypergraph = _ip_data.local_hypergraph();
    kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue =
      _ip_data.local_hypernode_fast_reset_flag_array();
    kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
      _ip_data.local_hyperedge_fast_reset_flag_array();

    parallel::scalable_vector<HypernodeID> start_nodes =
      PseudoPeripheralStartNodes<TypeTraits>::computeStartNodes(_ip_data, _context);

    // Insert each start node for each block into its corresponding queue
    hypernodes_in_queue.reset();
    hyperedges_in_queue.reset();
    parallel::scalable_vector<Queue> queues(_context.partition.k);
    ASSERT(start_nodes.size() == static_cast<size_t>(_context.partition.k));
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      queues[block].push(start_nodes[block]);
      markHypernodeAsInQueue(hypergraph, hypernodes_in_queue, start_nodes[block], block);
    }

    _ip_data.reset_unassigned_hypernodes();
    HypernodeID num_assigned_hypernodes = 0;
    // We grow the k blocks of the partition starting from each start node in
    // a BFS-fashion. The BFS queues for each block are visited in round-robin-fashion.
    // Once a block is on turn, it pops it first hypernode and pushes
    // all adjacent vertices into its queue.
    while ( num_assigned_hypernodes < hypergraph.initialNumNodes() ) {
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

    _ip_data.commit(InitialPartitioningAlgorithm::bfs);
    return nullptr;
  }

 private:
  bool fitsIntoBlock(HyperGraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.localPartWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.max_part_weights[block];
  }

  // ! Pushes all adjacent hypernodes (not visited before) of hypernode hn
  // ! into the BFS queue of the corresponding block.
  inline void pushIncidentHypernodesIntoQueue(const HyperGraph& hypergraph,
                                              const Context& context,
                                              Queue& queue,
                                              kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                              kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                              const HypernodeID hn,
                                              const PartitionID block) {
    ASSERT(hn != kInvalidHypernode && block != kInvalidPartition);
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      const HyperedgeID original_he_id = hypergraph.originalEdgeID(he);
      if ( !hyperedges_in_queue[block * hypergraph.initialNumEdges() + original_he_id] ) {
        if ( hypergraph.edgeSize(he) <= context.partition.hyperedge_size_threshold ) {
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
            if ( !hypernodes_in_queue[block * hypergraph.initialNumNodes() + original_pin_id] &&
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

  inline void markHypernodeAsInQueue(const HyperGraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                     const HypernodeID hn,
                                     const PartitionID block) {
    hypernodes_in_queue.set(block * hypergraph.initialNumNodes() +
      hypergraph.originalNodeID(hn), true);
  }

  inline void markHyperedgeAsInQueue(const HyperGraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                     const HyperedgeID he,
                                     const PartitionID block) {
    hyperedges_in_queue.set(block * hypergraph.initialNumEdges() +
      hypergraph.originalEdgeID(he), true);
  }

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
};

template <typename TypeTraits>
PartitionID BFSInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID BFSInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using BFSInitialPartitioner = BFSInitialPartitionerT<GlobalTypeTraits>;

} // namespace mt_kahypar
