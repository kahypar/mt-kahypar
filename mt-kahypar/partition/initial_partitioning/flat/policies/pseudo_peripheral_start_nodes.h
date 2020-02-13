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

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template<typename TypeTraits>
class PseudoPeripheralStartNodes {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;
  using StartNodes = parallel::scalable_vector<HypernodeID>;
  using Queue = parallel::scalable_queue<HypernodeID>;

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

 public:
  static inline StartNodes computeStartNodes(InitialPartitioningDataContainer& ip_data,
                                             const Context& context) {
    HyperGraph& hypergraph = ip_data.local_partitioned_hypergraph();
    kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue =
      ip_data.local_hypernode_fast_reset_flag_array();
    kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
      ip_data.local_hyperedge_fast_reset_flag_array();
    int cpu_id = sched_getcpu();

    StartNodes start_nodes;
    HypernodeID start_hn = hypergraph.globalNodeID(
      utils::Randomize::instance().getRandomInt(0, hypergraph.initialNumNodes() - 1, cpu_id));
    ASSERT(hypergraph.nodeIsEnabled(start_hn));
    start_nodes.push_back(start_hn);

    // We perform k - 1 BFS on the hypergraph to find k vertices that
    // are "far" away from each other. Each BFS adds a new hypernode to
    // list of start nodes. Each entry in start_nodes represents a start
    // node for a specific block of the partition. The new vertex added to
    // the list of start nodes is the one last touched by the current BFS.
    for ( PartitionID i = 0; i < context.partition.k - 1; ++i ) {
      Queue queue;
      hypernodes_in_queue.reset();
      hyperedges_in_queue.reset();
      initializeQueue(hypergraph, queue, start_nodes, hypernodes_in_queue);

      HypernodeID last_hypernode_touched = kInvalidHypernode;
      HypernodeID num_touched_hypernodes = 0;
      ASSERT(queue.size() > 0);
      while ( !queue.empty() ) {
        last_hypernode_touched = queue.front();
        queue.pop();
        ++num_touched_hypernodes;

        // Add all adjacent non-visited vertices of the current visited hypernode
        // to queue.
        for ( const HyperedgeID& he : hypergraph.incidentEdges(last_hypernode_touched) ) {
          const HyperedgeID original_he_id = hypergraph.originalEdgeID(he);
          if ( !hyperedges_in_queue[original_he_id] ) {
            if ( hypergraph.edgeSize(he) <= context.partition.hyperedge_size_threshold ) {
              for ( const HypernodeID& pin : hypergraph.pins(he) ) {
                const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
                if ( !hypernodes_in_queue[original_pin_id] ) {
                  queue.push(pin);
                  hypernodes_in_queue.set(original_pin_id, true);
                }
              }
            }
            hyperedges_in_queue.set(original_he_id, true);
          }
        }

        // In case the queue is empty and we have not visited all hypernodes, we
        // add non-visited vertex to the queue (can happen if the hypergraph is not connected)
        if ( queue.empty() && num_touched_hypernodes < hypergraph.initialNumNodes() ) {
          for ( const HypernodeID& hn : hypergraph.nodes() ) {
            const HypernodeID original_hn_id = hypergraph.originalNodeID(hn);
            if ( !hypernodes_in_queue[original_hn_id] ) {
              queue.push(hn);
              hypernodes_in_queue.set(original_hn_id, true);
            }
          }
        }
      }

      // Add last touched hypernode of the BFS as new start node for block i + 1
      start_nodes.push_back(last_hypernode_touched);
    }

    ASSERT(start_nodes.size() == static_cast<size_t>(context.partition.k));
    return start_nodes;
  }

 private:
  static inline void initializeQueue(HyperGraph& hypergraph, Queue& queue, StartNodes& start_nodes,
                                     kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue) {
    for ( const HypernodeID& hn : start_nodes ) {
      queue.push(hn);
      hypernodes_in_queue.set(hypergraph.originalNodeID(hn), true);
    }
  }
};

template <typename TypeTraits>
PartitionID PseudoPeripheralStartNodes<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID PseudoPeripheralStartNodes<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

} // namespace mt_kahypar
