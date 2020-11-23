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

class PseudoPeripheralStartNodes {
  using StartNodes = parallel::scalable_vector<HypernodeID>;
  using Queue = parallel::scalable_queue<HypernodeID>;

  static constexpr bool debug = false;

 public:
  static inline StartNodes computeStartNodes(InitialPartitioningDataContainer& ip_data,
                                             const Context& context,
                                             const PartitionID default_block,
                                             std::mt19937& rng) {
    PartitionedHypergraph& hypergraph = ip_data.local_partitioned_hypergraph();
    kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue =
      ip_data.local_hypernode_fast_reset_flag_array();
    kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
      ip_data.local_hyperedge_fast_reset_flag_array();

    ASSERT(hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes() >= ID(hypergraph.k()));
    StartNodes start_nodes;
    HypernodeID start_hn =
            std::uniform_int_distribution<HypernodeID>(0, hypergraph.initialNumNodes() -1 )(rng);
    if ( !hypergraph.nodeIsEnabled(start_hn) ) {
      start_hn = ip_data.get_unassigned_hypernode(default_block);
    }
    ASSERT(start_hn != kInvalidHypernode && hypergraph.nodeIsEnabled(start_hn));
    start_nodes.push_back(start_hn);

    // We perform k - 1 BFS on the hypergraph to find k vertices that
    // are "far" away from each other. Each BFS adds a new hypernode to
    // list of start nodes. Each entry in start_nodes represents a start
    // node for a specific block of the partition. The new vertex added to
    // the list of start nodes is the one last touched by the current BFS.
    const HypernodeID current_num_nodes =
      hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
    parallel::scalable_vector<HypernodeID> non_touched_hypernodes;
    for ( PartitionID i = 0; i < context.partition.k - 1; ++i ) {
      Queue queue;
      hypernodes_in_queue.reset();
      hyperedges_in_queue.reset();
      initializeQueue(queue, start_nodes, hypernodes_in_queue);

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
          if ( !hyperedges_in_queue[he] ) {
            if ( hypergraph.edgeSize(he) <= context.partition.ignore_hyperedge_size_threshold ) {
              for ( const HypernodeID& pin : hypergraph.pins(he) ) {
                if ( !hypernodes_in_queue[pin] ) {
                  queue.push(pin);
                  hypernodes_in_queue.set(pin, true);
                }
              }
            }
            hyperedges_in_queue.set(he, true);
          }
        }

        // In case the queue is empty and we have not visited all hypernodes.
        // Therefore, we choose one unvisited vertex at random.
        if ( queue.empty() && num_touched_hypernodes < current_num_nodes ) {
          for ( const HypernodeID& hn : hypergraph.nodes() ) {
            if ( !hypernodes_in_queue[hn] ) {
              non_touched_hypernodes.push_back(hn);
              hypernodes_in_queue.set(hn, true);
            }
          }
          const int rand_idx = utils::Randomize::instance().getRandomInt(
            0, non_touched_hypernodes.size() - 1, sched_getcpu());
          last_hypernode_touched = non_touched_hypernodes[rand_idx];
        }
      }

      // Add last touched hypernode of the BFS as new start node for block i + 1
      start_nodes.push_back(last_hypernode_touched);
      non_touched_hypernodes.clear();
    }

    ASSERT(start_nodes.size() == static_cast<size_t>(context.partition.k));
    return start_nodes;
  }

 private:
  static inline void initializeQueue(Queue& queue, StartNodes& start_nodes,
                                     kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue) {
    for ( const HypernodeID& hn : start_nodes ) {
      queue.push(hn);
      hypernodes_in_queue.set(hn, true);
    }
  }
};


} // namespace mt_kahypar
