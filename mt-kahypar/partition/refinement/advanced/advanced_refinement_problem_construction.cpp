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


#include "mt-kahypar/partition/refinement/advanced/advanced_refinement_problem_construction.h"

#include "tbb/parallel_for.h"

namespace mt_kahypar {

void AdvancedRefinementProblemConstruction::BFSData::reset() {
  current_distance = 0;
  last_queue_idx = 0;
  while ( !queue[0].empty() ) queue[0].pop();
  while ( !queue[1].empty() ) queue[1].pop();
  while ( !next_queue[0].empty() ) next_queue[0].pop();
  while ( !next_queue[1].empty() ) next_queue[1].pop();
  std::fill(visited_hn.begin(), visited_hn.end(), false);
  std::fill(visited_he.begin(), visited_he.end(), false);
}

void AdvancedRefinementProblemConstruction::BFSData::pop_hypernode(HypernodeID& hn) {
  if ( !is_empty() ) {
    // Pop vertices alternating from one of the two queues
    size_t idx = last_queue_idx++ % 2;
    if ( queue[idx].empty() ) {
      idx = last_queue_idx++ % 2;
    }
    ASSERT(!queue[idx].empty());
    hn = queue[idx].front();
    queue[idx].pop();
  }
}

void AdvancedRefinementProblemConstruction::BFSData::add_pins_of_hyperedge_to_queue(
  const HyperedgeID& he,
  const PartitionedHypergraph& phg,
  const AdvancedProblemStats& stats,
  const size_t max_bfs_distance) {
  if ( current_distance <= max_bfs_distance ) {
    if ( !visited_he[he] ) {
      for ( const HypernodeID& pin : phg.pins(he) ) {
        const PartitionID block = phg.partID(pin);
        if ( !stats.isLocked(block) &&
             (blocks.i == block || blocks.j == block) &&
             !visited_hn[pin] ) {
          next_queue[blocks.i == block ? 0 : 1].push(pin);
          visited_hn[pin] = true;
        }
      }
      visited_he[he] = true;
    }
  }
}

void AdvancedRefinementProblemConstruction::ConstructionData::initialize(
  const vec<BlockPairCutHyperedges>& initial_cut_hes,
  const AdvancedProblemStats& stats,
  const PartitionedHypergraph& phg) {
  // Initialize BFS Queues
  used_slots = initial_cut_hes.size();
  while ( bfs.size() <= used_slots ) bfs.emplace_back(_num_nodes, _num_edges);
  for ( size_t i = 0; i < used_slots; ++i ) {
    bfs[i].reset();
    bfs[i].blocks = initial_cut_hes[i].blocks;
    for ( const HyperedgeID& he : initial_cut_hes[i].cut_hes ) {
      bfs[i].add_pins_of_hyperedge_to_queue(
        he, phg, stats, std::numeric_limits<size_t>::max());
    }
    bfs[i].swap_with_next_queue();
  }
}

void AdvancedRefinementProblemConstruction::ConstructionData::pop_hypernode(HypernodeID& hn, size_t& idx) {
  if ( !is_empty() ) {
    // BFS Queues are visited in round-robin fashion
    bool found = false;
    while ( !found ) {
      idx = last_idx++ % used_slots;
      // If the current bfs queue is empty,
      // we swap it with queue for the next layer
      if ( bfs[idx].is_empty() ) {
        bfs[idx].swap_with_next_queue();
      }
      if ( !bfs[idx].is_empty() ) {
        bfs[idx].pop_hypernode(hn);
        ASSERT(hn != kInvalidHypernode);
        found = true;
      }
    }
  }
}

void AdvancedRefinementProblemConstruction::ConstructionData::clearBlock(const PartitionID block) {
  auto clear_queue = [](parallel::scalable_queue<HypernodeID>& queue) {
    while ( !queue.empty() ) queue.pop();
  };

  for ( size_t i = 0; i < used_slots; ++i ) {
    if ( bfs[i].blocks.i == block ) {
      clear_queue(bfs[i].queue[0]);
      clear_queue(bfs[i].next_queue[0]);
    } else if ( bfs[i].blocks.j == block ) {
      clear_queue(bfs[i].queue[1]);
      clear_queue(bfs[i].next_queue[1]);
    }
  }
}

vec<HypernodeID> AdvancedRefinementProblemConstruction::construct(const SearchID search_id,
                                                                  QuotientGraph& quotient_graph,
                                                                  AdvancedRefinerAdapter& refiner,
                                                                  const PartitionedHypergraph& phg) {
  vec<HypernodeID> nodes;

  ConstructionData& data = _local_data.local();
  AdvancedProblemStats& stats = _local_stats.local();
  stats.reset();
  for ( const BlockPair& blocks : quotient_graph.getBlockPairs(search_id) ) {
    stats.addBlock(blocks.i);
    stats.addBlock(blocks.j);
  }

  const size_t num_block_pairs = quotient_graph.numBlockPairs(search_id);
  // We vertices to the problem as long as the associated refiner notifies the
  // construction algorithm that the maximum problem size is reached
  while ( !refiner.isMaximumProblemSizeReached(search_id, stats) ) {

    // We initialize the BFS with a fixed number of cut hyperedges running
    // between the involved block associated with the search
    const vec<BlockPairCutHyperedges> initial_cut_hes =
      quotient_graph.requestCutHyperedges(search_id, num_block_pairs *
        _context.refinement.advanced.num_cut_edges_per_block_pair);
    data.initialize(initial_cut_hes, stats, phg);
    // Special case, if they are no cut hyperedges left
    // between the involved blocks
    if ( data.is_empty() ) break;

    // BFS
    while ( !data.is_empty() &&
            !refiner.isMaximumProblemSizeReached(search_id, stats) ) {
      size_t queue_idx = std::numeric_limits<size_t>::max();
      HypernodeID hn = kInvalidHypernode;
      data.pop_hypernode(hn, queue_idx);
      ASSERT(hn != kInvalidHypernode);
      ASSERT(queue_idx != std::numeric_limits<size_t>::max());

      PartitionID block = phg.partID(hn);
      if ( !stats.isLocked(block) ) {
        // Search aquires ownership of the vertex. Each vertex is only allowed to
        // be part of one search at any time.
        if ( acquire_vertex(search_id, hn) ) {
          block = phg.partID(hn);
          // Double-check if vertex is still part of the blocks associated
          // with the search.
          if ( stats.isBlockContained(block) ) {
            nodes.push_back(hn);
            stats.addNode(hn, phg);

            // Push all neighbors of the added vertex into the queue
            for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
              data.bfs[queue_idx].add_pins_of_hyperedge_to_queue(
                he, phg, stats,_context.refinement.advanced.max_bfs_distance);
              stats.addEdge(he);
            }
          } else {
            release_vertex(search_id, hn);
          }
        }
      } else {
        // Note the associated refiner can lock a specific block. If
        // a block is locked, then the construction algorithm is not allowed
        // to add any vertex part of that block to the problem. In that case,
        // we clear all queues containing vertices of that block.
        data.clearBlock(block);
      }
    }
  }

  return nodes;
}

void AdvancedRefinementProblemConstruction::releaseNodes(const SearchID search_id,
                                                         const vec<HypernodeID>& nodes) {
  tbb::parallel_for(0UL, nodes.size(), [&](const size_t i) {
    release_vertex(search_id, nodes[i]);
  });
}

} // namespace mt_kahypar