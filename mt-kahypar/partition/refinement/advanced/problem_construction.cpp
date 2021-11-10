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


#include "mt-kahypar/partition/refinement/advanced/problem_construction.h"

#include <unordered_map>

#include "tbb/parallel_for.h"

namespace mt_kahypar {

void ProblemConstruction::BFSData::clearQueue(const PartitionID block) {
  ASSERT(blocks.i == block || blocks.j == block);
  const size_t idx = blocks.i == block ? 0 : 1;
  while ( !queue[idx].empty() ) queue[idx].pop();
  while ( !next_queue[idx].empty() ) next_queue[idx].pop();
}

void ProblemConstruction::BFSData::clearQueues() {
  current_distance = 0;
  last_queue_idx = 0;
  clearQueue(blocks.i);
  clearQueue(blocks.j);
}

void ProblemConstruction::BFSData::reset() {
  clearQueues();
  std::fill(visited_hn.begin(), visited_hn.end(), false);
  std::fill(visited_he.begin(), visited_he.end(), false);
}

void ProblemConstruction::BFSData::pop_hypernode(HypernodeID& hn) {
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

void ProblemConstruction::BFSData::add_pins_of_hyperedge_to_queue(
  const HyperedgeID& he,
  const PartitionedHypergraph& phg,
  const ProblemStats& stats,
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

namespace {
  using assert_map = std::unordered_map<HyperedgeID, bool>;
}

Subhypergraph ProblemConstruction::construct(const SearchID search_id,
                                             QuotientGraph& quotient_graph,
                                             AdvancedRefinerAdapter& refiner,
                                             const PartitionedHypergraph& phg) {
  Subhypergraph sub_hg;

  BFSData& bfs = _local_bfs.local();
  ProblemStats& stats = _local_stats.local();
  bfs.reset();
  stats.reset();
  bfs.blocks = quotient_graph.getBlockPair(search_id);
  stats.addBlock(bfs.blocks.i);
  stats.addBlock(bfs.blocks.j);


  // We initialize the BFS with all cut hyperedges running
  // between the involved block associated with the search
  bfs.clearQueues();
  quotient_graph.doForAllCutHyperedgesOfSearch(search_id, [&](const HyperedgeID& he) {
    bfs.add_pins_of_hyperedge_to_queue(he, phg, stats,
      _context.refinement.advanced.max_bfs_distance);
  });
  bfs.swap_with_next_queue();

  // BFS
  while ( !bfs.is_empty() &&
          !refiner.isMaximumProblemSizeReached(search_id, stats) ) {
    HypernodeID hn = kInvalidHypernode;
    bfs.pop_hypernode(hn);
    ASSERT(hn != kInvalidHypernode);

    PartitionID block = phg.partID(hn);
    if ( !stats.isLocked(block) ) {
      // Search aquires ownership of the vertex. Each vertex is only allowed to
      // be part of one search at any time in non-overlapping mode
      if ( _context.refinement.advanced.use_overlapping_searches || acquire_vertex(search_id, hn) ) {
        block = phg.partID(hn);
        // Double-check if vertex is still part of the blocks associated
        // with the search.
        if ( stats.isBlockContained(block) ) {
          if ( stats.block(0) == block ) {
            sub_hg.nodes_of_block_0.push_back(hn);
          } else {
            ASSERT(stats.block(1) == block);
            sub_hg.nodes_of_block_1.push_back(hn);
          }
          stats.addNode(hn, block, phg);

          // Push all neighbors of the added vertex into the queue
          for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
            bfs.add_pins_of_hyperedge_to_queue(
              he, phg, stats, _context.refinement.advanced.max_bfs_distance);
            stats.addEdge(he, sub_hg.hes);
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
      bfs.clearQueue(block);
    }

    if ( bfs.is_empty() ) {
      bfs.swap_with_next_queue();
    }
  }

  // Check if all touched hyperedges are contained in subhypergraph
  ASSERT([&]() {
    assert_map expected_hes;
    for ( const HyperedgeID& he : sub_hg.hes ) {
      if ( expected_hes.count(he) > 0 ) {
        LOG << "Hyperedge" << he << "is contained multiple times in subhypergraph!";
        return false;
      }
      expected_hes[he] = true;
    }

    for ( const HypernodeID& hn : sub_hg.nodes_of_block_0 ) {
      for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
        if ( expected_hes.count(he) == 0 ) {
          LOG << "Hyperedge" << he << "not contained in subhypergraph!";
          return false;
        }
        expected_hes[he] = false;
      }
    }

    for ( const HypernodeID& hn : sub_hg.nodes_of_block_1 ) {
      for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
        if ( expected_hes.count(he) == 0 ) {
          LOG << "Hyperedge" << he << "not contained in subhypergraph!";
          return false;
        }
        expected_hes[he] = false;
      }
    }

    for ( const auto& entry : expected_hes ) {
      const HyperedgeID he = entry.first;
      const bool visited = !entry.second;
      if ( !visited ) {
        LOG << "HyperedgeID" << he << "should be not part of subhypergraph!";
        return false;
      }
    }
    return true;
  }(), "Subhypergraph construction failed!");

  DBG << "Search ID =" << search_id
      << "-" << stats;

  const PartitionID block_0 = stats.block(0);
  const PartitionID block_1 = stats.block(1);
  sub_hg.weight_of_block_0 = stats.nodeWeightOfBlock(block_0);
  sub_hg.weight_of_block_1 = stats.nodeWeightOfBlock(block_1);
  sub_hg.num_pins = stats.numPins();
  return sub_hg;
}

void ProblemConstruction::releaseNodes(const SearchID search_id,
                                       const Subhypergraph& sub_hg) {
  if ( !_context.refinement.advanced.use_overlapping_searches ) {
    for ( size_t i = 0; i < sub_hg.nodes_of_block_0.size(); ++i ) {
      release_vertex(search_id, sub_hg.nodes_of_block_0[i]);
    }
    for ( size_t i = 0; i < sub_hg.nodes_of_block_1.size(); ++i ) {
      release_vertex(search_id, sub_hg.nodes_of_block_1[i]);
    }
  }
}

} // namespace mt_kahypar