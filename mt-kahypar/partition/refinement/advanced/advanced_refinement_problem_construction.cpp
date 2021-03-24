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
  visited_hn.clear();
  visited_he.clear();
}

void AdvancedRefinementProblemConstruction::BFSData::pop_hypernode(HypernodeID& hn) {
  if ( !is_empty() ) {
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
  const size_t max_bfs_distance) {
  if ( current_distance <= max_bfs_distance ) {
    if ( !visited_he.contains(he) ) {
      for ( const HypernodeID& pin : phg.pins(he) ) {
        const PartitionID block = phg.partID(pin);
        if ( (blocks.i == block || blocks.j == block) &&
              !visited_hn.contains(pin) ) {
          next_queue[blocks.i == block ? 0 : 1].push(pin);
          visited_hn[pin] = ds::EmptyStruct { };
        }
      }
      visited_he[he] = ds::EmptyStruct { };
    }
  }
}

void AdvancedRefinementProblemConstruction::ConstructionData::initialize(
  const vec<BlockPairCutHyperedges>& initial_cut_hes,
  const PartitionedHypergraph& phg) {
  used_slots = initial_cut_hes.size();
  while ( bfs.size() <= used_slots ) bfs.emplace_back();

  // Initialize BFS Queues
  for ( size_t i = 0; i < used_slots; ++i ) {
    bfs[i].reset();
    bfs[i].blocks = initial_cut_hes[i].blocks;
    for ( const HyperedgeID& he : initial_cut_hes[i].cut_hes ) {
      bfs[i].add_pins_of_hyperedge_to_queue(
        he, phg, std::numeric_limits<size_t>::max());
    }
    bfs[i].swap_with_next_queue();
  }
}

void AdvancedRefinementProblemConstruction::ConstructionData::pop_hypernode(HypernodeID& hn, size_t& idx) {
  if ( !is_empty() ) {
    bool found = false;
    idx = last_idx++ % used_slots;
    while ( !found ) {
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

void AdvancedRefinementProblemConstruction::ProblemData::add_hypernode(
  const HypernodeID hn,
  const PartitionedHypergraph& phg) {
  const PartitionID block = phg.partID(hn);
  if ( used_blocks[block] == std::numeric_limits<size_t>::max() ) {
    used_blocks[block] = stats.used_blocks.size();
    stats.used_blocks.push_back(block);
    stats.num_nodes_in_blocks.push_back(0);
  }
  size_t idx = used_blocks[block];
  ++stats.num_nodes_in_blocks[idx];
  stats.num_pins += phg.nodeDegree(hn);
}

void AdvancedRefinementProblemConstruction::ProblemData::reset() {
  stats.num_nodes_in_blocks.clear();
  stats.used_blocks.clear();
  stats.num_edges = 0;
  stats.num_pins = 0;
  used_blocks.assign(used_blocks.size(), std::numeric_limits<size_t>::max());
  visited_hes.clear();
}

vec<HypernodeID> AdvancedRefinementProblemConstruction::construct(const SearchID search_id,
                                                                  QuotientGraph& quotient_graph,
                                                                  AdvancedRefinerAdapter& refiner,
                                                                  const PartitionedHypergraph& phg) {
  vec<HypernodeID> nodes;

  // Initialize BFS
  ProblemData& problem = _local_problem.local();
  ConstructionData& data = _local_data.local();
  problem.reset();
  const size_t num_block_pairs = quotient_graph.numBlockPairs(search_id);

  while ( !refiner.isMaximumProblemSizeReached(search_id, problem.stats) ) {
    const vec<BlockPairCutHyperedges> initial_cut_hes =
      quotient_graph.requestCutHyperedges(search_id, num_block_pairs *
        _context.refinement.advanced.num_cut_edges_per_block_pair);
    data.initialize(initial_cut_hes, phg);
    if ( data.is_empty() ) {
      break;
    }

    while ( !data.is_empty() &&
            !refiner.isMaximumProblemSizeReached(search_id, problem.stats) ) {
      size_t queue_idx = std::numeric_limits<size_t>::max();
      HypernodeID hn = kInvalidHypernode;
      data.pop_hypernode(hn, queue_idx);
      ASSERT(hn != kInvalidHypernode);
      ASSERT(queue_idx != std::numeric_limits<size_t>::max());

      if ( acquire_vertex(search_id, hn) ) {
        nodes.push_back(hn);
        problem.add_hypernode(hn, phg);

        for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
          data.bfs[queue_idx].add_pins_of_hyperedge_to_queue(
            he, phg, _context.refinement.advanced.max_bfs_distance);
          problem.add_hyperedge(he);
        }
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