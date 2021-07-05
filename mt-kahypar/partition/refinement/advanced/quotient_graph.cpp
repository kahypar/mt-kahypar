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


#include "mt-kahypar/partition/refinement/advanced/quotient_graph.h"

#include <queue>

#include "tbb/parallel_sort.h"

#include "mt-kahypar/datastructures/sparse_map.h"


namespace mt_kahypar {

bool operator==(const BlockPairStats& lhs, const BlockPairStats& rhs) {
  return lhs.initial_cut_he_weight == rhs.initial_cut_he_weight &&
         lhs.num_acitve_searches_on_block_pair == rhs.num_acitve_searches_on_block_pair &&
         lhs.max_active_searches_on_one_block == rhs.max_active_searches_on_one_block;
}

std::ostream& operator<<(std::ostream& out, const BlockPairStats& stats) {
  out << "Initial Cut Hyperedge Weight = " << stats.initial_cut_he_weight << ", "
      << "Num Active Searches on Block Pair = " << stats.num_acitve_searches_on_block_pair << ", "
      << "Max Active Searches on One Block = " << stats.max_active_searches_on_one_block;
  return out;
}

void QuotientGraph::QuotientGraphEdge::add_hyperedge(const HyperedgeID he,
                                                     const HyperedgeWeight weight) {
  cut_hes.push_back(CutHyperedge { he, CAtomic<bool>(false) });
  cut_he_weight += weight;
}

QuotientGraph::CutHyperedge& QuotientGraph::QuotientGraphEdge::pop_hyperedge() {
  ASSERT(isActive());
  return cut_hes[first_valid_entry++];
}

void QuotientGraph::QuotientGraphEdge::reset() {
  cut_hes.clear();
  ownership.store(INVALID_SEARCH_ID, std::memory_order_relaxed);
  is_in_queue.store(false, std::memory_order_relaxed);
  first_valid_entry = 0;
  stats.initial_cut_he_weight = 0;
  stats.num_acitve_searches_on_block_pair.store(0, std::memory_order_relaxed);
  stats.max_active_searches_on_one_block = 0;
  cut_he_weight.store(0, std::memory_order_relaxed);
  round.store(0, std::memory_order_relaxed);
}

SearchID QuotientGraph::requestNewSearch(AdvancedRefinerAdapter& refiner) {
  ASSERT(_phg);
  SearchID search_id = INVALID_SEARCH_ID;
  if ( !_block_scheduler.empty() ) {
    BlockPair blocks { kInvalidPartition, kInvalidPartition };
    bool success = popBlockPairFromQueue(blocks);
    _register_search_lock.lock();
    const SearchID tmp_search_id = _searches.size();
    if ( success && _quotient_graph[blocks.i][blocks.j].acquire(tmp_search_id) ) {
      ++_num_active_searches;
      // Create new search
      search_id = tmp_search_id;
      _searches.emplace_back();
      _searches[search_id].addBlockPair(blocks);
      _register_search_lock.unlock();

      // Associate refiner with search id
      success = refiner.registerNewSearch(search_id, *_phg);
      ASSERT(success); unused(success);

      // TODO: extend with additional block pairs if search requests
    } else {
      _register_search_lock.unlock();
    }
  }
  return search_id;
}

vec<BlockPairCutHyperedges> QuotientGraph::requestCutHyperedges(const SearchID search_id,
                                                                const size_t max_num_edges) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  vec<BlockPairCutHyperedges> block_pair_cut_hes;
  if ( !_searches[search_id].is_finalized ) {
    Search& search = _searches[search_id];
    for ( const BlockPair& blocks : search.block_pairs ) {
      block_pair_cut_hes.emplace_back();
      block_pair_cut_hes.back().blocks = blocks;
    }

    // Add cut hyperedges in round robin fashion
    size_t num_block_pairs = block_pair_cut_hes.size();
    size_t num_active_block_pairs = num_block_pairs;
    size_t pos = 0;
    size_t num_edges = 0;
    while ( num_active_block_pairs > 0 && num_edges < max_num_edges ) {
      const size_t i = pos % num_block_pairs;
      if ( i == 0 ) num_active_block_pairs = num_block_pairs;
      const BlockPair& blocks = block_pair_cut_hes[i].blocks;
      QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
      if ( qg_edge.isActive() ) {
        CutHyperedge& cut_he = qg_edge.pop_hyperedge();
        if ( cut_he.acquire() ) {
          // Note, we only consider hyperedges that contains pins of both blocks.
          // There might be some edges which were initially cut, but were removed due
          // to vertex moves. Thus, we remove them lazily here.
          const HyperedgeID he = cut_he.he;
          qg_edge.cut_he_weight -= _phg->edgeWeight(he);
          if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
              _phg->pinCountInPart(he, blocks.j) > 0 ) {
            block_pair_cut_hes[i].cut_hes.push_back(he);
            _searches[search_id].used_cut_hes[i].push_back(he);
            ++num_edges;
          }
        }
      } else {
        --num_active_block_pairs;
      }
      ++pos;
    }
  }
  return block_pair_cut_hes;
}

bool QuotientGraph::popBlockPairFromQueue(BlockPair& blocks) {
  blocks.i = kInvalidPartition;
  blocks.j = kInvalidPartition;
  while ( _block_scheduler.try_pop(blocks) ) {
    _quotient_graph[blocks.i][blocks.j].markAsNotInQueue();
    const bool not_too_many_concurrent_searches =
      _num_active_searches_on_blocks[blocks.i] < _context.refinement.advanced.max_concurrency_per_block &&
      _num_active_searches_on_blocks[blocks.j] < _context.refinement.advanced.max_concurrency_per_block;
    if ( !not_too_many_concurrent_searches ) {
      pushBlockPairIntoQueue(blocks);
      blocks.i = kInvalidPartition;
      blocks.j = kInvalidPartition;
    } else {
      break;
    }
  }
  return blocks.i != kInvalidPartition && blocks.j != kInvalidPartition;
}

bool QuotientGraph::pushBlockPairIntoQueue(const BlockPair& blocks) {
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  const bool is_promising_block_pair =
        !_context.refinement.advanced.skip_unpromising_blocks ||
        ( qg_edge.round == 0 || qg_edge.num_improvements > 0 );
  if ( qg_edge.isActive() && is_promising_block_pair && qg_edge.markAsInQueue() ) {
    _block_scheduler.push(blocks);
    return true;
  } else {
    return false;
  }
}

void QuotientGraph::addNewCutHyperedge(const HyperedgeID he,
                                       const PartitionID block) {
  ASSERT(_phg);
  ASSERT(_phg->pinCountInPart(he, block) > 0);
  // Add hyperedge he as a cut hyperedge to each block pair that contains 'block'
  for ( const PartitionID& other_block : _phg->connectivitySet(he) ) {
    if ( other_block != block ) {
      _quotient_graph[std::min(block, other_block)][std::max(block, other_block)]
        .add_hyperedge(he, _phg->edgeWeight(he));
    }
  }
}

size_t QuotientGraph::acquireUsedCutHyperedges(const SearchID& search_id, const vec<bool>& used_hes) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  size_t additional_cut_hes = 0;
  for ( size_t i = 0; i < _searches[search_id].block_pairs.size(); ++i ) {
    const BlockPair& blocks = _searches[search_id].block_pairs[i];
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    const size_t start_idx = qg_edge.first_valid_entry;
    const size_t end_idx = qg_edge.cut_hes.size();
    for ( size_t j = start_idx; j < end_idx; ++j ) {
      CutHyperedge& cut_he = qg_edge.cut_hes[j];
      if ( used_hes[cut_he.he] && cut_he.acquire() ) {
        qg_edge.cut_he_weight -= _phg->edgeWeight(cut_he.he);
        _searches[search_id].used_cut_hes[i].push_back(cut_he.he);
        ++additional_cut_hes;
      }
    }
  }
  return additional_cut_hes;
}

void QuotientGraph::finalizeConstruction(const SearchID search_id) {
  ASSERT(search_id < _searches.size());
  _searches[search_id].is_finalized = true;
  for ( const BlockPair& blocks : _searches[search_id].block_pairs ) {
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    qg_edge.release(search_id);
    pushBlockPairIntoQueue(blocks);
  }
}

void QuotientGraph::finalizeSearch(const SearchID search_id,
                                   const bool success) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  ASSERT(_searches[search_id].is_finalized);

  ASSERT(_searches[search_id].block_pairs.size() == _searches[search_id].used_cut_hes.size());
  for ( size_t i = 0; i < _searches[search_id].block_pairs.size(); ++i ) {
    const BlockPair& blocks = _searches[search_id].block_pairs[i];
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    ++qg_edge.round;
    if ( success ) {
      // If the search improves the quality of the partition, we reinsert
      // all hyperedges that were used by the search and are still cut.
      ++qg_edge.num_improvements;
      for ( const HyperedgeID& he : _searches[search_id].used_cut_hes[i] ) {
        if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
            _phg->pinCountInPart(he, blocks.j) > 0 ) {
          qg_edge.add_hyperedge(he, _phg->edgeWeight(he));
        }
      }
    }
    // In case the block pair becomes active,
    // we reinsert it into the queue
    pushBlockPairIntoQueue(blocks);
  }
  --_num_active_searches;
}

void QuotientGraph::initialize(const PartitionedHypergraph& phg) {
  _phg = &phg;

  // Reset internal members
  resetQuotientGraphEdges(phg);
  _block_scheduler.clear();
  _num_active_searches.store(0, std::memory_order_relaxed);
  _searches.clear();
  _num_active_searches_on_blocks.assign(
    _context.partition.k, CAtomic<size_t>(0));

  // Find all cut hyperedges between the blocks
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
    const HyperedgeWeight edge_weight = phg.edgeWeight(he);
    for ( const PartitionID i : phg.connectivitySet(he) ) {
      for ( const PartitionID j : phg.connectivitySet(he) ) {
        if ( i < j ) {
          _quotient_graph[i][j].add_hyperedge(he, edge_weight);
        }
      }
    }
  });

  // Initalize block scheduler queue
  std::vector<BlockPair> active_blocks;
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].stats.initial_cut_he_weight = _quotient_graph[i][j].cut_he_weight;
      if ( _quotient_graph[i][j].isActive() ) {
        active_blocks.push_back( BlockPair { i, j } );
      }
    }
  }
  std::sort(active_blocks.begin(), active_blocks.end(),
    [&](const BlockPair& lhs, const BlockPair& rhs) {
      return _quotient_graph[lhs.i][lhs.j].stats.initial_cut_he_weight >
        _quotient_graph[rhs.i][rhs.j].stats.initial_cut_he_weight;
    });
  for ( const BlockPair& blocks : active_blocks ) {
    pushBlockPairIntoQueue(blocks);
  }

  // Sort cut hyperedges of each block
  if ( _context.refinement.advanced.sort_cut_hes ) {
    tbb::parallel_for(0, _context.partition.k, [&](const PartitionID i) {
      tbb::parallel_for(i + 1, _context.partition.k, [&, i](const PartitionID j) {
        BFSData& bfs_data = _local_bfs.local();
        sortCutHyperedges(i, j, bfs_data);
      });
    });
  }
}

void QuotientGraph::resetQuotientGraphEdges(const PartitionedHypergraph& phg) {
  const bool skip_small_cuts = !isInputHypergraph(phg) && _context.refinement.advanced.skip_small_cuts;
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].reset();
      _quotient_graph[i][j].skip_small_cuts = skip_small_cuts;
    }
  }
}

void QuotientGraph::sortCutHyperedges(const PartitionID i,
                                      const PartitionID j,
                                      BFSData& bfs_data) {
  ASSERT(_phg);
  ASSERT(i < j);
  ASSERT(0 <= i && i < _context.partition.k);
  ASSERT(0 <= j && j < _context.partition.k);
  bfs_data.reset();
  int current_distance = 0;

  // BFS that traverses all hyperedges reachable from the
  // corresponding seed hyperedge
  auto bfs = [&](const HyperedgeID seed) {
    std::queue<HyperedgeID> q;
    std::queue<HyperedgeID> next_q;
    q.push(seed);
    ASSERT(_phg->pinCountInPart(seed, i) > 0 && _phg->pinCountInPart(seed, j) > 0);
    bfs_data.distance[seed] = current_distance;

    while ( !q.empty() ) {
      const HyperedgeID he = q.front();
      q.pop();

      for ( const HypernodeID& pin : _phg->pins(he) ) {
        const PartitionID block = _phg->partID(pin);
        if ( ( block == i || block == j ) && !bfs_data.visited_hns[pin] ) {
          for ( const HyperedgeID& inc_he : _phg->incidentEdges(pin) ) {
            if ( bfs_data.distance[inc_he] == -1 ) {
              next_q.push(inc_he);
              bfs_data.distance[inc_he] = current_distance + 1;
            }
          }
          bfs_data.visited_hns[pin] = true;
        }
      }

      if ( q.empty() ) {
        q.swap(next_q);
        ++current_distance;
      }
    }
  };

  // Start BFS
  for ( const CutHyperedge& cut_he : _quotient_graph[i][j].cut_hes ) {
    const HyperedgeID he = cut_he.he;
    if ( bfs_data.distance[he] == -1 ) {
      bfs(he);
    }
  }

  // Sort all cut hyperedges according to their distance label
  std::sort(_quotient_graph[i][j].cut_hes.begin(), _quotient_graph[i][j].cut_hes.end(),
    [&](const CutHyperedge& lhs, const CutHyperedge& rhs) {
      ASSERT(bfs_data.distance[lhs.he] != -1);
      ASSERT(bfs_data.distance[rhs.he] != -1);
      const int distance_lhs = bfs_data.distance[lhs.he];
      const int distance_rhs = bfs_data.distance[rhs.he];
      return distance_lhs < distance_rhs || (distance_lhs == distance_rhs && lhs.he < rhs.he);
    });
}

} // namespace mt_kahypar