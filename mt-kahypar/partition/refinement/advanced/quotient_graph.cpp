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
  return lhs.num_active_searches == rhs.num_active_searches &&
    lhs.is_one_block_underloaded == rhs.is_one_block_underloaded &&
    lhs.cut_he_weight == rhs.cut_he_weight;
}

std::ostream& operator<<(std::ostream& out, const BlockPairStats& stats) {
  out << "Num Active Searches = " << stats.num_active_searches << ", "
      << "Is One Block Underloaded = " << std::boolalpha << stats.is_one_block_underloaded << ", "
      << "Cut Hyperedge Weight = " << stats.cut_he_weight;
  return out;
}

void QuotientGraph::QuotientGraphEdge::add_hyperedge(const HyperedgeID he,
                                                     const HyperedgeWeight weight) {
  cut_hes.push_back(he);
  stats.cut_he_weight += weight;
}

HyperedgeID QuotientGraph::QuotientGraphEdge::pop_hyperedge() {
  ASSERT(is_active());
  return cut_hes[first_valid_entry++];
}

void QuotientGraph::QuotientGraphEdge::reset(const bool is_one_block_underloaded) {
  cut_hes.clear();
  first_valid_entry = 0;
  stats.num_active_searches.store(0, std::memory_order_relaxed);
  stats.is_one_block_underloaded = is_one_block_underloaded;
  stats.cut_he_weight.store(0, std::memory_order_relaxed);
}

SearchID QuotientGraph::requestNewSearch() {
  SearchID search_id = INVALID_SEARCH_ID;
  if ( !_block_scheduler.empty() ) {
    _heap_lock.lock();
    if ( !_block_scheduler.empty() ) {
      ++_num_active_searches;
      // Retrieve block pair for next search
      const size_t part_index = _block_scheduler.top();
      const BlockPair blocks = blockPairFromIndex(part_index);
      _block_scheduler.pop();
      // Create new search
      search_id = _searches.size();
      _searches.emplace_back(blocks);

      // Increment number of active searches for each edge in the quotient graph
      // that contains either block blocks.i or blocks.j
      auto increment_num_active_searches = [&](const PartitionID i, const PartitionID j) {
        ++_quotient_graph[i][j].stats.num_active_searches;
        const size_t part_index = index(_quotient_graph[i][j]);
        if ( _block_scheduler.contains(part_index) ) {
          const BlockPairStats tmp_stats = _quotient_graph[i][j].stats;
          _block_scheduler.updateKey(part_index, tmp_stats);
        }
      };
      for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
        for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
          if ( i == blocks.i || j == blocks.i || i == blocks.j || j == blocks.j ) {
            increment_num_active_searches(i, j);
          }
        }
      }
      increment_num_active_searches(blocks.i, blocks.j);
    }
    _heap_lock.unlock();
  }
  return search_id;
}

vec<HyperedgeID> QuotientGraph::requestCutHyperedges(const SearchID search_id,
                                                     const size_t max_num_edges) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  vec<HyperedgeID> cut_hes;
  if ( !_searches[search_id].is_finalized ) {
    BlockPair blocks = _searches[search_id].blocks;
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    while ( qg_edge.is_active() &&
            cut_hes.size() < max_num_edges ) {
      const HyperedgeID he = qg_edge.pop_hyperedge();
      qg_edge.stats.cut_he_weight -= _phg->edgeWeight(he);
      // Note, we only consider hyperedges that contains pins of both blocks.
      // There might be some edges which were initially cut, but were removed due
      // to vertex moves. Thus, we remove them lazily here.
      if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
          _phg->pinCountInPart(he, blocks.j) > 0 ) {
        cut_hes.push_back(he);
        _searches[search_id].used_cut_hes.push_back(he);
      }
    }
  }
  return cut_hes;
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

void QuotientGraph::finalizeConstruction(const SearchID search_id) {
  ASSERT(search_id < _searches.size());
  _searches[search_id].is_finalized = true;
  BlockPair blocks = _searches[search_id].blocks;
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  if ( qg_edge.is_active() ) {
    _heap_lock.lock();
    const size_t part_index = index(qg_edge);
    ASSERT(!_block_scheduler.contains(part_index));
    const BlockPairStats tmp_stats = qg_edge.stats;
    _block_scheduler.push(part_index, tmp_stats);
    _heap_lock.unlock();
  }
}

void QuotientGraph::finalizeSearch(const SearchID search_id,
                                   const bool success) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  ASSERT(_searches[search_id].is_finalized);
  BlockPair blocks = _searches[search_id].blocks;
  // Decrement number of active searches for each edge in the quotient graph
  // that contains either block blocks.i or blocks.j
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      if ( i == blocks.i || j == blocks.i || i == blocks.j || j == blocks.j ) {
        ASSERT(_quotient_graph[i][j].stats.num_active_searches > 0);
        --_quotient_graph[i][j].stats.num_active_searches;
      }
    }
  }
  ASSERT(_quotient_graph[blocks.i][blocks.j].stats.num_active_searches > 0);
  --_quotient_graph[blocks.i][blocks.j].stats.num_active_searches;

  if ( success ) {
    // If the search improves the quality of the partition, we reinsert
    // all hyperedges that were used by the search and are still cut.
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    for ( const HyperedgeID& he : _searches[search_id].used_cut_hes ) {
      if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
           _phg->pinCountInPart(he, blocks.j) > 0 ) {
        qg_edge.add_hyperedge(he, _phg->edgeWeight(he));
      }
    }
  }

  // Update all entries of the heap
  fullHeapUpdate();
  --_num_active_searches;
}

void QuotientGraph::initialize(const PartitionedHypergraph& phg) {
  _phg = &phg;

  // Reset internal members
  resetQuotientGraphEdges();
  _block_scheduler.clear();
  _num_active_searches.store(0, std::memory_order_relaxed);
  _searches.clear();

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

  // Initalize block scheduler heap
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      if ( _quotient_graph[i][j].is_active() ) {
        _block_scheduler.push(
          index(_quotient_graph[i][j]), _quotient_graph[i][j].stats);
      }
    }
  }

  // Sort cut hyperedges of each block
  tbb::parallel_for(0, _context.partition.k, [&](const PartitionID i) {
    tbb::parallel_for(i + 1, _context.partition.k, [&, i](const PartitionID j) {
      sortCutHyperedges(i, j);
    });
  });
}

void QuotientGraph::fullHeapUpdate() {
  ASSERT(_phg);
  _heap_lock.lock();
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      QuotientGraphEdge& qg_edge = _quotient_graph[i][j];
      const size_t part_index = index(qg_edge);
      qg_edge.stats.is_one_block_underloaded = isOneBlockUnderloaded(i, j);
      const BlockPairStats tmp_stats = qg_edge.stats;
      if ( _block_scheduler.contains(part_index) ) {
        _block_scheduler.updateKey(part_index, tmp_stats);
      } else if ( qg_edge.is_active() && qg_edge.stats.num_active_searches == 0 ) {
        _block_scheduler.push(part_index, tmp_stats);
      }
    }
  }
  _heap_lock.unlock();
}

void QuotientGraph::resetQuotientGraphEdges() {
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].reset(isOneBlockUnderloaded(i, j));
    }
  }
}

void QuotientGraph::sortCutHyperedges(const PartitionID i, const PartitionID j) {
  ASSERT(_phg);
  ASSERT(i < j);
  ASSERT(0 <= i && i < _context.partition.k);
  ASSERT(0 <= j && j < _context.partition.k);
  ds::DynamicSparseSet<HypernodeID> visited_hns;
  ds::DynamicSparseMap<HyperedgeID, int> distance;
  int current_distance = 0;

  // BFS that traverses all hyperedges reachable from the
  // corresponding seed hyperedge
  auto bfs = [&](const HyperedgeID seed) {
    std::queue<HyperedgeID> q;
    std::queue<HyperedgeID> next_q;
    q.push(seed);
    ASSERT(_phg->pinCountInPart(seed, i) > 0 && _phg->pinCountInPart(seed, j) > 0);
    ASSERT(!distance.contains(seed));
    distance[seed] = current_distance;

    while ( !q.empty() ) {
      const HyperedgeID he = q.front();
      q.pop();

      for ( const HypernodeID& pin : _phg->pins(he) ) {
        const PartitionID block = _phg->partID(pin);
        if ( ( block == i || block == j ) && !visited_hns.contains(pin) ) {
          for ( const HyperedgeID& inc_he : _phg->incidentEdges(pin) ) {
            if ( !distance.contains(inc_he) ) {
              next_q.push(inc_he);
              distance[inc_he] = current_distance + 1;
            }
          }
          visited_hns[pin] = ds::EmptyStruct { };
        }
      }

      if ( q.empty() ) {
        q.swap(next_q);
        ++current_distance;
      }
    }
  };

  // Start BFS
  for ( const HyperedgeID& he : _quotient_graph[i][j].cut_hes ) {
    if ( !distance.contains(he) ) {
      bfs(he);
    }
  }

  // Sort all cut hyperedges according to their distance label
  tbb::parallel_sort(_quotient_graph[i][j].cut_hes.begin(), _quotient_graph[i][j].cut_hes.end(),
    [&](const HyperedgeID& lhs, const HyperedgeID& rhs) {
      ASSERT(distance.contains(lhs));
      ASSERT(distance.contains(rhs));
      const int distance_lhs = distance[lhs];
      const int distance_rhs = distance[rhs];
      return distance_lhs < distance_rhs || (distance_lhs == distance_rhs && lhs < rhs);
    });
}

} // namespace mt_kahypar