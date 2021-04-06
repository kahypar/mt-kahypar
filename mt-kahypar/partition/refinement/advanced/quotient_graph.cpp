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

SearchID QuotientGraph::requestNewSearch(AdvancedRefinerAdapter& refiner) {
  ASSERT(_phg);
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
      _searches.emplace_back();
      _searches.back().addBlockPair(blocks);

      // Associate refiner with search id
      const bool success = refiner.registerNewSearch(search_id, *_phg);
      ASSERT(success); unused(success);

      // Extend with additional
      const PartitionID num_blocks = refiner.maxNumberOfBlocks(search_id);
      if ( num_blocks > 2 ) {
        vec<BlockPair> additional_blocks;
        if ( isOneBlockUnderloaded(blocks.i, blocks.j) ) {
          // Find Path including edge (blocks.i, blocks.j)
          const PartitionID overloaded_block =
            _phg->partWeight(blocks.i) >= _context.partition.perfect_balance_part_weights[blocks.i] ?
              blocks.i : blocks.j;
          const PartitionID underloaded_block = overloaded_block == blocks.i ? blocks.j : blocks.i;
          additional_blocks = extendWithPath(overloaded_block, underloaded_block, num_blocks - 2);
        } else {
          // Find Cycle including edge (blocks.i, blocks.j)
          additional_blocks = extendWithCycle(blocks, num_blocks - 2);

          // Verify that the additional blocks form a cycle
          ASSERT([&]() {
            if ( additional_blocks.size() > 0 ) {
              vec<BlockPair> cycle = { blocks };
              for ( const BlockPair& add_block : additional_blocks ) {
                cycle.push_back(add_block);
              }
              if ( !verifyCycle(cycle) ) {
                for ( const BlockPair& blocks : cycle ) {
                  LOG << V(blocks.i) << V(blocks.j);
                }
                return false;
              }
            }
            return true;
          }(), "Additional block pairs do not form a cycle!");
        }

        // Add additional blocks to search and remove them from block scheduler heap
        for ( const BlockPair& additional_block : additional_blocks ) {
          _searches.back().addBlockPair(additional_block);
          const size_t part_index =
            index(_quotient_graph[additional_block.i][additional_block.j]);
          ASSERT(_block_scheduler.contains(part_index));
          _block_scheduler.remove(part_index);
        }
      }

      // Increment number of active searches for each edge in the quotient graph
      // that contains one of the used blocks in the search
      auto increment_num_active_searches = [&](const PartitionID i, const PartitionID j) {
        ++_quotient_graph[i][j].stats.num_active_searches;
        const size_t part_index = index(_quotient_graph[i][j]);
        if ( _block_scheduler.contains(part_index) ) {
          const BlockPairStats tmp_stats = _quotient_graph[i][j].stats;
          _block_scheduler.updateKey(part_index, tmp_stats);
        }
      };

      vec<bool> used_blocks(_context.partition.k, false);
      for ( const BlockPair& blocks : _searches.back().block_pairs ) {
        used_blocks[blocks.i] = true;
        used_blocks[blocks.j] = true;
        increment_num_active_searches(blocks.i, blocks.j);
      }
      for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
        for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
          if ( used_blocks[i] || used_blocks[j] ) {
            increment_num_active_searches(i, j);
          }
        }
      }
    }
    _heap_lock.unlock();
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
      if ( qg_edge.is_active() ) {
        const HyperedgeID he = qg_edge.pop_hyperedge();
        qg_edge.stats.cut_he_weight -= _phg->edgeWeight(he);
        // Note, we only consider hyperedges that contains pins of both blocks.
        // There might be some edges which were initially cut, but were removed due
        // to vertex moves. Thus, we remove them lazily here.
        if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
             _phg->pinCountInPart(he, blocks.j) > 0 ) {
          block_pair_cut_hes[i].cut_hes.push_back(he);
          _searches[search_id].used_cut_hes[i].push_back(he);
          ++num_edges;
        }
      } else {
        --num_active_block_pairs;
      }
      ++pos;
    }
  }
  return block_pair_cut_hes;
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
  for ( const BlockPair& blocks : _searches[search_id].block_pairs ) {
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
}

void QuotientGraph::finalizeSearch(const SearchID search_id,
                                   const bool success) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  ASSERT(_searches[search_id].is_finalized);
  // Decrement number of active searches for each edge in the quotient graph
  // that contains one of the used blocks in the search
  vec<bool> used_blocks(_context.partition.k, false);
  for ( const BlockPair& blocks : _searches[search_id].block_pairs ) {
    used_blocks[blocks.i] = true;
    used_blocks[blocks.j] = true;
    ASSERT(_quotient_graph[blocks.i][blocks.j].stats.num_active_searches > 0);
    --_quotient_graph[blocks.i][blocks.j].stats.num_active_searches;
  }
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      if ( used_blocks[i] || used_blocks[j] ) {
        ASSERT(_quotient_graph[i][j].stats.num_active_searches > 0, V(i) << V(j));
        --_quotient_graph[i][j].stats.num_active_searches;
      }
    }
  }

  if ( success ) {
    // If the search improves the quality of the partition, we reinsert
    // all hyperedges that were used by the search and are still cut.
    ASSERT(_searches[search_id].block_pairs.size() == _searches[search_id].used_cut_hes.size());
    for ( size_t i = 0; i < _searches[search_id].block_pairs.size(); ++i ) {
      const BlockPair& blocks = _searches[search_id].block_pairs[i];
      QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
      for ( const HyperedgeID& he : _searches[search_id].used_cut_hes[i] ) {
        if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
            _phg->pinCountInPart(he, blocks.j) > 0 ) {
          qg_edge.add_hyperedge(he, _phg->edgeWeight(he));
        }
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
  if ( _context.refinement.advanced.sort_cut_hes ) {
    tbb::parallel_for(0, _context.partition.k, [&](const PartitionID i) {
      tbb::parallel_for(i + 1, _context.partition.k, [&, i](const PartitionID j) {
        BFSData& bfs_data = _local_bfs.local();
        sortCutHyperedges(i, j, bfs_data);
      });
    });
  }
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

vec<BlockPair> QuotientGraph::extendWithPath(const PartitionID overloaded_block,
                                             const PartitionID underloaded_block,
                                             const size_t num_additional_blocks) {
  vec<BlockPair> best_path;
  double best_rating = 0.0;
  vec<BlockPair> current_path;
  double current_rating = 0.0;
  std::function<void(const PartitionID)> dfs = [&](const PartitionID i) {
    // We find to find a path that includes as much
    // cut hyperedges as possible
    if ( current_rating > best_rating ) {
      best_path = current_path;
      best_rating = current_rating;
    }

    if ( current_path.size() == num_additional_blocks ) return;
    for ( PartitionID j = 0; j < _context.partition.k; ++j ) {
      if ( i != j && j != underloaded_block && j != overloaded_block &&
           _phg->partWeight(j) >= _context.partition.perfect_balance_part_weights[j] ) {
        const QuotientGraphEdge& qg_edge = _quotient_graph[std::min(i,j)][std::max(i,j)];
        const size_t part_index = index(qg_edge);
        // Perform DFS on all blocks that are currently contained in
        // the block scheduler heap
        if ( _block_scheduler.contains(part_index) ) {
          current_path.push_back(BlockPair { std::min(i,j), std::max(i,j) });
          // Note, weight of all cut hyperedges could change, since it can be
          // modified concurrently
          const double rating =
            static_cast<double>(qg_edge.stats.cut_he_weight) /
              std::max(qg_edge.stats.num_active_searches.load(std::memory_order_relaxed), 1);
          current_rating += rating;

          // Continue DFS on block j
          dfs(j);

          current_path.pop_back();
          current_rating -= rating;
        }
      }
    }
  };

  // Start DFS from the overloaded and underloaded block
  dfs(overloaded_block);
  dfs(underloaded_block);
  return best_path;
}

vec<BlockPair> QuotientGraph::extendWithCycle(const BlockPair initial_blocks,
                                              const size_t num_additional_blocks) {
  vec<BlockPair> best_cycle;
  double best_rating = 0.0;
  vec<BlockPair> current_cycle;
  double current_rating = 0.0;
  std::function<void(const PartitionID,
                     const PartitionID,
                     const PartitionID)>
    dfs = [&](const PartitionID start,
              const PartitionID end,
              const PartitionID i) {
    // We find to find a path that includes as much
    // cut hyperedges as possible
    if ( i == end && current_rating > best_rating ) {
      best_cycle = current_cycle;
      best_rating = current_rating;
      return;
    } else if ( current_cycle.size() == num_additional_blocks + 1 ) return;

    for ( PartitionID j = 0; j < _context.partition.k; ++j ) {
      if ( i != j && j != start && ( i != start || j != end ) ) {
        const QuotientGraphEdge& qg_edge = _quotient_graph[std::min(i,j)][std::max(i,j)];
        const size_t part_index = index(qg_edge);
        // Perform DFS on all blocks that are currently contained in
        // the block scheduler heap
        if ( _block_scheduler.contains(part_index) ) {
          current_cycle.push_back(BlockPair { std::min(i,j), std::max(i,j) });
          // Note, weight of all cut hyperedges could change, since it can be
          // modified concurrently
          const double rating =
            static_cast<double>(qg_edge.stats.cut_he_weight) /
              std::max(qg_edge.stats.num_active_searches.load(std::memory_order_relaxed), 1);
          current_rating += rating;

          // Continue DFS on block j
          dfs(start, end, j);

          current_cycle.pop_back();
          current_rating -= rating;
        }
      }
    }
  };

  // Find cycle from initial_blocks.i to initial_blocks.j
  dfs(initial_blocks.i, initial_blocks.j, initial_blocks.i);
  // Find cycle from initial_blocks.j to initial_blocks.i
  dfs(initial_blocks.j, initial_blocks.i, initial_blocks.j);

  return best_cycle;
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
  for ( const HyperedgeID& he : _quotient_graph[i][j].cut_hes ) {
    if ( bfs_data.distance[he] == -1 ) {
      bfs(he);
    }
  }

  // Sort all cut hyperedges according to their distance label
  std::sort(_quotient_graph[i][j].cut_hes.begin(), _quotient_graph[i][j].cut_hes.end(),
    [&](const HyperedgeID& lhs, const HyperedgeID& rhs) {
      ASSERT(bfs_data.distance[lhs] != -1);
      ASSERT(bfs_data.distance[rhs] != -1);
      const int distance_lhs = bfs_data.distance[lhs];
      const int distance_rhs = bfs_data.distance[rhs];
      return distance_lhs < distance_rhs || (distance_lhs == distance_rhs && lhs < rhs);
    });
}

} // namespace mt_kahypar