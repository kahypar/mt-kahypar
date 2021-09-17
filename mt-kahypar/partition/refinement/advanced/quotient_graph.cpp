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

void QuotientGraph::QuotientGraphEdge::add_hyperedge(const HyperedgeID he,
                                                     const HyperedgeWeight weight) {
  cut_hes.push_back(he);
  cut_he_weight += weight;
}

HyperedgeID QuotientGraph::QuotientGraphEdge::pop_hyperedge() {
  ASSERT(cut_he_weight > 0);
  return cut_hes[first_valid_entry++];
}

void QuotientGraph::QuotientGraphEdge::reset() {
  cut_hes.clear();
  ownership.store(INVALID_SEARCH_ID, std::memory_order_relaxed);
  is_in_queue.store(false, std::memory_order_relaxed);
  first_valid_entry = 0;
  initial_num_cut_hes = 0;
  initial_cut_he_weight = 0;
  cut_he_weight.store(0, std::memory_order_relaxed);
}

bool QuotientGraph::ActiveBlockSchedulingRound::popBlockPairFromQueue(BlockPair& blocks) {
  blocks.i = kInvalidPartition;
  blocks.j = kInvalidPartition;
  const size_t current_size = _unscheduled_blocks.unsafe_size();
  size_t current_idx = 0;
  while ( current_idx < current_size && _unscheduled_blocks.try_pop(blocks) ) {
    _quotient_graph[blocks.i][blocks.j].markAsNotInQueue();
    const bool not_too_many_concurrent_searches =
       _num_active_searches_on_blocks[blocks.i] < _context.refinement.advanced.max_concurrency_per_block &&
       _num_active_searches_on_blocks[blocks.j] < _context.refinement.advanced.max_concurrency_per_block;
    if ( !not_too_many_concurrent_searches ) {
      const bool success = pushBlockPairIntoQueue(blocks);
      _remaining_blocks -= (1 - success);
      blocks.i = kInvalidPartition;
      blocks.j = kInvalidPartition;
    } else {
      break;
    }
    ++current_idx;
  }
  return blocks.i != kInvalidPartition && blocks.j != kInvalidPartition;
}


void QuotientGraph::ActiveBlockSchedulingRound::finalizeSearch(const BlockPair& blocks,
                                                               const HyperedgeWeight improvement,
                                                               bool& block_0_becomes_active,
                                                               bool& block_1_becomes_active) {
  _round_improvement += improvement;
  --_remaining_blocks;
  if ( improvement > 0 ) {
    _active_blocks_lock.lock();
    block_0_becomes_active = !_active_blocks[blocks.i];
    block_1_becomes_active = !_active_blocks[blocks.j];
    _active_blocks[blocks.i] = true;
    _active_blocks[blocks.j] = true;
    _active_blocks_lock.unlock();
  }
}

bool QuotientGraph::ActiveBlockSchedulingRound::pushBlockPairIntoQueue(const BlockPair& blocks) {
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  if ( qg_edge.markAsInQueue() ) {
    _unscheduled_blocks.push(blocks);
    ++_remaining_blocks;
    return true;
  } else {
    return false;
  }
}

void QuotientGraph::ActiveBlockScheduler::initialize(const vec<uint8_t>& active_blocks,
                                                     const bool is_input_hypergraph) {
  reset();
  _is_input_hypergraph = is_input_hypergraph;

  HyperedgeWeight best_total_improvement = 1;
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      best_total_improvement = std::max(best_total_improvement,
        _quotient_graph[i][j].total_improvement.load(std::memory_order_relaxed));
    }
  }

  auto should_accept_stable_block_pair = [&](const PartitionID i, const PartitionID j) {
    const double relative_improvement_to_best =
      _quotient_graph[i][j].total_improvement.load() /
      static_cast<double>(best_total_improvement);
    return relative_improvement_to_best >
      _context.refinement.advanced.stable_block_relative_improvement_threshold;
  };

  vec<BlockPair> active_block_pairs;
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      if ( isActiveBlockPair(i, j, 0) &&
          ( active_blocks[i] || active_blocks[j] || should_accept_stable_block_pair(i, j) ) ) {
        active_block_pairs.push_back( BlockPair { i, j } );
      }
    }
  }

  if ( active_block_pairs.size() > 0 ) {
    std::sort(active_block_pairs.begin(), active_block_pairs.end(),
      [&](const BlockPair& lhs, const BlockPair& rhs) {
        return _quotient_graph[lhs.i][lhs.j].total_improvement >
          _quotient_graph[rhs.i][rhs.j].total_improvement ||
          ( _quotient_graph[lhs.i][lhs.j].total_improvement ==
            _quotient_graph[rhs.i][rhs.j].total_improvement &&
            _quotient_graph[lhs.i][lhs.j].cut_he_weight >
            _quotient_graph[rhs.i][rhs.j].cut_he_weight );
      });
    _rounds.emplace_back(_context, _quotient_graph, _num_active_searches_on_blocks);
    ++_num_rounds;
    for ( const BlockPair& blocks : active_block_pairs ) {
      DBG << "Schedule blocks (" << blocks.i << "," << blocks.j << ") in round 1 ("
          << "Total Improvement =" << _quotient_graph[blocks.i][blocks.j].total_improvement << ","
          << "Cut Weight =" << _quotient_graph[blocks.i][blocks.j].cut_he_weight << ")";
      _rounds.back().pushBlockPairIntoQueue(blocks);
    }
  }
}

bool QuotientGraph::ActiveBlockScheduler::popBlockPairFromQueue(BlockPair& blocks, size_t& round) {
  bool success = false;
  round = _first_active_round;
  while ( !_terminate && round < _num_rounds ) {
    success = _rounds[round].popBlockPairFromQueue(blocks);
    if ( success ) {
      break;
    }
    ++round;
  }

  if ( success && round == _num_rounds - 1 ) {
    _round_lock.lock();
    if ( round == _num_rounds - 1 ) {
      // There must always be a next round available such that we can
      // reschedule block pairs that become active.
      _rounds.emplace_back(_context, _quotient_graph, _num_active_searches_on_blocks);
      ++_num_rounds;
    }
    _round_lock.unlock();
  }

  return success;
}

void QuotientGraph::ActiveBlockScheduler::finalizeSearch(const BlockPair& blocks,
                                                         const size_t round,
                                                         const HyperedgeWeight improvement) {
  ASSERT(round < _rounds.size());
  --_num_active_searches_on_blocks[blocks.i];
  --_num_active_searches_on_blocks[blocks.j];
  bool block_0_becomes_active = false;
  bool block_1_becomes_active = false;
  _rounds[round].finalizeSearch(blocks, improvement,
    block_0_becomes_active, block_1_becomes_active);

  if ( block_0_becomes_active ) {
    // If blocks.i becomes active, we push all adjacent blocks into the queue of the next round
    ASSERT(round + 1 < _rounds.size());
    for ( PartitionID j = blocks.i + 1; j < _context.partition.k; ++j ) {
      if ( isActiveBlockPair(blocks.i, j, round + 1) ) {
        DBG << "Schedule blocks (" << blocks.i << "," << j << ") in round" << (round + 2) << " ("
            << "Total Improvement =" << _quotient_graph[blocks.i][j].total_improvement << ","
            << "Cut Weight =" << _quotient_graph[blocks.i][j].cut_he_weight << ")";
        _rounds[round + 1].pushBlockPairIntoQueue(BlockPair { blocks.i, j });
      }
    }
  }

  if ( block_1_becomes_active ) {
    // If blocks.j becomes active, we push all adjacent blocks into the queue of the next round
    ASSERT(round + 1 < _rounds.size());
    for ( PartitionID j = blocks.j + 1; j < _context.partition.k; ++j ) {
      if ( isActiveBlockPair(blocks.j, j, round + 1) ) {
        DBG << "Schedule blocks (" << blocks.j << "," << j << ") in round" << (round + 2) << " ("
            << "Total Improvement =" << _quotient_graph[blocks.j][j].total_improvement << ","
            << "Cut Weight =" << _quotient_graph[blocks.j][j].cut_he_weight << ")";
        _rounds[round + 1].pushBlockPairIntoQueue(BlockPair { blocks.j, j });
      }
    }
  }

  if ( round == _first_active_round && _rounds[round].numRemainingBlocks() == 0 ) {
    _round_lock.lock();
    // We consider a round as finished, if the previous round is also finished and there
    // are no remaining blocks in the queue of that round.
    while ( _first_active_round < _rounds.size() &&
            _rounds[_first_active_round].numRemainingBlocks() == 0 ) {
      DBG << GREEN << "Round" << (_first_active_round + 1) << "terminates with improvement"
          << _rounds[_first_active_round].roundImprovement() << "("
          << "Minimum Required Improvement =" << _min_improvement_per_round << ")" << END;
      // We require that minimum improvement per round must be greater than a threshold,
      // otherwise we terminate early
      _terminate = _rounds[_first_active_round].roundImprovement() < _min_improvement_per_round;
      ++_first_active_round;
    }
    _round_lock.unlock();
  }
}


bool QuotientGraph::ActiveBlockScheduler::isActiveBlockPair(const PartitionID i,
                                                            const PartitionID j,
                                                            const size_t round) const {
  const bool skip_small_cuts = !_is_input_hypergraph &&
    _context.refinement.advanced.skip_small_cuts;
  const bool contains_enough_cut_hes =
    (skip_small_cuts && _quotient_graph[i][j].cut_he_weight > 10) ||
    (!skip_small_cuts && _quotient_graph[i][j].cut_he_weight > 0);
  const bool is_promising_blocks_pair =
    !_context.refinement.advanced.skip_unpromising_blocks ||
      ( round == 0 || _quotient_graph[i][j].num_improvements_found > 0 );
  return contains_enough_cut_hes && is_promising_blocks_pair;
}

SearchID QuotientGraph::requestNewSearch(AdvancedRefinerAdapter& refiner) {
  ASSERT(_phg);
  SearchID search_id = INVALID_SEARCH_ID;
  BlockPair blocks { kInvalidPartition, kInvalidPartition };
  size_t round = 0;
  bool success = _active_block_scheduler.popBlockPairFromQueue(blocks, round);
  _register_search_lock.lock();

  const SearchID tmp_search_id = _searches.size();
  if ( success && _quotient_graph[blocks.i][blocks.j].acquire(tmp_search_id) ) {
    ++_num_active_searches;
    // Create new search
    search_id = tmp_search_id;
    _searches.emplace_back(blocks, round);
    _active_block_scheduler.startSearch(blocks);
    _register_search_lock.unlock();

    // Associate refiner with search id
    success = refiner.registerNewSearch(search_id, *_phg);
    ASSERT(success); unused(success);
  } else {
    _register_search_lock.unlock();
    if ( success ) {
      _active_block_scheduler.finalizeSearch(blocks, round, 0);
    }
  }
  return search_id;
}

BlockPairCutHyperedges QuotientGraph::requestCutHyperedges(const SearchID search_id,
                                                           const size_t max_num_edges) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  BlockPairCutHyperedges block_pair_cut_hes;
  if ( !_searches[search_id].is_finalized ) {
    Search& search = _searches[search_id];
    const BlockPair blocks = search.blocks;
    QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
    block_pair_cut_hes.blocks = blocks;

    size_t num_edges = 0;
    while ( num_edges < max_num_edges && qg_edge.cut_he_weight > 0 ) {
      // Note, we only consider hyperedges that contains pins of both blocks.
      // There might be some edges which were initially cut, but were removed due
      // to vertex moves. Thus, we remove them lazily here.
      const HyperedgeID he = qg_edge.pop_hyperedge();
      qg_edge.cut_he_weight -= _phg->edgeWeight(he);
      if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
          _phg->pinCountInPart(he, blocks.j) > 0 ) {
        block_pair_cut_hes.cut_hes.push_back(he);
        _searches[search_id].used_cut_hes.push_back(he);
        ++num_edges;
      }
    }
  }
  return block_pair_cut_hes;
}

size_t QuotientGraph::acquireUsedCutHyperedges(const SearchID& search_id, vec<bool>& used_hes) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  ASSERT(!_searches[search_id].is_finalized);
  size_t additional_cut_hes = 0;
  const BlockPair& blocks = _searches[search_id].blocks;
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  const size_t used_hes_size = _searches[search_id].used_cut_hes.size();
  const size_t start_idx = qg_edge.first_valid_entry;
  const size_t end_idx = qg_edge.cut_hes.size();
  for ( size_t i = start_idx; i < end_idx; ++i ) {
    const HyperedgeID he = qg_edge.cut_hes[i];
    if ( used_hes[he] ) {
      used_hes[he] = false;
      qg_edge.cut_he_weight -= _phg->edgeWeight(he);
      _searches[search_id].used_cut_hes.push_back(he);
      // Note that only one thread constructs a problem at any particular point in time.
      // This function is called before construction terminates. Therefore, this
      // operation is thread safe.
      std::swap(qg_edge.cut_hes[qg_edge.first_valid_entry++], qg_edge.cut_hes[i]);
      ++additional_cut_hes;
    }
  }
  // Flag used cut hes again
  for ( size_t i = used_hes_size; i < _searches[search_id].used_cut_hes.size(); ++i ) {
    used_hes[_searches[search_id].used_cut_hes[i]] = true;
  }
  return additional_cut_hes;
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
  const BlockPair& blocks = _searches[search_id].blocks;
  _quotient_graph[blocks.i][blocks.j].release(search_id);
}

void QuotientGraph::finalizeSearch(const SearchID search_id,
                                   const HyperedgeWeight total_improvement) {
  ASSERT(_phg);
  ASSERT(search_id < _searches.size());
  ASSERT(_searches[search_id].is_finalized);

  const BlockPair& blocks = _searches[search_id].blocks;
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  if ( total_improvement > 0 ) {
    // If the search improves the quality of the partition, we reinsert
    // all hyperedges that were used by the search and are still cut.
    ++qg_edge.num_improvements_found;
    qg_edge.total_improvement += total_improvement;
    for ( const HyperedgeID& he : _searches[search_id].used_cut_hes ) {
      if ( _phg->pinCountInPart(he, blocks.i) > 0 &&
          _phg->pinCountInPart(he, blocks.j) > 0 ) {
        qg_edge.add_hyperedge(he, _phg->edgeWeight(he));
      }
    }
  }
  // In case the block pair becomes active,
  // we reinsert it into the queue
  _active_block_scheduler.finalizeSearch(
    blocks, _searches[search_id].round, total_improvement);
  --_num_active_searches;
}

void QuotientGraph::initialize(const PartitionedHypergraph& phg) {
  _phg = &phg;

  // Reset internal members
  resetQuotientGraphEdges();
  _num_active_searches.store(0, std::memory_order_relaxed);
  _searches.clear();

  // Find all cut hyperedges between the blocks
  tbb::enumerable_thread_specific<HyperedgeID> local_num_hes(0);
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
    ++local_num_hes.local();
    const HyperedgeWeight edge_weight = phg.edgeWeight(he);
    for ( const PartitionID i : phg.connectivitySet(he) ) {
      for ( const PartitionID j : phg.connectivitySet(he) ) {
        if ( i < j ) {
          _quotient_graph[i][j].add_hyperedge(he, edge_weight);
        }
      }
    }
  });
  const HyperedgeID tmp_num_edges = local_num_hes.combine(std::plus<HyperedgeID>());
  const bool is_same_hypergraph = _current_num_edges == tmp_num_edges;
  _current_num_edges = tmp_num_edges;

  vec<uint8_t> active_blocks(_context.partition.k, false);
  if ( is_same_hypergraph &&
       _context.refinement.advanced.skip_stable_blocks ) {
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID prev_id = _partition_snapshot[hn];
      const PartitionID cur_id = phg.partID(hn);
      if ( prev_id != kInvalidPartition && prev_id != cur_id ) {
        active_blocks[prev_id] = true;
        active_blocks[cur_id] = true;
      }
    });
  } else {
    active_blocks.assign(_context.partition.k, true);
  }

  // Initalize block scheduler queue
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].initial_cut_he_weight = _quotient_graph[i][j].cut_he_weight;
      _quotient_graph[i][j].initial_num_cut_hes =  _quotient_graph[i][j].cut_hes.size();
    }
  }
  _active_block_scheduler.initialize(active_blocks, isInputHypergraph());

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

size_t QuotientGraph::maximumRequiredRefiners() const {
  const size_t current_active_block_pairs =
    _active_block_scheduler.numRemainingBlocks() + _num_active_searches + 1;
  return std::min(current_active_block_pairs, _context.shared_memory.num_threads);
}

void QuotientGraph::resetQuotientGraphEdges() {
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].reset();
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