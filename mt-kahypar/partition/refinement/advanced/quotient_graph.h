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

#pragma once

#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_vector.h"
#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/advanced/refiner_adapter.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {

struct BlockPair {
  PartitionID i = kInvalidPartition;
  PartitionID j = kInvalidPartition;
};

struct BlockPairCutHyperedges {
  BlockPairCutHyperedges() :
    blocks(),
    cut_hes() { }

  BlockPair blocks;
  vec<HyperedgeID> cut_hes;
};

struct BlockPairStats {
  BlockPairStats() :
    initial_cut_he_weight(0),
    num_acitve_searches_on_block_pair(0),
    max_active_searches_on_one_block(0) { }

  BlockPairStats(bool /* sentinel */) :
    initial_cut_he_weight(std::numeric_limits<HyperedgeWeight>::max()),
    num_acitve_searches_on_block_pair(std::numeric_limits<size_t>::min()),
    max_active_searches_on_one_block(std::numeric_limits<size_t>::min()) { }

  HyperedgeWeight initial_cut_he_weight;
  CAtomic<size_t> num_acitve_searches_on_block_pair;
  size_t max_active_searches_on_one_block;
};

struct BlockPairStatsComparator {
  bool operator()(const BlockPairStats& lhs, const BlockPairStats& rhs) const {
    return lhs.max_active_searches_on_one_block > rhs.max_active_searches_on_one_block ||
      ( lhs.max_active_searches_on_one_block == rhs.max_active_searches_on_one_block &&
        lhs.num_acitve_searches_on_block_pair > rhs.num_acitve_searches_on_block_pair ) ||
      ( lhs.max_active_searches_on_one_block == rhs.max_active_searches_on_one_block &&
        lhs.num_acitve_searches_on_block_pair == rhs.num_acitve_searches_on_block_pair &&
        lhs.initial_cut_he_weight < rhs.initial_cut_he_weight );
  }
};

bool operator==(const BlockPairStats& lhs, const BlockPairStats& rhs);

std::ostream& operator<<(std::ostream& out, const BlockPairStats& stats);

class QuotientGraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  struct CutHyperedge {
    HyperedgeID he;
    CAtomic<bool> acquired;

    bool acquire() {
      bool expected = false;
      bool desired = true;
      return acquired.compare_exchange_strong(expected, desired);
    }
  };

  // ! Represents an edge of the quotient graph
  struct QuotientGraphEdge {
    QuotientGraphEdge() :
      skip_small_cuts(false),
      blocks(),
      ownership(INVALID_SEARCH_ID),
      is_in_queue(false),
      first_valid_entry(0),
      cut_hes(),
      stats(),
      cut_he_weight(0),
      round(0),
      num_improvements(0) { }

    void add_hyperedge(const HyperedgeID he,
                       const HyperedgeWeight weight);

    CutHyperedge& pop_hyperedge();

    void reset();

    bool isActive() const {
      return ( !skip_small_cuts && cut_he_weight > 0 ) ||
             ( skip_small_cuts && cut_he_weight > 10 );
    }

    bool isAcquired() const {
      return ownership.load() != INVALID_SEARCH_ID;
    }

    bool acquire(const SearchID search_id) {
      SearchID expected = INVALID_SEARCH_ID;
      SearchID desired = search_id;
      return ownership.compare_exchange_strong(expected, desired);
    }

    void release(const SearchID search_id) {
      unused(search_id);
      ASSERT(ownership.load() == search_id);
      ownership.store(INVALID_SEARCH_ID);
    }

    bool markAsInQueue() {
      bool expected = false;
      bool desired = true;
      return is_in_queue.compare_exchange_strong(expected, desired);
    }

    bool markAsNotInQueue() {
      bool expected = true;
      bool desired = false;
      return is_in_queue.compare_exchange_strong(expected, desired);
    }

    bool skip_small_cuts;
    BlockPair blocks;
    CAtomic<SearchID> ownership;
    CAtomic<bool> is_in_queue;
    size_t first_valid_entry;
    tbb::concurrent_vector<CutHyperedge> cut_hes;
    BlockPairStats stats;
    CAtomic<HyperedgeWeight> cut_he_weight;
    CAtomic<size_t> round;
    CAtomic<size_t> num_improvements;
  };

  // Contains information required by a local search
  struct Search {
    explicit Search() :
      block_pairs(),
      used_cut_hes(),
      is_finalized(false) { }

    void addBlockPair(const BlockPair& blocks) {
      block_pairs.emplace_back(blocks);
      used_cut_hes.emplace_back();
    }

    // ! Blocks on which this search operates on
    vec<BlockPair> block_pairs;
    // ! Used cut hyperedges
    vec<vec<HyperedgeID>> used_cut_hes;
    // ! Flag indicating if construction of the corresponding search
    // ! is finalized
    bool is_finalized;
  };

  struct BFSData {
    BFSData(const HypernodeID num_nodes,
            const HyperedgeID num_edges) :
      visited_hns(num_nodes, false),
      distance(num_edges, -1) {}

    void reset() {
      std::fill(visited_hns.begin(), visited_hns.end(), false);
      std::fill(distance.begin(), distance.end(), -1);
    }

    vec<bool> visited_hns;
    vec<int> distance;
  };

public:
  static constexpr SearchID INVALID_SEARCH_ID = std::numeric_limits<SearchID>::max();

  explicit QuotientGraph(const Hypergraph& hg,
                         const Context& context) :
    _phg(nullptr),
    _context(context),
    _initial_num_nodes(hg.initialNumNodes()),
    _quotient_graph(context.partition.k,
      vec<QuotientGraphEdge>(context.partition.k)),
    _register_search_lock(),
    _block_scheduler(),
    _num_active_searches(0),
    _searches(),
    _num_active_searches_on_blocks(context.partition.k, CAtomic<size_t>(0)),
    _local_bfs(hg.initialNumNodes(), hg.initialNumEdges()) {
    for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
      for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
        _quotient_graph[i][j].blocks.i = i;
        _quotient_graph[i][j].blocks.j = j;
      }
    }
  }

  QuotientGraph(const QuotientGraph&) = delete;
  QuotientGraph(QuotientGraph&&) = delete;

  QuotientGraph & operator= (const QuotientGraph &) = delete;
  QuotientGraph & operator= (QuotientGraph &&) = delete;

  /**
   * Returns a new search id which is associated with a certain number
   * of block pairs. The corresponding search can request hyperedges
   * with the search id that are cut between the corresponding blocks
   * associated with the search. If there are currently no block pairs
   * available then INVALID_SEARCH_ID is returned.
   */
  SearchID requestNewSearch(AdvancedRefinerAdapter& refiner);

  // ! Returns the block pairs on which the corresponding search operates on
  vec<BlockPair> getBlockPairs(const SearchID search_id) const {
    ASSERT(search_id < _searches.size());
    return _searches[search_id].block_pairs;
  }

  // ! Number of block pairs used by the corresponding search
  size_t numBlockPairs(const SearchID search_id) const {
    ASSERT(search_id < _searches.size());
    return _searches[search_id].block_pairs.size();
  }

  /**
   * Requests cut hyperedges that contains the blocks
   * associated with the corresponding search.
   */
  vec<BlockPairCutHyperedges> requestCutHyperedges(const SearchID search_id,
                                                   const size_t max_num_edges);

  /**
   * Notifies the quotient graph that hyperedge he contains
   * a new block, which was previously not contained. The thread
   * that increases the pin count of hyperedge he in the corresponding
   * block to 1 is responsible to call this function.
   */
  void addNewCutHyperedge(const HyperedgeID he,
                          const PartitionID block);

  size_t acquireUsedCutHyperedges(const SearchID& search_id, const vec<bool>& used_hes);

  /**
   * Notify the quotient graph that the construction of the corresponding
   * search is completed. The corresponding block pairs associated with the
   * search are made available again for other searches.
   */
  void finalizeConstruction(const SearchID search_id);

  /**
   * Notify the quotient graph that the corrseponding search terminated.
   * If the search improves the quality of the partition (success == true),
   * we reinsert all hyperedges that were used throughout the construction
   * and are still cut between the corresponding block.
   */
  void finalizeSearch(const SearchID search_id,
                      const bool success);

  // ! Initializes the quotient graph. This includes to find
  // ! all cut hyperedges between all block pairs
  void initialize(const PartitionedHypergraph& phg);

  bool terminate() const {
    return _block_scheduler.empty() && _num_active_searches == 0;
  }

  // ! Only for testing
  HyperedgeWeight getCutHyperedgeWeightOfBlockPair(const PartitionID i, const PartitionID j) const {
    ASSERT(i < j);
    ASSERT(0 <= i && i < _context.partition.k);
    ASSERT(0 <= j && j < _context.partition.k);
    return _quotient_graph[i][j].cut_he_weight;
  }

 private:

  void resetQuotientGraphEdges(const PartitionedHypergraph& phg);

  bool popBlockPairFromQueue(BlockPair& blocks);

  bool pushBlockPairIntoQueue(const BlockPair& blocks);

  /**
   * The idea is to sort the cut hyperedges of a block pair (i,j)
   * in increasing order of their distance to each other. Meaning that
   * we perform a BFS from a randomly selected start cut hyperedge and
   * expand along cut hyperedges that contains pins of both blocks.
   * The BFS distance determines the order of the cut hyperedges.
   */
  void sortCutHyperedges(const PartitionID i,
                         const PartitionID j,
                         BFSData& bfs_data);

  bool isInputHypergraph(const PartitionedHypergraph& phg) const {
    return phg.initialNumNodes() == _initial_num_nodes;
  }

  const PartitionedHypergraph* _phg;
  const Context& _context;
  const HypernodeID _initial_num_nodes;

  // ! Each edge contains stats and the cut hyperedges
  // ! of the block pair which its represents.
  vec<vec<QuotientGraphEdge>> _quotient_graph;

  SpinLock _register_search_lock;
  // ! Queue that contains all block pairs.
  tbb::concurrent_queue<BlockPair> _block_scheduler;

  // ! Number of active searches
  CAtomic<size_t> _num_active_searches;
  // ! Information about searches that are currently running
  tbb::concurrent_vector<Search> _searches;

  // ! Number of active searches on each block
  vec<CAtomic<size_t>> _num_active_searches_on_blocks;

  // ! BFS data required to sort cut hyperedges
  tbb::enumerable_thread_specific<BFSData> _local_bfs;
};

}  // namespace kahypar
