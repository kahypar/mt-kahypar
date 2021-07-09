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

class QuotientGraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  // ! Represents an edge of the quotient graph
  struct QuotientGraphEdge {
    QuotientGraphEdge() :
      skip_small_cuts(false),
      blocks(),
      ownership(INVALID_SEARCH_ID),
      is_in_queue(false),
      cut_hes(),
      first_valid_entry(0),
      initial_num_cut_hes(0),
      initial_cut_he_weight(0),
      cut_he_weight(0),
      num_improvements_found(0) { }

    // ! Adds a cut hyperedge to this quotient graph edge
    void add_hyperedge(const HyperedgeID he,
                       const HyperedgeWeight weight);

    // ! Pops a cut hyperedge from this quotient graph edge
    HyperedgeID pop_hyperedge();

    void reset();

    // ! A quotient graph edge is considered as active, if there
    // ! are cut hyperedges left and if there weight is greater than
    // ! some threshold. Only active quotient graph edges are available
    // ! for refinement.
    bool isActive() const {
      return ( !skip_small_cuts && cut_he_weight > 0 ) ||
             ( skip_small_cuts && cut_he_weight > 10 );
    }

    // ! Returns true, if quotient graph edge is acquired by a search
    bool isAcquired() const {
      return ownership.load() != INVALID_SEARCH_ID;
    }

    // ! Tries to acquire quotient graph edge with corresponding search id
    bool acquire(const SearchID search_id) {
      SearchID expected = INVALID_SEARCH_ID;
      SearchID desired = search_id;
      return ownership.compare_exchange_strong(expected, desired);
    }

    // ! Releases quotient graph edge
    void release(const SearchID search_id) {
      unused(search_id);
      ASSERT(ownership.load() == search_id);
      ownership.store(INVALID_SEARCH_ID);
    }

    // ! Marks quotient graph edge as in queue. Queued edges are scheduled
    // ! for refinement.
    bool markAsInQueue() {
      bool expected = false;
      bool desired = true;
      return is_in_queue.compare_exchange_strong(expected, desired);
    }

    // ! Marks quotient graph edge as nnot in queue
    bool markAsNotInQueue() {
      bool expected = true;
      bool desired = false;
      return is_in_queue.compare_exchange_strong(expected, desired);
    }

    // ! Returns false, if all cut hyperedges initial contained are used
    // ! at least once by a refinement algorithm
    bool isFirstRound() const {
      return first_valid_entry < initial_num_cut_hes;
    }

    // ! If true, then block pairs with cut less than 10 are skipped
    bool skip_small_cuts;
    // ! Block pair this quotient graph edge represents
    BlockPair blocks;
    // ! Atomic that contains the search currently constructing
    // ! a problem on this block pair
    CAtomic<SearchID> ownership;
    // ! True, if block is contained in block scheduler queue
    CAtomic<bool> is_in_queue;
    // ! Cut hyperedges of block pair
    tbb::concurrent_vector<HyperedgeID> cut_hes;
    // ! Position of the first valid cut hyperedge in cut_hes
    size_t first_valid_entry;
    // ! Initial number of cut hyperedges
    size_t initial_num_cut_hes;
    // ! Initial weight of all cut hyperedges
    HyperedgeWeight initial_cut_he_weight;
    // ! Current weight of all cut hyperedges
    CAtomic<HyperedgeWeight> cut_he_weight;
    // ! Number of improvements found on this block pair
    CAtomic<size_t> num_improvements_found;
  };

  // Contains information required by a local search
  struct Search {
    explicit Search(const BlockPair& blocks) :
      blocks(blocks),
      used_cut_hes(),
      is_finalized(false) { }

    // ! Block pair on which this search operates on
    BlockPair blocks;
    // ! Used cut hyperedges
    vec<HyperedgeID> used_cut_hes;
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

  // ! Returns the block pair on which the corresponding search operates on
  BlockPair getBlockPair(const SearchID search_id) const {
    ASSERT(search_id < _searches.size());
    return _searches[search_id].blocks;
  }

  // ! Number of block pairs used by the corresponding search
  size_t numBlockPairs(const SearchID) const {
    return 1;
  }

  /**
   * Requests cut hyperedges that contains the blocks
   * associated with the corresponding search.
   */
  BlockPairCutHyperedges requestCutHyperedges(const SearchID search_id,
                                              const size_t max_num_edges);

  /**
   * During problem construction we might acquire additional cut hyperedges
   * not requested by the search. This function associates those hyperedges
   * with the search and flags them as used.
   */
  size_t acquireUsedCutHyperedges(const SearchID& search_id, vec<bool>& used_hes);


  /**
   * Notifies the quotient graph that hyperedge he contains
   * a new block, which was previously not contained. The thread
   * that increases the pin count of hyperedge he in the corresponding
   * block to 1 is responsible to call this function.
   */
  void addNewCutHyperedge(const HyperedgeID he,
                          const PartitionID block);

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

  size_t maximumRequiredRefiners() const;

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
