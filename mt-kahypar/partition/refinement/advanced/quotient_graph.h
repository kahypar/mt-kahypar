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

class QuotientGraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  // ! Represents an edge of the quotient graph
  struct QuotientGraphEdge {
    QuotientGraphEdge() :
      blocks(),
      ownership(INVALID_SEARCH_ID),
      cut_hes(),
      initial_num_cut_hes(0),
      cut_he_weight(0),
      num_improvements_found(0),
      total_improvement(0) { }

    // ! Adds a cut hyperedge to this quotient graph edge
    void add_hyperedge(const HyperedgeID he,
                       const HyperedgeWeight weight);

    void reset();

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

    // ! Block pair this quotient graph edge represents
    BlockPair blocks;
    // ! Atomic that contains the search currently constructing
    // ! a problem on this block pair
    CAtomic<SearchID> ownership;
    // ! Cut hyperedges of block pair
    tbb::concurrent_vector<HyperedgeID> cut_hes;
    // ! Initial number of cut hyperedges
    size_t initial_num_cut_hes;
    // ! Current weight of all cut hyperedges
    CAtomic<HyperedgeWeight> cut_he_weight;
    // ! Number of improvements found on this block pair
    CAtomic<size_t> num_improvements_found;
    // ! Total improvement found on this block pair
    CAtomic<HyperedgeWeight> total_improvement;
  };

  /**
   * Maintains the block pair of a round of the active block scheduling strategy
   */
  class ActiveBlockSchedulingRound {

   public:
    explicit ActiveBlockSchedulingRound(const Context& context) :
      _context(context),
      _unscheduled_blocks(),
      _round_improvement(0),
      _active_blocks_lock(),
      _active_blocks(context.partition.k, false),
      _remaining_blocks(0) { }

    // ! Pops a block pair from the queue.
    // ! Returns true, if a block pair was successfully popped from the queue.
    // ! The corresponding block pair will be stored in blocks.
    bool popBlockPairFromQueue(BlockPair& blocks);

    // ! Pushes a block pair into the queue.
    // ! Note, that a block pair is only allowed to be contained in one queue
    // ! (there are multiple active rounds).
    void pushBlockPairIntoQueue(const BlockPair& blocks);

    // ! Signals that the search on the corresponding block pair terminated.
    void finalizeSearch(const BlockPair& blocks,
                        const HyperedgeWeight improvement,
                        bool& block_0_becomes_active,
                        bool& block_1_becomes_active);

    HyperedgeWeight roundImprovement() const {
      return _round_improvement.load(std::memory_order_relaxed);
    }

    size_t numRemainingBlocks() const {
      return _remaining_blocks;
    }

   const Context& _context;
    // ! Queue that contains all unscheduled block pairs of the current round
    tbb::concurrent_queue<BlockPair> _unscheduled_blocks;
    // ! Current improvement made in this round
    CAtomic<HyperedgeWeight> _round_improvement;
    // Active blocks for next round
    SpinLock _active_blocks_lock;
    vec<uint8_t> _active_blocks;
    // Remaining active block pairs in the current round.
    CAtomic<size_t> _remaining_blocks;
  };

  /**
   * Implements the active block scheduling strategy.
   * The active block scheduling strategy proceeds in rounds. In each round,
   * all active edges of the quotient graph are scheduled for refinement.
   * A edge is called active, if at least of the blocks is active and a block
   * is called active if a refinement involving that block in the previous round
   * leads to an improvement. In the sequential active block scheduling strategy
   * the rounds acts as synchronization barriers. However, to achieve better scalibility
   * we immediatly schedule an edge in the next round once we find an improvement.
   * Thus, there can be multiple active searches that process block pairs from different
   * rounds. However, block pairs from earlier rounds have an higher priority to be
   * scheduled.
   */
  class ActiveBlockScheduler {

   public:
    explicit ActiveBlockScheduler(const Context& context,
                                  vec<vec<QuotientGraphEdge>>& quotient_graph) :
      _context(context),
      _quotient_graph(quotient_graph),
      _num_rounds(0),
      _rounds(),
      _min_improvement_per_round(0),
      _terminate(false),
      _round_lock(),
      _first_active_round(0),
      _is_input_hypergraph(false) { }

    // ! Initialize the first round of the active block scheduling strategy
    void initialize(const vec<uint8_t>& active_blocks,
                    const bool is_input_hypergraph);

    // ! Pops a block pair from the queue.
    // ! Returns true, if a block pair was successfully popped from the queue.
    // ! The corresponding block pair and the round to which this blocks corresponds
    // ! to are stored in blocks and round.
    bool popBlockPairFromQueue(BlockPair& blocks, size_t& round);

    // ! Signals that the search on the corresponding block pair terminated.
    // ! If one the two blocks become active, we immediatly schedule all edges
    // ! adjacent in the quotient graph in the next round of active block scheduling
    void finalizeSearch(const BlockPair& blocks,
                        const size_t round,
                        const HyperedgeWeight improvement);

    size_t numRemainingBlocks() const {
      size_t num_remaining_blocks = 0;
      for ( size_t i = _first_active_round; i < _rounds.size(); ++i ) {
        num_remaining_blocks += _rounds[i].numRemainingBlocks();
      }
      return num_remaining_blocks;
    }

    void setObjective(const HyperedgeWeight objective) {
      _min_improvement_per_round =
        _context.refinement.advanced.min_relative_improvement_per_round * objective;
    }

   private:

    void reset() {
      _num_rounds.store(0);
      _rounds.clear();
      _first_active_round = 0;
      _terminate = false;
    }

    bool isActiveBlockPair(const PartitionID i,
                           const PartitionID j) const;

    const Context& _context;
    // ! Quotient graph
    vec<vec<QuotientGraphEdge>>& _quotient_graph;
    // Contains all active block scheduling rounds
    CAtomic<size_t> _num_rounds;
    tbb::concurrent_vector<ActiveBlockSchedulingRound> _rounds;
    // ! Minimum improvement per round to continue with next round
    HyperedgeWeight _min_improvement_per_round;
    // ! If true, then search is immediatly terminated
    bool _terminate;
    // ! First Active Round
    SpinLock _round_lock;
    size_t _first_active_round;
    // ! Indicate if the current hypergraph represents the input hypergraph
    bool _is_input_hypergraph;
  };

  // Contains information required by a local search
  struct Search {
    explicit Search(const BlockPair& blocks, const size_t round) :
      blocks(blocks),
      round(round),
      is_finalized(false) { }

    // ! Block pair on which this search operates on
    BlockPair blocks;
    // ! Round of active block scheduling
    size_t round;
    // ! Flag indicating if construction of the corresponding search
    // ! is finalized
    bool is_finalized;
  };

public:
  static constexpr SearchID INVALID_SEARCH_ID = std::numeric_limits<SearchID>::max();

  explicit QuotientGraph(const Hypergraph& hg,
                         const Context& context) :
    _phg(nullptr),
    _context(context),
    _initial_num_edges(hg.initialNumEdges()),
    _current_num_edges(kInvalidHyperedge),
    _quotient_graph(context.partition.k,
      vec<QuotientGraphEdge>(context.partition.k)),
    _register_search_lock(),
    _active_block_scheduler(context, _quotient_graph),
    _num_active_searches(0),
    _searches(),
    _partition_snapshot() {
    for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
      for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
        _quotient_graph[i][j].blocks.i = i;
        _quotient_graph[i][j].blocks.j = j;
      }
    }

    if ( doPartitionSnapshot() ) {
      _partition_snapshot.assign(hg.initialNumNodes(), kInvalidPartition);
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

  template<typename F>
  void doForAllCutHyperedgesOfSearch(const SearchID search_id, const F& f) {
    const BlockPair& blocks = _searches[search_id].blocks;
    std::random_shuffle(_quotient_graph[blocks.i][blocks.j].cut_hes.begin(),
                        _quotient_graph[blocks.i][blocks.j].cut_hes.end());
    for ( const HyperedgeID& he : _quotient_graph[blocks.i][blocks.j].cut_hes ) {
      if ( _phg->pinCountInPart(he, blocks.i) > 0 && _phg->pinCountInPart(he, blocks.j) > 0 ) {
        f(he);
      }
    }
  }


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
                      const HyperedgeWeight total_improvement);

  // ! Initializes the quotient graph. This includes to find
  // ! all cut hyperedges between all block pairs
  void initialize(const PartitionedHypergraph& phg);

  void setObjective(const HyperedgeWeight objective) {
    _active_block_scheduler.setObjective(objective);
  }

  size_t maximumRequiredRefiners() const;

  // ! Only for testing
  HyperedgeWeight getCutHyperedgeWeightOfBlockPair(const PartitionID i, const PartitionID j) const {
    ASSERT(i < j);
    ASSERT(0 <= i && i < _context.partition.k);
    ASSERT(0 <= j && j < _context.partition.k);
    return _quotient_graph[i][j].cut_he_weight;
  }

  void storePartition(const PartitionedHypergraph& phg) {
    if ( doPartitionSnapshot() ) {
      phg.doParallelForAllNodes([&](const HypernodeID hn) {
        _partition_snapshot[hn] = phg.partID(hn);
      });
    }
  }

 private:

  void resetQuotientGraphEdges();

  bool isInputHypergraph() const {
    return _current_num_edges == _initial_num_edges;
  }

  bool doPartitionSnapshot() const {
    return ( ( _context.partition.paradigm == Paradigm::nlevel &&
             _context.refinement.global_fm.refine_until_no_improvement ) ||
           ( _context.partition.paradigm == Paradigm::multilevel &&
             _context.refinement.refine_until_no_improvement ) ) &&
           _context.refinement.advanced.skip_stable_blocks;
  }

  const PartitionedHypergraph* _phg;
  const Context& _context;
  const HypernodeID _initial_num_edges;
  HypernodeID _current_num_edges;

  // ! Each edge contains stats and the cut hyperedges
  // ! of the block pair which its represents.
  vec<vec<QuotientGraphEdge>> _quotient_graph;

  SpinLock _register_search_lock;
  // ! Queue that contains all block pairs.
  ActiveBlockScheduler _active_block_scheduler;

  // ! Number of active searches
  CAtomic<size_t> _num_active_searches;
  // ! Information about searches that are currently running
  tbb::concurrent_vector<Search> _searches;

  // ! Snapshot of the partition to detect changes, if called
  // ! multiple times (refine until no improvement)
  vec<PartitionID> _partition_snapshot;
};

}  // namespace kahypar
