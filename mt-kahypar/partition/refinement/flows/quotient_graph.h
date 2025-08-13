/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

static constexpr SearchID INVALID_SEARCH_ID = std::numeric_limits<SearchID>::max();

struct BlockPair {
  PartitionID i = kInvalidPartition;
  PartitionID j = kInvalidPartition;
};


// ! Represents an edge of the quotient graph
struct QuotientGraphEdge {
  QuotientGraphEdge() :
    blocks(),
    ownership(INVALID_SEARCH_ID),
    is_in_queue(false),
    cut_hes(),
    num_cut_hes(0),
    cut_he_weight(0),
    num_improvements_found(0),
    total_improvement(0) { }

  // ! Adds a cut hyperedge to this quotient graph edge
  void add_hyperedge(const HyperedgeID he, const HyperedgeWeight weight);

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

  bool isInQueue() const {
    return is_in_queue.load(std::memory_order_relaxed);
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

  // ! Block pair this quotient graph edge represents
  BlockPair blocks;
  // ! Atomic that contains the search currently constructing
  // ! a problem on this block pair
  CAtomic<SearchID> ownership;
  // ! True, if block is contained in block scheduler queue
  CAtomic<bool> is_in_queue;
  // ! Cut hyperedges of block pair
  tbb::concurrent_vector<HyperedgeID> cut_hes;
  // ! Number of cut hyperedges
  CAtomic<size_t> num_cut_hes;
  // ! Current weight of all cut hyperedges
  CAtomic<HyperedgeWeight> cut_he_weight;
  // ! Number of improvements found on this block pair
  CAtomic<size_t> num_improvements_found;
  // ! Total improvement found on this block pair
  CAtomic<HyperedgeWeight> total_improvement;
};


class QuotientGraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

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
  explicit QuotientGraph(const HyperedgeID num_hyperedges, const Context& context);

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
  template<typename TypeTraits>
  SearchID requestNewSearch(typename TypeTraits::PartitionedHypergraph& phg,
                            FlowRefinerAdapter<TypeTraits>& refiner,
                            BlockPair blocks,
                            size_t round) {
    SearchID search_id = INVALID_SEARCH_ID;

    _register_search_lock.lock();
    const SearchID tmp_search_id = _searches.size();
    if ( _quotient_graph[blocks.i][blocks.j].acquire(tmp_search_id) ) {
      ++_num_active_searches;
      // Create new search
      search_id = tmp_search_id;
      _searches.emplace_back(blocks, round);
      _register_search_lock.unlock();

      // Associate refiner with search id
      bool success = refiner.registerNewSearch(search_id, phg);
      ASSERT(success); unused(success);
    } else {
      _register_search_lock.unlock();
    }
    return search_id;
  }

  // ! Returns the block pair on which the corresponding search operates on
  BlockPair getBlockPair(const SearchID search_id) const;

  size_t getRound(const SearchID search_id) const;

  vec<vec<QuotientGraphEdge>>& getGraph();

  size_t numActiveSearches() const;

  template<typename PartitionedHypergraph, typename F>
  void doForAllCutHyperedgesOfSearch(const PartitionedHypergraph& phg, const SearchID search_id, const F& f) {
    const BlockPair& blocks = _searches[search_id].blocks;
    const size_t num_cut_hes = _quotient_graph[blocks.i][blocks.j].num_cut_hes.load();
    std::shuffle(_quotient_graph[blocks.i][blocks.j].cut_hes.begin(),
                 _quotient_graph[blocks.i][blocks.j].cut_hes.begin() + num_cut_hes,
                 utils::Randomize::instance().getGenerator());
    for ( size_t i = 0; i < num_cut_hes; ++i ) {
      const HyperedgeID he = _quotient_graph[blocks.i][blocks.j].cut_hes[i];
      if ( phg.pinCountInPart(he, blocks.i) > 0 && phg.pinCountInPart(he, blocks.j) > 0 ) {
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
  template<typename PartitionedHypergraph>
  void addNewCutHyperedges(const PartitionedHypergraph& phg, const vec<std::pair<HyperedgeID, PartitionID>>& new_cut_hes);

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
  template<typename PartitionedHypergraph>
  void initialize(const PartitionedHypergraph& phg);

  // ! Only for testing
  HyperedgeWeight getCutHyperedgeWeightOfBlockPair(const PartitionID i, const PartitionID j) const;

  void changeNumberOfBlocks(const PartitionID new_k);

  bool isInputHypergraph() const;

 private:
  void resetQuotientGraphEdges();


  const Context& _context;
  const HypernodeID _initial_num_edges;
  HypernodeID _current_num_edges;

  // ! Each edge contains stats and the cut hyperedges
  // ! of the block pair which its represents.
  vec<vec<QuotientGraphEdge>> _quotient_graph;

  SpinLock _register_search_lock;

  // ! Number of active searches
  CAtomic<size_t> _num_active_searches;
  // ! Information about searches that are currently running
  tbb::concurrent_vector<Search> _searches;
};

}  // namespace kahypar
