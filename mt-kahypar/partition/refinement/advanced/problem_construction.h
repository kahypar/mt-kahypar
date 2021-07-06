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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/advanced/refiner_adapter.h"
#include "mt-kahypar/partition/refinement/advanced/quotient_graph.h"
#include "mt-kahypar/partition/refinement/advanced/problem_stats.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {

class ProblemConstruction {

  static constexpr bool debug = false;

  /**
   * Contains data required to grow two region around
   * the cut of two blocks of the partition.
   */
  struct BFSData {
    explicit BFSData(const HypernodeID num_nodes,
                     const HyperedgeID num_edges) :
      current_distance(0),
      last_queue_idx(0),
      queue(2),
      next_queue(2),
      visited_hn(num_nodes, false),
      visited_he(num_edges, false) { }

    void clearQueues();

    void reset();

    void pop_hypernode(HypernodeID& hn);

    void add_pins_of_hyperedge_to_queue(const HyperedgeID& he,
                                        const PartitionedHypergraph& phg,
                                        const ProblemStats& stats,
                                        const size_t max_bfs_distance);

    bool is_empty() const {
      return queue[0].empty() && queue[1].empty();
    }

    bool is_next_empty() const {
      return next_queue[0].empty() && next_queue[1].empty();
    }

    void swap_with_next_queue() {
      if ( !is_next_empty() ) {
        std::swap(queue, next_queue);
        ++current_distance;
      }
    }

    BlockPair blocks;
    size_t current_distance;
    size_t last_queue_idx;
    vec<parallel::scalable_queue<HypernodeID>> queue;
    vec<parallel::scalable_queue<HypernodeID>> next_queue;
    vec<bool> visited_hn;
    vec<bool> visited_he;
  };

  /**
   * Manages BFS data of all block pairs on which simultanously
   * grow regions. The search can pop hypernodes from the queues,
   * which are selected in round-robin fashion.
   */
  struct ConstructionData {
    explicit ConstructionData(const HypernodeID num_nodes,
                              const HyperedgeID num_edges) :
      _num_nodes(num_nodes),
      _num_edges(num_edges),
      used_slots(0),
      last_idx(0),
      bfs() { }

    void reset();

    void initialize(const vec<BlockPairCutHyperedges>& initial_cut_hes,
                    const ProblemStats& stats,
                    const PartitionedHypergraph& phg);

    void pop_hypernode(HypernodeID& hn, size_t& idx);

    void clearBlock(const PartitionID block);

    bool is_empty() {
      bool empty = true;
      for ( size_t i = 0; i < used_slots; ++i ) {
        empty &= ( bfs[i].is_empty() && bfs[i].is_next_empty() );
      }
      return empty;
    }

    const HypernodeID _num_nodes;
    const HyperedgeID _num_edges;
    size_t used_slots;
    size_t last_idx;
    vec<BFSData> bfs;
  };

 public:
  explicit ProblemConstruction(const Hypergraph& hg,
                               const Context& context) :
    _context(context),
    _vertex_ownership(hg.initialNumNodes(),
      CAtomic<SearchID>(QuotientGraph::INVALID_SEARCH_ID)),
    _local_data(hg.initialNumNodes(), hg.initialNumEdges()),
    _local_stats(hg.initialNumEdges(), context.partition.k) { }

  ProblemConstruction(const ProblemConstruction&) = delete;
  ProblemConstruction(ProblemConstruction&&) = delete;

  ProblemConstruction & operator= (const ProblemConstruction &) = delete;
  ProblemConstruction & operator= (ProblemConstruction &&) = delete;

  vec<HypernodeID> construct(const SearchID search_id,
                             QuotientGraph& quotient_graph,
                             AdvancedRefinerAdapter& refiner,
                             const PartitionedHypergraph& phg);

  void releaseNodes(const SearchID search_id,
                    const vec<HypernodeID>& nodes);

 private:
  bool acquire_vertex(const SearchID& search_id, const HypernodeID& hn) {
    ASSERT(static_cast<size_t>(hn) < _vertex_ownership.size());
    SearchID expected_id = QuotientGraph::INVALID_SEARCH_ID;
    return _vertex_ownership[hn].compare_exchange_strong(
      expected_id, search_id, std::memory_order_relaxed);
  }

  void release_vertex(SearchID search_id, const HypernodeID& hn) {
    ASSERT(static_cast<size_t>(hn) < _vertex_ownership.size());
    ASSERT(_vertex_ownership[hn] == search_id);
    _vertex_ownership[hn].compare_exchange_strong(
      search_id, QuotientGraph::INVALID_SEARCH_ID, std::memory_order_relaxed);
  }

  const Context& _context;

  // ! Keeps track which vertices belong to which search
  // ! Each vertex is only allowed to be part of one search at any time
  vec<CAtomic<SearchID>> _vertex_ownership;

  // ! Contains data required for BFS construction algorithm
  tbb::enumerable_thread_specific<ConstructionData> _local_data;

  // ! Contains statistic about the currently constructed problem
  tbb::enumerable_thread_specific<ProblemStats> _local_stats;
};

}  // namespace kahypar
