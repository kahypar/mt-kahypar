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

    void clearQueue(const PartitionID block);

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

 public:
  explicit ProblemConstruction(const Hypergraph& hg,
                               const Context& context) :
    _context(context),
    _vertex_ownership(hg.initialNumNodes(),
      CAtomic<SearchID>(QuotientGraph::INVALID_SEARCH_ID)),
    _local_bfs(hg.initialNumNodes(), hg.initialNumEdges()),
    _local_stats(hg.initialNumEdges(), context.partition.k) { }

  ProblemConstruction(const ProblemConstruction&) = delete;
  ProblemConstruction(ProblemConstruction&&) = delete;

  ProblemConstruction & operator= (const ProblemConstruction &) = delete;
  ProblemConstruction & operator= (ProblemConstruction &&) = delete;

  Subhypergraph construct(const SearchID search_id,
                          QuotientGraph& quotient_graph,
                          AdvancedRefinerAdapter& refiner,
                          const PartitionedHypergraph& phg);

  void releaseNodes(const SearchID search_id,
                    const Subhypergraph& sub_hg);

 private:
  bool acquire_vertex(const SearchID& search_id, const HypernodeID& hn) {
    ASSERT(static_cast<size_t>(hn) < _vertex_ownership.size());
    SearchID expected_id = QuotientGraph::INVALID_SEARCH_ID;
    return _vertex_ownership[hn].compare_exchange_strong(
      expected_id, search_id, std::memory_order_relaxed);
  }

  void release_vertex(SearchID search_id, const HypernodeID& hn) {
    ASSERT(static_cast<size_t>(hn) < _vertex_ownership.size());
    ASSERT(_context.refinement.advanced.use_overlapping_searches ||
           _vertex_ownership[hn] == search_id);
    _vertex_ownership[hn].compare_exchange_strong(
      search_id, QuotientGraph::INVALID_SEARCH_ID, std::memory_order_relaxed);
  }

  const Context& _context;

  // ! Keeps track which vertices belong to which search
  // ! Each vertex is only allowed to be part of one search at any time
  vec<CAtomic<SearchID>> _vertex_ownership;

  // ! Contains data required for BFS construction algorithm
  tbb::enumerable_thread_specific<BFSData> _local_bfs;

  // ! Contains statistic about the currently constructed problem
  tbb::enumerable_thread_specific<ProblemStats> _local_stats;
};

}  // namespace kahypar
