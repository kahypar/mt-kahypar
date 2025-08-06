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


#include "mt-kahypar/partition/refinement/flows/quotient_graph.h"

#include <queue>

#include <tbb/parallel_sort.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/sparse_map.h"


namespace mt_kahypar {

void QuotientGraphEdge::add_hyperedge(const HyperedgeID he, const HyperedgeWeight weight) {
  cut_hes.push_back(he);
  cut_he_weight += weight;
  ++num_cut_hes;
}

void QuotientGraphEdge::reset() {
  cut_hes.clear();
  ownership.store(INVALID_SEARCH_ID, std::memory_order_relaxed);
  is_in_queue.store(false, std::memory_order_relaxed);
  num_cut_hes.store(0, std::memory_order_relaxed);
  cut_he_weight.store(0, std::memory_order_relaxed);
}


QuotientGraph::QuotientGraph(const HyperedgeID num_hyperedges, const Context& context) :
  _context(context),
  _initial_num_edges(num_hyperedges),
  _current_num_edges(kInvalidHyperedge),
  _quotient_graph(context.partition.k,
    vec<QuotientGraphEdge>(context.partition.k)),
  _register_search_lock(),
  _num_active_searches(0),
  _searches() {
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].blocks.i = i;
      _quotient_graph[i][j].blocks.j = j;
    }
  }
}

BlockPair QuotientGraph::getBlockPair(const SearchID search_id) const {
  ASSERT(search_id < _searches.size());
  return _searches[search_id].blocks;
}

size_t QuotientGraph::getRound(const SearchID search_id) const {
  ASSERT(search_id < _searches.size());
  return _searches[search_id].round;
}

vec<vec<QuotientGraphEdge>>& QuotientGraph::getGraph() {
  return _quotient_graph;
}

size_t QuotientGraph::numActiveSearches() const {
  return _num_active_searches;
}

template<typename PartitionedHypergraph>
void QuotientGraph::addNewCutHyperedges(const PartitionedHypergraph& phg, const vec<std::pair<HyperedgeID, PartitionID>>& new_cut_hes) {
  for ( const auto& [he, block] : new_cut_hes ) {
    ASSERT(block != kInvalidPartition);
    ASSERT(phg.pinCountInPart(he, block) > 0);
    // Add hyperedge he as a cut hyperedge to each block pair that contains 'block'
    for ( const PartitionID& other_block : phg.connectivitySet(he) ) {
      if ( other_block != block ) {
        _quotient_graph[std::min(block, other_block)][std::max(block, other_block)]
          .add_hyperedge(he, phg.edgeWeight(he));
      }
    }
  }
}

void QuotientGraph::finalizeConstruction(const SearchID search_id) {
  ASSERT(search_id < _searches.size());
  _searches[search_id].is_finalized = true;
  const BlockPair& blocks = _searches[search_id].blocks;
  _quotient_graph[blocks.i][blocks.j].release(search_id);
}

void QuotientGraph::finalizeSearch(const SearchID search_id, const HyperedgeWeight total_improvement) {
  ASSERT(search_id < _searches.size());
  ASSERT(_searches[search_id].is_finalized);

  const BlockPair& blocks = _searches[search_id].blocks;
  QuotientGraphEdge& qg_edge = _quotient_graph[blocks.i][blocks.j];
  if ( total_improvement > 0 ) {
    // If the search improves the quality of the partition, we reinsert
    // all hyperedges that were used by the search and are still cut.
    ++qg_edge.num_improvements_found;
    qg_edge.total_improvement += total_improvement;
  }
  --_num_active_searches;
}

template<typename PartitionedHypergraph>
void QuotientGraph::initialize(const PartitionedHypergraph& phg) {
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
  _current_num_edges = local_num_hes.combine(std::plus<HyperedgeID>());
}

HyperedgeWeight QuotientGraph::getCutHyperedgeWeightOfBlockPair(const PartitionID i, const PartitionID j) const {
  ASSERT(i < j);
  ASSERT(0 <= i && i < _context.partition.k);
  ASSERT(0 <= j && j < _context.partition.k);
  return _quotient_graph[i][j].cut_he_weight;
}

void QuotientGraph::changeNumberOfBlocks(const PartitionID new_k) {
  // Reset improvement history as the number of blocks had changed
  for ( size_t i = 0; i < _quotient_graph.size(); ++i ) {
    for ( size_t j = 0; j < _quotient_graph.size(); ++j ) {
      _quotient_graph[i][j].num_improvements_found.store(0, std::memory_order_relaxed);
      _quotient_graph[i][j].total_improvement.store(0, std::memory_order_relaxed);
    }
  }

  if ( static_cast<size_t>(new_k) > _quotient_graph.size() ) {
    _quotient_graph.clear();
    _quotient_graph.assign(new_k, vec<QuotientGraphEdge>(new_k));
  }
}

void QuotientGraph::resetQuotientGraphEdges() {
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    for ( PartitionID j = i + 1; j < _context.partition.k; ++j ) {
      _quotient_graph[i][j].reset();
    }
  }
}

bool QuotientGraph::isInputHypergraph() const {
  return _current_num_edges == _initial_num_edges;
}


namespace {
#define ADD_NEW_CUT_HYPEREDGE(X) void QuotientGraph::addNewCutHyperedges(const X& phg, const vec<std::pair<HyperedgeID, PartitionID>>& new_cut_hes)
#define INITIALIZE(X) void QuotientGraph::initialize(const X& phg)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(ADD_NEW_CUT_HYPEREDGE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(INITIALIZE)

} // namespace mt_kahypar
