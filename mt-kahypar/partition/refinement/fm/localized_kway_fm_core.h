/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"

namespace mt_kahypar {


template<typename TypeTraits, typename GainTypes>
class LocalizedKWayFM {

  static constexpr size_t MAP_SIZE_LARGE = 16384;
  static constexpr size_t MAP_SIZE_MOVE_DELTA = 8192;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using DeltaGainCache = typename GainTypes::DeltaGainCache;
  using DeltaPartitionedHypergraph = typename PartitionedHypergraph::template DeltaPartition<DeltaGainCache::requires_connectivity_set>;
  using AttributedGains = typename GainTypes::AttributedGains;

public:
  explicit LocalizedKWayFM(const Context& context,
                           const HypernodeID numNodes,
                           FMSharedData& sharedData,
                           GainCache& gainCache) :
    context(context),
    thisSearch(0),
    deltaPhg(context),
    neighborDeduplicator(numNodes, 0),
    fm_strategy(context, sharedData, runStats),
    gain_cache(gainCache),
    delta_gain_cache(gainCache),
    sharedData(sharedData) {
    const bool top_level = context.type == ContextType::main;
    delta_gain_cache.initialize(top_level ? MAP_SIZE_LARGE : MAP_SIZE_MOVE_DELTA);
  }


  bool findMoves(PartitionedHypergraph& phg, size_t taskID, size_t numSeeds);

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

  void changeNumberOfBlocks(const PartitionID new_k);

  FMStats stats;

private:

  // ! Performs localized FM local search on the delta partitioned hypergraph.
  // ! Moves made by this search are not immediately visible to other concurrent local searches.
  // ! The best prefix of moves is applied to the global partitioned hypergraph after the search finishes.
  //void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg, FMSharedData& sharedData);


  template<bool use_delta, bool has_fixed_vertices>
  void internalFindMoves(PartitionedHypergraph& phg);

  template<bool has_fixed_vertices, typename PHG, typename CACHE>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void acquireOrUpdateNeighbors(PHG& phg, CACHE& gain_cache, const Move& move);


  // ! Makes moves applied on delta hypergraph visible on the global partitioned hypergraph.
  std::pair<Gain, size_t> applyBestLocalPrefixToSharedPartition(PartitionedHypergraph& phg,
                                                                const size_t best_index_locally_observed,
                                                                const Gain best_improvement_locally_observed,
                                                                const bool apply_delta_improvement);

  // ! Rollback to the best improvement found during local search in case we applied moves
  // ! directly on the global partitioned hypergraph.
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex);

 private:

  const Context& context;

  // ! Unique search id associated with the current local search
  SearchID thisSearch;

  // ! Local data members required for one localized search run
  //FMLocalData localData;
  vec< std::pair<Move, MoveID> > localMoves;

  // ! Wrapper around the global partitioned hypergraph, that allows
  // ! to perform moves non-visible for other local searches
  DeltaPartitionedHypergraph deltaPhg;

  // ! Used after a move. Stores whether a neighbor of the just moved vertex has already been updated.
  vec<HypernodeID> neighborDeduplicator;
  HypernodeID deduplicationTime = 0;

  // ! Stores hyperedges whose pins's gains may have changed after vertex move
  vec<HyperedgeID> edgesWithGainChanges;

  FMStats runStats;

  GainCacheStrategy fm_strategy;

  GainCache& gain_cache;

  DeltaGainCache delta_gain_cache;

  FMSharedData& sharedData;
};

}