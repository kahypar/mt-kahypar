/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
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

#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>
#include <mt-kahypar/partition/refinement/async_refiners_common.h>

#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"

namespace mt_kahypar {


template<typename FMStrategy>
class AsyncKWayFM {
private:
    static constexpr bool enable_heavy_assert = true;

public:
  explicit AsyncKWayFM(const Context& context, HypernodeID numNodes, AsyncFMSharedData& sharedData,
                       ds::GroupLockManager * const lock_manager) :
      context(context),
      k(context.partition.k),
      deltaPhg(context.partition.k),
      neighborDeduplicator(numNodes, 0),
      fm_strategy(context, numNodes, sharedData, runStats),
      sharedData(sharedData),
      uncontraction_locks(lock_manager),
      contraction_group_id(ds::invalidGroupID),
      _km1_delta(0),
      attempted_to_move(numNodes),
      num_nodes_in_pq(0),
      max_num_nodes_in_pq(context.refinement.fm.async_max_num_nodes_in_pq),
      max_num_moves(context.refinement.fm.async_max_num_moves)
          {}

  // ! Finds a sequence of moves, applies the best prefix and returns the actual km1 improvement of that prefix.
  // ! Returns 0 if no moves were performed.
  Gain findMoves(PartitionedHypergraph& phg, const vec<HypernodeID>& refinement_nodes);

  void memoryConsumption(utils::MemoryTreeNode* parent) const ;

  void resetForGroup(const ds::ContractionGroupID groupID) {
    contraction_group_id = groupID;
  }

  FMStats stats;

private:

  // ! Performs localized FM local search on the delta partitioned hypergraph.
  // ! Moves made by this search are not immediately visible to other concurrent local searches.
  // ! The best prefix of moves is applied to the global partitioned hypergraph after the search finishes.
  //void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg, FMSharedData& sharedData);


  void internalFindMoves(PartitionedHypergraph& phg);

  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void acquireOrUpdateNeighbors(PHG& phg, const Move& move);


  // ! Makes moves applied on delta hypergraph visible on the global partitioned hypergraph.
  std::pair<Gain, size_t> applyBestLocalPrefixToSharedPartition(PartitionedHypergraph& phg,
                                                                const size_t best_index_locally_observed,
                                                                const Gain best_improvement_locally_observed,
                                                                bool apply_all_moves);

  void revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex);

  void releaseMoveLocksForLocalMovedNodes();

  void lockUncontractionLockWithWaiting(const HypernodeID& hn);

 private:

  const Context& context;

  // ! Number of blocks
  PartitionID k;

  // ! Moves found by this FM search
  vec< Move > localMoves;

  // ! Wrapper around the global partitioned hypergraph, that allows
  // ! to perform moves non-visible for other local searches
  ds::DeltaPartitionedHypergraph<PartitionedHypergraph> deltaPhg;

  // ! Used after a move. Stores whether a neighbor of the just moved vertex has already been updated.
  vec<HypernodeID> neighborDeduplicator;
  HypernodeID deduplicationTime = 0;

  // ! Stores hyperedges whose pins's gains may have changed after vertex move
  vec<HyperedgeID> edgesWithGainChanges;

  FMStats runStats;

  FMStrategy fm_strategy;

  AsyncFMSharedData& sharedData;

  // ! Used to hold a lock on nodes while they are being moved on the hypergraph.
  // ! Prevents concurrent uncontractions on those moved nodes.
  ds::GroupLockManager* const uncontraction_locks;

  // ! ID of the contraction group that the current local search is based on. Identifies a local search run.
  ds::ContractionGroupID contraction_group_id;

  // ! Used to keep track of the total gain of all moves in this local search run
  Gain _km1_delta;

  // ! Nodes which have been attempted to move in this local search run
  kahypar::ds::FastResetFlagArray<> attempted_to_move;

  size_t num_nodes_in_pq;
  const size_t max_num_nodes_in_pq;

  const size_t max_num_moves;


};

}