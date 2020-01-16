/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

template<typename HyperGraph>
class ZeroGainCache {

 private:
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }


  static constexpr PartitionID kInvalidPartition = -1;
  static constexpr size_t kInvalidPosition = std::numeric_limits<size_t>::max();

  struct VertexCacheEntry {
    VertexCacheEntry(const PartitionID k) :
      from(kInvalidPartition),
      valid_to(k, false) { }

    void reset(const PartitionID block) {
      ASSERT(from != block);
      from = block;
      std::fill(valid_to.begin(), valid_to.end(), false);
    }

    PartitionID from;
    parallel::scalable_vector<bool> valid_to;
  };

  template<typename T>
  using vector = parallel::scalable_vector<T>;
  using Cache = vector<vector<vector<HypernodeID>>>;
  using CachePtr = vector<VertexCacheEntry>;

 public:
  explicit ZeroGainCache(const HypernodeID num_nodes,
                         const Context& context) :
    _context(context),
    _cache(context.partition.k, vector<vector<HypernodeID>>(context.partition.k)),
    _cache_entry(num_nodes, VertexCacheEntry(context.partition.k)) { }

  void insert(const HyperGraph& hg,
              const HypernodeID hn,
              const PartitionID from,
              const PartitionID to) {
    const HypernodeID original_id = hg.originalNodeID(hn);
    ASSERT(original_id < _cache_entry.size());
    ASSERT(from != kInvalidPartition && from < _context.partition.k);
    ASSERT(to != kInvalidPartition && to < _context.partition.k);

    VertexCacheEntry& cache_entry = _cache_entry[original_id];
    // In case there exists a cache entry for vertex hn in an other block,
    // the block of the vertex changed and the current entry is therefore
    // invalid => reset corresponding entry
    if ( cache_entry.from != from ) {
      cache_entry.reset(from);
    }
    ASSERT(cache_entry.from == from);

    // If not inserted before, than we insert zero gain move
    if ( !cache_entry.valid_to[to] ) {
      _cache[from][to].emplace_back(hn);
      cache_entry.valid_to[to] = true;
    }
  }

  bool isMovePossible(HyperGraph& hg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to) {
    ASSERT(from != kInvalidPartition && from < _context.partition.k);
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    const HypernodeWeight weight_buffer =
      _context.partition.max_part_weights[from] - hg.localPartWeight(from);
    if ( weight_buffer < 0 ) {
      // Rebalancing not possible since balance constraint already violated
      return false;
    }
    // The weight of all vertices used for rebalancing must be greater or
    // equal than the hypernode weight, but smaller or equal than
    // the upper bound. If the rebalancing weight would be greater
    // than the upper bound we would violate the balance constraint
    // in the block 'from'
    const HypernodeWeight hn_weight = hg.nodeWeight(hn);
    const HypernodeWeight upper_bound = weight_buffer + hn_weight;
    HypernodeWeight rebalancing_weight = 0;

    size_t end = _cache[to][from].size();
    for ( size_t pos = 0; pos < end; ++pos ) {
      const HypernodeID u = _cache[to][from][pos];
      if ( hg.partID(u) != to ) {
        // Lazy removal of vertices
        // In that case, cache entry indicates that moving hypernode u
        // from block to to block from is a zero gain move, but u is no
        // longer part of part to.
        std::swap(_cache[to][from][pos--], _cache[to][from][--end]);
        _cache[to][from].pop_back();
        continue;
      }

      const HypernodeWeight u_weight = hg.nodeWeight(u);
      if ( rebalancing_weight + u_weight <= upper_bound ) {
        rebalancing_weight += u_weight;
        if ( rebalancing_weight >= hn_weight ) {
          break;
        }
      }
    }

    return hn_weight <= rebalancing_weight && rebalancing_weight <= upper_bound;
  }

  bool performMove(HyperGraph& hg,
                   const HypernodeID hn,
                   const PartitionID from,
                   const PartitionID to,
                   Gain& delta) {
    ASSERT(from != kInvalidPartition && from < _context.partition.k);
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    auto objective_delta = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
                              if ( _context.partition.objective == kahypar::Objective::cut ) {
                                delta += HyperGraph::cutDelta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
                              } else {
                                delta += HyperGraph::km1Delta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
                              }
                           };

    const HypernodeWeight weight_buffer =
      _context.partition.max_part_weights[from] - hg.localPartWeight(from);
    if ( weight_buffer < 0 ) {
      // Rebalancing not possible since balance constraint already violated
      return false;
    }
    // The weight of all vertices used for rebalancing must be greater or
    // equal than the hypernode weight, but smaller or equal than
    // the upper bound. If the rebalancing weight would be greater
    // than the upper bound we would violate the balance constraint
    // in the block 'from'
    const HypernodeWeight hn_weight = hg.nodeWeight(hn);
    const HypernodeWeight upper_bound = weight_buffer + hn_weight;
    HypernodeWeight rebalancing_weight = 0;

    parallel::scalable_vector<HypernodeID> moved_hypernodes;
    size_t end = _cache[to][from].size();
    for ( size_t pos = 0; pos < end; ++pos ) {
      const HypernodeID u = _cache[to][from][pos];
      if ( hg.partID(u) != to ) {
        // Lazy removal of vertices
        // In that case, cache entry indicates that moving hypernode u
        // from block to to block from is a zero gain move, but u is no
        // longer part of part to.
        std::swap(_cache[to][from][pos--], _cache[to][from][--end]);
        _cache[to][from].pop_back();
        continue;
      }

      const HypernodeWeight u_weight = hg.nodeWeight(u);
      if ( rebalancing_weight + u_weight <= upper_bound ) {
        Gain delta_before = delta;
        if ( hg.changeNodePart(u, to, from, objective_delta) ) {
          if ( delta - delta_before == 0 ) {
            // In case the move of vertex u from block 'to' to block
            // 'from' succeeds and it is still a zero gain move, we
            // add its weight to the rebalancing weight
            rebalancing_weight += u_weight;
            moved_hypernodes.emplace_back(u);
          } else {
            // In case, vertex u is no longer a zero gain move, we
            // revert the move immediatly
            hg.changeNodePart(u, from, to, objective_delta);
          }

          // Remove vertex u from cache, because it was either moved
          // successfully to an other block or it is no longer a
          // zero gain move.
          const HypernodeID original_id = hg.originalNodeID(u);
          ASSERT(original_id < _cache_entry.size());
          _cache_entry[original_id].valid_to[from] = false;
          std::swap(_cache[to][from][pos--], _cache[to][from][--end]);
          _cache[to][from].pop_back();

          if ( rebalancing_weight >= hn_weight ) {
            break;
          }
        }
      }
    }

    bool success = hn_weight <= rebalancing_weight && rebalancing_weight <= upper_bound;
    if ( !success || !hg.changeNodePart(hn, from, to, objective_delta) ) {
      // In case, we were not successful, because either the rebalancing weight
      // is smaller than the vertex weight or we would violate the balance
      // constraint, we immediatly revert all moves
      revertMoves(hg, moved_hypernodes, to, from, delta);
    }
    return success;
  }

 private:
  FRIEND_TEST(AZeroGainCache, InsertsAZeroGainMove);
  FRIEND_TEST(AZeroGainCache, InsertsTwoZeroGainMoves);
  FRIEND_TEST(AZeroGainCache, InsertsTwoZeroGainMovesToDifferentBlocks);
  FRIEND_TEST(AZeroGainCache, InsertsTwoZeroGainMovesForSameVertex1);
  FRIEND_TEST(AZeroGainCache, InsertsTwoZeroGainMovesForSameVertex2);
  FRIEND_TEST(AZeroGainCache, ReinsertZeroGainMoveAfterChangeNodePart1);
  FRIEND_TEST(AZeroGainCache, ReinsertZeroGainMoveAfterChangeNodePart2);

  void revertMoves(HyperGraph& hg,
                   const parallel::scalable_vector<HypernodeID>& moved_hypernodes,
                   const PartitionID from,
                   const PartitionID to,
                   Gain& delta) {
    ASSERT(from != kInvalidPartition && from < _context.partition.k);
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    auto objective_delta = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
                              if ( _context.partition.objective == kahypar::Objective::cut ) {
                                delta += HyperGraph::cutDelta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
                              } else {
                                delta += HyperGraph::km1Delta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
                              }
                           };

    for ( const HypernodeID hn : moved_hypernodes ) {
      const PartitionID part_id = hg.partID(hn);
      if ( part_id == from && hg.changeNodePart(hn, from, to, objective_delta) ) {
        // Revert move and reinsert vertex hn to cache
        const HypernodeID original_id = hg.originalNodeID(hn);
        ASSERT(original_id < _cache_entry.size());
        ASSERT(!_cache_entry[original_id].valid_to[from]);
        _cache[to][from].emplace_back(hn);
        _cache_entry[original_id].valid_to[from] = true;
      }
    }
  }

  const Context _context;
  Cache _cache;
  CachePtr _cache_entry;
};

}  // namespace kahypar
