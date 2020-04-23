/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <algorithm>
#include <limits>
#include <stack>
#include <vector>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/sparse_map.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
template <typename TypeTraits = Mandatory,
          typename ScorePolicy = Mandatory,
          typename HeavyNodePenaltyPolicy = Mandatory,
          typename AcceptancePolicy = Mandatory>
class MultilevelVertexPairRater {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using LargeTmpRatingMap = ds::SparseMap<HypernodeID, RatingType>;
  using CacheEfficientRatingMap = ds::FixedSizeSparseMap<HypernodeID, RatingType>;
  using ThreadLocalLargeTmpRatingMap = tbb::enumerable_thread_specific<LargeTmpRatingMap>;
  using ThreadLocalCacheEfficientRatingMap = tbb::enumerable_thread_specific<CacheEfficientRatingMap>;

 private:
  static constexpr bool debug = false;

  class VertexPairRating {
   public:
    VertexPairRating(HypernodeID trgt, RatingType val, bool is_valid) :
      target(trgt),
      value(val),
      valid(is_valid) { }

    VertexPairRating() :
      target(std::numeric_limits<HypernodeID>::max()),
      value(std::numeric_limits<RatingType>::min()),
      valid(false) { }

    VertexPairRating(const VertexPairRating&) = delete;
    VertexPairRating & operator= (const VertexPairRating &) = delete;

    VertexPairRating(VertexPairRating&&) = default;
    VertexPairRating & operator= (VertexPairRating &&) = delete;

    HypernodeID target;
    RatingType value;
    bool valid;
  };

  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

 public:
  using Rating = VertexPairRating;

  MultilevelVertexPairRater(HyperGraph& hypergraph,
                           const Context& context) :
    _context(context),
    _current_num_nodes(hypergraph.initialNumNodes()),
    _local_cache_efficient_rating_map(0.0),
    _local_large_rating_map([&] {
      return construct_large_tmp_rating_map();
    }),
    _local_visited_representatives(hypergraph.initialNumNodes()),
    _already_matched(hypergraph.initialNumNodes()) { }

  MultilevelVertexPairRater(const MultilevelVertexPairRater&) = delete;
  MultilevelVertexPairRater & operator= (const MultilevelVertexPairRater &) = delete;

  MultilevelVertexPairRater(MultilevelVertexPairRater&&) = delete;
  MultilevelVertexPairRater & operator= (MultilevelVertexPairRater &&) = delete;

  VertexPairRating rate(const HyperGraph& hypergraph,
                        const HypernodeID u,
                        const parallel::scalable_vector<HypernodeID>& cluster_ids,
                        const parallel::scalable_vector<AtomicWeight>& cluster_weight,
                        const HypernodeWeight max_allowed_node_weight) {
    if ( ratingsFitIntoSmallSparseMap(hypergraph, u) ) {
      return rate(hypergraph, u, _local_cache_efficient_rating_map.local(),
        cluster_ids, cluster_weight, max_allowed_node_weight);
    } else {
      LargeTmpRatingMap& large_tmp_rating_map = _local_large_rating_map.local();
      large_tmp_rating_map.setMaxSize(_current_num_nodes);
      return rate(hypergraph, u, large_tmp_rating_map,
        cluster_ids, cluster_weight, max_allowed_node_weight);
    }
  }

  // ! Several threads will mark matches in parallel. However, since
  // ! we only set the corresponding value to true this function is
  // ! thread-safe.
  void markAsMatched(const HypernodeID original_id) {
    _already_matched.set(original_id, true);
  }

  // ! Note, this function is not thread safe
  void resetMatches() {
    _already_matched.reset();
  }

  void setCurrentNumberOfNodes(const HypernodeID current_num_nodes) {
    _current_num_nodes = current_num_nodes;
  }

 private:
  template<typename RatingMap>
  VertexPairRating rate(const HyperGraph& hypergraph,
                        const HypernodeID u,
                        RatingMap& tmp_ratings,
                        const parallel::scalable_vector<HypernodeID>& cluster_ids,
                        const parallel::scalable_vector<AtomicWeight>& cluster_weight,
                        const HypernodeWeight max_allowed_node_weight) {
    kahypar::ds::FastResetFlagArray<>& visited_representatives = _local_visited_representatives.local();
    const HypernodeID original_u_id = hypergraph.originalNodeID(u);
    const HypernodeWeight weight_u = cluster_weight[original_u_id];
    for ( const HyperedgeID& he : hypergraph.incidentEdges(u) ) {
      ASSERT(hypergraph.edgeSize(he) > 1, V(he));
      if ( hypergraph.edgeSize(he) < _context.partition.hyperedge_size_threshold ) {
        const RatingType score = ScorePolicy::score(hypergraph, he);
        for ( const HypernodeID& v : hypergraph.pins(he) ) {
          const HypernodeID original_v_id = hypergraph.originalNodeID(v);
          const HypernodeID representative = cluster_ids[original_v_id];
          ASSERT(representative < hypergraph.initialNumNodes());
          if ( !visited_representatives[representative] ) {
            tmp_ratings[representative] += score;
            visited_representatives.set(representative, true);
          }
        }
        visited_representatives.reset();
      }
    }

    int cpu_id = sched_getcpu();
    const PartitionID community_u_id = hypergraph.communityID(u);
    RatingType max_rating = std::numeric_limits<RatingType>::min();
    HypernodeID target = std::numeric_limits<HypernodeID>::max();
    HypernodeID target_id = std::numeric_limits<HypernodeID>::max();
    for (auto it = tmp_ratings.end() - 1; it >= tmp_ratings.begin(); --it) {
      const HypernodeID tmp_target_id = it->key;
      const HypernodeID tmp_target = hypergraph.globalNodeID(tmp_target_id);
      const HypernodeWeight target_weight = cluster_weight[tmp_target_id];

      if ( tmp_target != u && weight_u + target_weight <= max_allowed_node_weight ) {
        HypernodeWeight penalty = HeavyNodePenaltyPolicy::penalty(weight_u, target_weight);
        penalty = penalty == 0 ? std::max(std::max(weight_u, target_weight), 1) : penalty;
        const RatingType tmp_rating = it->value / static_cast<double>(penalty);

        DBG << "r(" << u << "," << tmp_target << ")=" << tmp_rating;
        if ( community_u_id == hypergraph.communityID(tmp_target) &&
            AcceptancePolicy::acceptRating(tmp_rating, max_rating,
                                           target_id, tmp_target_id,
                                           cpu_id, _already_matched) ) {
          max_rating = tmp_rating;
          target_id = tmp_target_id;
          target = tmp_target;
        }
      }
    }

    VertexPairRating ret;
    if (max_rating != std::numeric_limits<RatingType>::min()) {
      ASSERT(target != std::numeric_limits<HypernodeID>::max(), "invalid contraction target");
      ret.value = max_rating;
      ret.target = target;
      ret.valid = true;
    }
    tmp_ratings.clear();
    return ret;
  }


  inline bool ratingsFitIntoSmallSparseMap(const HyperGraph& hypergraph,
                                           const HypernodeID u)  {
    // In case the current number of nodes is smaller than size
    // of the cache-efficient sparse map, the large tmp rating map
    // consumes less memory
    if ( _current_num_nodes < CacheEfficientRatingMap::MAP_SIZE ) {
      return false;
    }
    // Compute estimation for the upper bound of neighbors of u
    HypernodeID ub_neighbors_u = 0;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(u) ) {
      const HypernodeID edge_size = hypergraph.edgeSize(he);
      // Ignore large hyperedges
      ub_neighbors_u += edge_size < _context.partition.hyperedge_size_threshold ? edge_size : 0;
      // If the number of estimated neighbors is greater than MAP_SIZE / 3, we
      // use the large sparse map. The division by 3 also ensures that the fill grade
      // of the small sparse map would be small enough such that linear probing
      // is fast.
      if ( ub_neighbors_u > CacheEfficientRatingMap::MAP_SIZE / 3UL ) {
        return false;
      }
    }
    return true;
  }

  LargeTmpRatingMap construct_large_tmp_rating_map() {
    return LargeTmpRatingMap(_current_num_nodes);
  }

  const Context& _context;
  HypernodeID _current_num_nodes;
  ThreadLocalCacheEfficientRatingMap _local_cache_efficient_rating_map;
  ThreadLocalLargeTmpRatingMap _local_large_rating_map;
  ThreadLocalFastResetFlagArray _local_visited_representatives;
  kahypar::ds::FastResetFlagArray<> _already_matched;
};
}  // namespace mt_kahypar
