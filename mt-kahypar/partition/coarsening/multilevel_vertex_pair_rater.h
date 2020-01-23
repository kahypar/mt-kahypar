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
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/concurrent_union_find.h"

namespace mt_kahypar {
template <typename TypeTraits = Mandatory,
          typename ScorePolicy = Mandatory,
          typename HeavyNodePenaltyPolicy = Mandatory,
          typename AcceptancePolicy = Mandatory>
class MultilevelVertexPairRater {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using UnionFind = ds::ConcurrentUnionFind<HyperGraph>;
  using TmpRatingMap = kahypar::ds::SparseMap<HypernodeID, RatingType>;
  using ThreadLocalTmpRatingMap = tbb::enumerable_thread_specific<TmpRatingMap>;

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

 public:
  using Rating = VertexPairRating;

  MultilevelVertexPairRater(HyperGraph& hypergraph,
                           const Context& context,
                           UnionFind& uf) :
    _hg(hypergraph),
    _context(context),
    _uf(uf),
    _local_tmp_ratings(_hg.initialNumNodes()),
    _already_matched(_hg.initialNumNodes()) { }

  MultilevelVertexPairRater(const MultilevelVertexPairRater&) = delete;
  MultilevelVertexPairRater & operator= (const MultilevelVertexPairRater &) = delete;

  MultilevelVertexPairRater(MultilevelVertexPairRater&&) = delete;
  MultilevelVertexPairRater & operator= (MultilevelVertexPairRater &&) = delete;

  VertexPairRating rate(const HypernodeID u) {
    TmpRatingMap& tmp_ratings = _local_tmp_ratings.local();
    const HypernodeID original_u_id = _hg.originalNodeID(u);
    const HypernodeWeight weight_u = _uf.weight(original_u_id);
    for ( const HyperedgeID& he : _hg.incidentEdges(u) ) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      if ( _hg.edgeSize(he) < _context.partition.hyperedge_size_threshold ) {
        const RatingType score = ScorePolicy::score(_hg, he);
        for ( const HypernodeID& v : _hg.pins(he) ) {
          const HypernodeID original_v_id = _hg.originalNodeID(v);
          if ( u != v && belowThresholdNodeWeight(original_u_id,
               original_v_id, weight_u, _uf.weight(original_v_id)) ) {
            ASSERT(original_v_id < tmp_ratings.size());
            tmp_ratings[original_v_id] += score;
          }
        }
      }
    }

    int cpu_id = sched_getcpu();
    const PartitionID community_u_id = _hg.communityID(u);
    RatingType max_rating = std::numeric_limits<RatingType>::min();
    HypernodeID target = std::numeric_limits<HypernodeID>::max();
    HypernodeID target_id = std::numeric_limits<HypernodeID>::max();
    for (auto it = tmp_ratings.end() - 1; it >= tmp_ratings.begin(); --it) {
      const HypernodeID tmp_target_id = it->key;
      const HypernodeID tmp_target = _hg.globalNodeID(tmp_target_id);
      const HypernodeWeight target_weight = _uf.isSameSet(original_u_id, tmp_target_id) ? 1 : _uf.weight(tmp_target);

      HypernodeWeight penalty = HeavyNodePenaltyPolicy::penalty(weight_u, target_weight);
      penalty = penalty == 0 ? std::max(std::max(weight_u, target_weight), 1) : penalty;
      const RatingType tmp_rating = it->value / static_cast<double>(penalty);

      DBG << "r(" << u << "," << tmp_target << ")=" << tmp_rating;
      if ( community_u_id == _hg.communityID(tmp_target) &&
           AcceptancePolicy::acceptRating(tmp_rating, max_rating,
                                          target_id, tmp_target_id,
                                          cpu_id, _already_matched) ) {
        max_rating = tmp_rating;
        target_id = tmp_target_id;
        target = tmp_target;
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

  // ! Several threads will mark matches in parallel. However, since
  // ! we only set the corresponding value to true this function is
  // ! thread-safe.
  void markAsMatched(const HypernodeID hn) {
    _already_matched.set(_hg.originalNodeID(hn), true);
  }

  // ! Note, this function is not thread safe
  void resetMatches() {
    _already_matched.reset();
  }

 private:
  inline bool belowThresholdNodeWeight(const HypernodeID original_u_id,
                                       const HypernodeID original_v_id,
                                       const HypernodeWeight weight_u,
                                       const HypernodeWeight weight_v) const {
    // In case, if u and v are already in the same set (which means that they are
    // already contracted togehter), the weight is always below the threshold, otherwise
    // we perform an explicit check
    return !_uf.isSameSet(original_u_id, original_v_id) ?
      weight_v + weight_u <= thresholdNodeWeight(original_v_id) : true;
  }

  inline HypernodeWeight thresholdNodeWeight(const HypernodeID v) const {
    return _uf.containsHighDegreeVertex(v) ?
     _context.coarsening.max_allowed_high_degree_node_weight :
     _context.coarsening.max_allowed_node_weight;
  }

  HyperGraph& _hg;
  const Context& _context;
  UnionFind& _uf;
  ThreadLocalTmpRatingMap _local_tmp_ratings;
  kahypar::ds::FastResetFlagArray<> _already_matched;
};
}  // namespace mt_kahypar
