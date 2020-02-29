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

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
template <typename TypeTraits = Mandatory,
          typename ScorePolicy = Mandatory,
          typename HeavyNodePenaltyPolicy = Mandatory,
          typename AcceptancePolicy = Mandatory>
class CommunityVertexPairRater {
  using HyperGraph = typename TypeTraits::HyperGraph;

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
  using HypernodeMapping = std::shared_ptr<std::vector<HypernodeID> >;

  CommunityVertexPairRater(HyperGraph& hypergraph,
                           const Context& context,
                           const PartitionID community_id,
                           const parallel::scalable_vector<HypernodeID>& community_node_mapping) :
    _hg(hypergraph),
    _context(context),
    _community_id(community_id),
    _community_node_mapping(community_node_mapping),
    _tmp_ratings(_hg.numCommunityHypernodes(_community_id)),
    _already_matched(_hg.numCommunityHypernodes(_community_id)) { }

  CommunityVertexPairRater(const CommunityVertexPairRater&) = delete;
  CommunityVertexPairRater & operator= (const CommunityVertexPairRater &) = delete;

  CommunityVertexPairRater(CommunityVertexPairRater&&) = delete;
  CommunityVertexPairRater & operator= (CommunityVertexPairRater &&) = delete;

  VertexPairRating rate(const HypernodeID u) {
    ASSERT(_hg.communityID(u) == _community_id);
    const HypernodeWeight weight_u = _hg.nodeWeight(u);
    for (const HyperedgeID& he : _hg.multiPinIncidentEdges(u, _community_id)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      if (_hg.edgeSize(he, _community_id) < _context.partition.hyperedge_size_threshold) {
        const RatingType score = ScorePolicy::score(_hg, he, _community_id);
        for (const HypernodeID& v : _hg.pins(he, _community_id)) {
          if (u != v && belowThresholdNodeWeight(weight_u, _hg.nodeWeight(v))) {
            ASSERT(_hg.communityID(v) == _community_id);
            ASSERT(_community_node_mapping[_hg.communityNodeId(v)] == v);
            _tmp_ratings[_hg.communityNodeId(v)] += score;
          }
        }
      }
    }

    int cpu_id = sched_getcpu();
    RatingType max_rating = std::numeric_limits<RatingType>::min();
    HypernodeID target = std::numeric_limits<HypernodeID>::max();
    HypernodeID community_target = std::numeric_limits<HypernodeID>::max();
    for (auto it = _tmp_ratings.end() - 1; it >= _tmp_ratings.begin(); --it) {
      const HypernodeID tmp_community_target = it->key;
      const HypernodeID tmp_target = _community_node_mapping[tmp_community_target];
      ASSERT(_hg.communityID(tmp_target) == _community_id);
      ASSERT(tmp_community_target == _hg.communityNodeId(tmp_target));
      const HypernodeWeight target_weight = _hg.nodeWeight(tmp_target);

      HypernodeWeight penalty = HeavyNodePenaltyPolicy::penalty(weight_u, target_weight);
      penalty = penalty == 0 ? std::max(std::max(weight_u, target_weight), 1) : penalty;
      const RatingType tmp_rating = it->value / static_cast<double>(penalty);

      DBG << "r(" << u << "," << tmp_target << ")=" << tmp_rating;
      if (AcceptancePolicy::acceptRating(tmp_rating, max_rating,
                                         community_target, tmp_community_target,
                                         cpu_id, _already_matched)) {
        max_rating = tmp_rating;
        community_target = tmp_community_target;
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
    _tmp_ratings.clear();
    DBG << "rating=(" << ret.value << "," << ret.target << "," << ret.valid << ")";
    return ret;
  }

  void markAsMatched(const HypernodeID hn) {
    _already_matched.set(_hg.communityNodeId(hn), true);
  }

  void resetMatches() {
    _already_matched.reset();
  }

 private:
  inline bool belowThresholdNodeWeight(const HypernodeWeight weight_u,
                                       const HypernodeWeight weight_v) const {
    return weight_v + weight_u <= _context.coarsening.max_allowed_node_weight;
  }

  HyperGraph& _hg;
  const Context& _context;
  const PartitionID _community_id;
  const parallel::scalable_vector<HypernodeID>& _community_node_mapping;
  kahypar::ds::SparseMap<HypernodeID, RatingType> _tmp_ratings;
  kahypar::ds::FastResetFlagArray<> _already_matched;
};
}  // namespace mt_kahypar
