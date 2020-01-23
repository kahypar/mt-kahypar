/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <string>

#include "tbb/concurrent_queue.h"
#include "tbb/task_group.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <typename TypeTraits,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class MultilevelCoarsenerT : public ICoarsener {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  using UnionFind = ds::ConcurrentUnionFind<HyperGraph>;
  using Rater = MultilevelVertexPairRater<TypeTraits,
                                          ScorePolicy,
                                          HeavyNodePenaltyPolicy,
                                          AcceptancePolicy>;
  using Rating = typename Rater::Rating;


  static constexpr bool debug = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsenerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _uf(hypergraph),
    _rater(hypergraph, context, _uf),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) { }

  MultilevelCoarsenerT(const MultilevelCoarsenerT&) = delete;
  MultilevelCoarsenerT(MultilevelCoarsenerT&&) = delete;
  MultilevelCoarsenerT & operator= (const MultilevelCoarsenerT &) = delete;
  MultilevelCoarsenerT & operator= (MultilevelCoarsenerT &&) = delete;

  ~MultilevelCoarsenerT() = default;

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void coarsenImpl() override {

  }

  bool uncoarsenImpl(std::unique_ptr<IRefiner>&) override {
    return true;
  }

  HyperGraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  UnionFind _uf;
  Rater _rater;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
using MultilevelCoarsener = MultilevelCoarsenerT<GlobalTypeTraits, ScorePolicy,
                                                 HeavyNodePenaltyPolicy, AcceptancePolicy>;
}  // namespace mt_kahypar
