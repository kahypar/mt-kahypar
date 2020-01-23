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


  static constexpr bool debug = true;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsenerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _uf(hypergraph),
    _rater(hypergraph, context, _uf),
    // _progress_bar(hypergraph.initialNumNodes(), 0, false),
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
    int pass_nr = 0;
    parallel::scalable_vector<HypernodeID> current_vertices;
    while ( _hg.initialNumNodes() > _context.coarsening.contraction_limit ) {
      DBG << V(pass_nr) << V(_hg.initialNumNodes()) << V(_hg.initialNumEdges());

      // Random shuffle vertices of current hypergraph
      current_vertices.clear();
      for ( const HypernodeID& hn : _hg.nodes() ) {
        current_vertices.emplace_back(hn);
      }
      utils::Randomize::instance().parallelShuffleVector(current_vertices);

      // We iterate in parallel over all vertices of the hypergraph and compute its contraction
      // partner. The vertices are processed on the numa node which they are placed on.
      // Matched vertices are linked in a concurrent union find data structure, that also aggregates
      // weights of the resulting clusters and keep track of the number of nodes left, if we would
      // contract all matched vertices.
      _rater.resetMatches();
      const HypernodeID num_hns_before_pass = _uf.numDistinctSets();
      const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit();
      DBG << V(_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
        tbb::parallel_for(0UL, _hg.initialNumNodes(node), [&, node](const HypernodeID id) {
          const HypernodeID hn = StreamingHyperGraph::get_global_node_id(node, id);
          if ( _uf.numDistinctSets() > hierarchy_contraction_limit &&
               _hg.nodeIsEnabled(hn) &&
               !_hg.isHighDegreeVertex(hn) ) {
            const Rating rating = _rater.rate(hn);
            if ( rating.target != kInvalidHypernode ) {
              _rater.markAsMatched(hn);
              _rater.markAsMatched(rating.target);
              _uf.link(_hg.originalNodeID(hn), _hg.originalNodeID(rating.target));
            }
          }
        });
      });
      DBG << V(_uf.numDistinctSets());

      if ( num_hns_before_pass == _uf.numDistinctSets() ) {
        break;
      }

      // Compute community structure that is given by the representatives of each
      // node in the union find data structure.
      parallel::scalable_vector<HypernodeID> communities(_hg.initialNumNodes());
      tbb::parallel_for(0UL, _hg.initialNumNodes(), [&](const HypernodeID id) {
        communities[id] = _uf.find(id);
      });

      // TODO(heuer): Contract Communities here

      break;
      _uf.reset(_hg);
      ++pass_nr;
    }
  }

  bool uncoarsenImpl(std::unique_ptr<IRefiner>&) override {
    return true;
  }

  HypernodeID hierarchyContractionLimit() const {
    return std::max( static_cast<HypernodeID>( static_cast<double>(_hg.initialNumNodes()) /
      _context.coarsening.multilevel_shrink_factor ), _context.coarsening.contraction_limit );
  }

  HyperGraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  UnionFind _uf;
  Rater _rater;
  // utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
using MultilevelCoarsener = MultilevelCoarsenerT<GlobalTypeTraits, ScorePolicy,
                                                 HeavyNodePenaltyPolicy, AcceptancePolicy>;
}  // namespace mt_kahypar
