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
#include "tbb/parallel_for.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel_coarsener_base.h"
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
class MultilevelCoarsenerT : public ICoarsenerT<TypeTraits>,
                             private MultilevelCoarsenerBase<TypeTraits> {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  using Base = MultilevelCoarsenerBase<TypeTraits>;
  using UnionFind = ds::ConcurrentUnionFind<HyperGraph>;
  using Rater = MultilevelVertexPairRater<TypeTraits,
                                          ScorePolicy,
                                          HeavyNodePenaltyPolicy,
                                          AcceptancePolicy>;
  using Rating = typename Rater::Rating;

  using Refiner = IRefinerT<TypeTraits>;

  static constexpr bool debug = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsenerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    Base(hypergraph, context, task_group_id),
    _uf(hypergraph),
    _rater(hypergraph, context, _uf),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += hypergraph.numRemovedHypernodes();
    if ( _context.coarsening.use_adaptive_max_allowed_node_weight &&
          hypergraph.totalWeight() !=
          static_cast<HypernodeWeight>(hypergraph.initialNumNodes()) ) {
      // If we have a weighted instance and adaptive maximum node weight is
      // enabled we adapt the maximum allowed node such that it is greater
      // than the heaviest node of the hypergraph.
      const HypernodeWeight max_vertex_weight = tbb::parallel_reduce(
      tbb::blocked_range<HypernodeID>(ID(0), hypergraph.initialNumNodes()), 0,
      [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID id = range.begin(); id < range.end(); ++id) {
          const HypernodeID hn = hypergraph.globalNodeID(id);
          if ( hypergraph.nodeIsEnabled(hn) ) {
            weight = std::max(weight, hypergraph.nodeWeight(hn));
          }
        }
        return weight;
      }, [](const HypernodeWeight lhs, const HypernodeWeight rhs) {
        return std::max(lhs, rhs);
      });
      double node_weight_multiplier = std::pow(2.0, std::ceil(std::log2(
        static_cast<double>(max_vertex_weight) / static_cast<double>(_max_allowed_node_weight))));
      if ( node_weight_multiplier > 1.0 ) {
        _max_allowed_node_weight = increaseMaximumAllowedNodeWeight(node_weight_multiplier);
      }
    }
  }

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
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar ) {
      _progress_bar.enable();
    }

    int pass_nr = 0;
    const HypernodeID initial_num_nodes = Base::currentNumNodes();
    parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> current_vertices(TBB::instance().num_used_numa_nodes());
    while ( Base::currentNumNodes() > _context.coarsening.contraction_limit ) {
      HyperGraph& current_hg = Base::currentHypergraph();
      DBG << V(pass_nr)
          << V(current_hg.initialNumNodes())
          << V(current_hg.initialNumEdges())
          << V(current_hg.initialNumPins());

      // Random shuffle vertices of current hypergraph
      utils::Timer::instance().start_timer("shuffle_vertices", "Shuffle Vertices");
      for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
        current_vertices[node].resize(current_hg.initialNumNodes(node));
      }
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID id) {
        const HypernodeID hn = current_hg.globalNodeID(id);
        const int node = common::get_numa_node_of_vertex(hn);
        const HypernodeID local_id = common::get_local_position_of_vertex(hn);
        ASSERT(local_id < current_vertices[node].size());
        current_vertices[node][local_id] = hn;
      });

      if ( _enable_randomization ) {
        for ( int node = 0; node < TBB::instance().num_used_numa_nodes(); ++node ) {
          utils::Randomize::instance().localizedParallelShuffleVector(
            current_vertices[node], 0UL, current_vertices[node].size(), _context.shared_memory.shuffle_block_size);
        }
      }
      utils::Timer::instance().stop_timer("shuffle_vertices");

      // We iterate in parallel over all vertices of the hypergraph and compute its contraction
      // partner. The vertices are processed on the numa node which they are placed on.
      // Matched vertices are linked in a concurrent union find data structure, that also aggregates
      // weights of the resulting clusters and keep track of the number of nodes left, if we would
      // contract all matched vertices.
      utils::Timer::instance().start_timer("parallel_clustering", "Parallel Clustering");
      _rater.resetMatches();
      const HypernodeID num_hns_before_pass = current_hg.initialNumNodes();
      const HypernodeID num_pins_before_pass = current_hg.initialNumPins();
      const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
      DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
        tbb::parallel_for(ID(0), current_hg.initialNumNodes(node), [&, node](const HypernodeID id) {
          ASSERT(id < current_vertices[node].size());
          const HypernodeID hn = current_vertices[node][id];
          // We perform rating if ...
          //  1.) The contraction limit of the current level is not reached
          //  2.) Vertex hn is enabled
          //  3.) Vertex hn is not a high degree vertex
          //  4.) Vertex hn is not matched before
          if ( _uf.numDistinctSets() > hierarchy_contraction_limit &&
               current_hg.nodeIsEnabled(hn) &&
               ( !_enable_randomization || !_rater.isMatched(current_hg, hn) ) ) {
            const Rating rating = _rater.rate(current_hg, hn, _max_allowed_node_weight);
            if ( rating.target != kInvalidHypernode ) {
              // Check that if we contract both vertices, that the resulting weight
              // is less than the maximum allowed node weight.
              const HypernodeID original_hn_id = current_hg.originalNodeID(hn);
              const HypernodeID original_target_id = current_hg.originalNodeID(rating.target);
              const bool is_same_set = _uf.isSameSet(original_hn_id, original_target_id);
              const HypernodeWeight contracted_weight = is_same_set ? _uf.weight(original_hn_id) :
                _uf.weight(original_hn_id) + _uf.weight(original_target_id);
              if ( contracted_weight <= _max_allowed_node_weight &&
                   _uf.numDistinctSets() > hierarchy_contraction_limit ) {
                _uf.link(original_hn_id, original_target_id);
                _rater.markAsMatched(current_hg, hn);
                _rater.markAsMatched(current_hg, rating.target);
              }
            }
          }
        });
      });
      utils::Timer::instance().stop_timer("parallel_clustering");
      DBG << V(num_hns_before_pass) << V(_uf.numDistinctSets());

      const double reduction_vertices_percentage =
        static_cast<double>(num_hns_before_pass) /
        static_cast<double>(_uf.numDistinctSets());
      if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
        break;
      }
      _progress_bar += (num_hns_before_pass - _uf.numDistinctSets());

      // Compute community structure that is given by the representatives of each
      // node in the union find data structure.
      utils::Timer::instance().start_timer("parallel_multilevel_contraction", "Parallel Multilevel Contraction");
      parallel::scalable_vector<HypernodeID> communities(current_hg.initialNumNodes());
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID id) {
        communities[id] = _uf.find(id);
      });
      Base::performMultilevelContraction(std::move(communities));
      utils::Timer::instance().stop_timer("parallel_multilevel_contraction");

      if ( _context.coarsening.use_adaptive_max_allowed_node_weight ) {
        // If the reduction ratio of the number of vertices or pins is below
        // a certain threshold, we increase the maximum allowed node weight by
        // a factor of two. Idea behind this is that if we are not able to reduce
        // the number of nodes or pins by a significant ratio, then some vertices
        // reach their maximum allowed node weight and are not able to contract
        // with other nodes, which prevents some high score contractions.
        const double reduction_pins_percentage =
          static_cast<double>(num_pins_before_pass) /
          static_cast<double>(Base::currentHypergraph().initialNumPins());
        const bool reduction_vertices_below_threshold = reduction_vertices_percentage <
          _context.coarsening.adaptive_node_weight_shrink_factor_threshold;
        const bool reduction_pins_below_threshold = reduction_pins_percentage <
          _context.coarsening.adaptive_node_weight_shrink_factor_threshold;
        if ( ( reduction_vertices_below_threshold && reduction_pins_below_threshold ) ||
             ( !reduction_vertices_below_threshold && reduction_pins_below_threshold ) ) {
          _max_allowed_node_weight = increaseMaximumAllowedNodeWeight(2.0);
        }
        DBG << V(reduction_vertices_percentage)
            << V(reduction_pins_percentage)
            << V(_max_allowed_node_weight);
      }

      _uf.reset(Base::currentHypergraph());
      ++pass_nr;
    }
    _progress_bar += (initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();
    Base::finalize();
  }

  PartitionedHyperGraph&& uncoarsenImpl(std::unique_ptr<Refiner>& label_propagation) override {
    return Base::doUncoarsen(label_propagation);
  }

  HyperGraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHyperGraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

  HypernodeID hierarchyContractionLimit(const HyperGraph& hypergraph) const {
    return std::max( static_cast<HypernodeID>( static_cast<double>(hypergraph.initialNumNodes() -
      hypergraph.numRemovedHypernodes()) / _context.coarsening.maximum_shrink_factor ),
      _context.coarsening.contraction_limit );
  }

  HypernodeWeight increaseMaximumAllowedNodeWeight(const double multiplier) {
    HypernodeWeight max_part_weight = 0;
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      max_part_weight = std::max(max_part_weight,
        _context.partition.max_part_weights[block]);
    }
    return std::min( multiplier * static_cast<double>(_max_allowed_node_weight),
      std::max( max_part_weight / _context.coarsening.max_allowed_weight_fraction,
        static_cast<double>(_context.coarsening.max_allowed_node_weight ) ) );
  }

  using Base::_context;
  using Base::_task_group_id;
  UnionFind _uf;
  Rater _rater;
  HypernodeWeight _max_allowed_node_weight;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
using MultilevelCoarsener = MultilevelCoarsenerT<GlobalTypeTraits, ScorePolicy,
                                                 HeavyNodePenaltyPolicy, AcceptancePolicy>;
}  // namespace mt_kahypar
