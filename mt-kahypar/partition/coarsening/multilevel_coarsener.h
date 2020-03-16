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
#include "mt-kahypar/parallel/atomic_wrapper.h"
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
  using Rater = MultilevelVertexPairRater<TypeTraits,
                                          ScorePolicy,
                                          HeavyNodePenaltyPolicy,
                                          AcceptancePolicy>;
  using Rating = typename Rater::Rating;
  using Refiner = IRefinerT<TypeTraits>;

  enum class MatchingState : uint8_t {
    UNMATCHED = 0,
    MATCHING_IN_PROGRESS = 1,
    MATCHED = 2
  };

  #define STATE(X) static_cast<uint8_t>(X)
  using AtomicMatchingState = parallel::IntegralAtomicWrapper<uint8_t>;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsenerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    Base(hypergraph, context, task_group_id),
    _rater(hypergraph, context),
    _matching_state(),
    _cluster_weight(),
    _matching_partner(),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _num_matchings(0),
    _num_conflicts(0),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += hypergraph.numRemovedHypernodes();

    // Initialize internal data structures parallel
    tbb::parallel_invoke([&] {
      _matching_state.resize(hypergraph.initialNumNodes());
    }, [&] {
      _cluster_weight.resize(hypergraph.initialNumNodes());
    }, [&] {
      _matching_partner.resize(hypergraph.initialNumNodes());
    });

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

  ~MultilevelCoarsenerT() {
    parallel::parallel_free(_matching_state,
      _cluster_weight, _matching_partner);
  };

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

      parallel::scalable_vector<HypernodeID> cluster_ids(current_hg.initialNumNodes());
      current_hg.doParallelForAllNodes(_task_group_id, [&](const HypernodeID& hn) {
        const int node = common::get_numa_node_of_vertex(hn);
        const HypernodeID local_id = common::get_local_position_of_vertex(hn);
        ASSERT(local_id < current_vertices[node].size());
        current_vertices[node][local_id] = hn;
        // Reset clustering
        const HypernodeID original_id = current_hg.originalNodeID(hn);
        _matching_state[original_id] = STATE(MatchingState::UNMATCHED);
        _cluster_weight[original_id] = current_hg.nodeWeight(hn);
        _matching_partner[original_id] = original_id;
        cluster_ids[original_id] = original_id;
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
      const HypernodeID num_hns_before_pass = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
      const HypernodeID num_pins_before_pass = current_hg.initialNumPins();
      const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
      tbb::enumerable_thread_specific<HypernodeID> contracted_nodes(0);
      DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
      TBB::instance().execute_parallel_on_all_numa_nodes(_task_group_id, [&](const int node) {
        tbb::parallel_for(ID(0), current_hg.initialNumNodes(node), [&, node](const HypernodeID id) {
          ASSERT(id < current_vertices[node].size());
          const HypernodeID hn = current_vertices[node][id];
          const HypernodeID u = current_hg.originalNodeID(hn);
          // We perform rating if ...
          //  1.) The contraction limit of the current level is not reached
          //  2.) Vertex hn is not matched before
          if ( _matching_state[u] == STATE(MatchingState::UNMATCHED) ) {
            const HypernodeID current_num_nodes = num_hns_before_pass -
              contracted_nodes.combine(std::plus<HypernodeID>());
            if ( current_num_nodes > hierarchy_contraction_limit ) {
              ASSERT(current_hg.nodeIsEnabled(hn));
              const Rating rating = _rater.rate(current_hg, hn,
                cluster_ids, _cluster_weight, _max_allowed_node_weight);
              if ( rating.target != kInvalidHypernode ) {
                const HypernodeID v = current_hg.originalNodeID(rating.target);
                matchVertices(current_hg, u, v, cluster_ids,
                  contracted_nodes.local(), _num_conflicts.local());
                ++_num_matchings.local();
              }
            }
          }
        });
      });
      utils::Timer::instance().stop_timer("parallel_clustering");
      const HypernodeID current_num_nodes = num_hns_before_pass -
        contracted_nodes.combine(std::plus<HypernodeID>());
      DBG << V(current_num_nodes);

      HEAVY_COARSENING_ASSERT([&] {
        parallel::scalable_vector<HypernodeWeight> expected_weights(current_hg.initialNumNodes());
        // Verify that clustering is correct
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          const HypernodeID u = current_hg.originalNodeID(hn);
          const HypernodeID root_u = cluster_ids[u];
          if ( root_u != cluster_ids[root_u] ) {
            LOG << "Hypernode" << u << "is part of cluster" << root_u << ", but cluster"
                << root_u << "is also part of cluster" << cluster_ids[root_u];
            return false;
          }
          expected_weights[root_u] += current_hg.nodeWeight(hn);
        }

        // Verify that cluster weights are aggregated correct
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          const HypernodeID u = current_hg.originalNodeID(hn);
          const HypernodeID root_u = cluster_ids[u];
          if ( root_u == u && expected_weights[u] != _cluster_weight[u] ) {
            LOG << "The expected weight of cluster" << u << "is" << expected_weights[u]
                << ", but currently it is" << _cluster_weight[u];
            return false;
          }
        }
        return true;
      }(), "Parallel clustering computed invalid cluster ids and weights");

      const double reduction_vertices_percentage =
        static_cast<double>(num_hns_before_pass) /
        static_cast<double>(current_num_nodes);
      if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
        break;
      }
      _progress_bar += (num_hns_before_pass - current_num_nodes);

      utils::Timer::instance().start_timer("parallel_multilevel_contraction", "Parallel Multilevel Contraction");
      // Perform parallel contraction
      Base::performMultilevelContraction(std::move(cluster_ids));
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
      ++pass_nr;
    }
    _progress_bar += (initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();

    // Aggregate total number of matchings
    int64_t num_matchings = 0;
    for ( const size_t& c : _num_matchings ) {
      num_matchings += c;
    }
    // Aggregate total number of conflicts
    int64_t num_conflicts = 0;
    for ( const size_t& c : _num_conflicts ) {
      num_conflicts += c;
    }
    utils::Stats::instance().add_stat("coarsening_num_conflicts", num_matchings);
    utils::Stats::instance().add_stat("coarsening_num_matchings", num_matchings);
    DBG << V(num_conflicts) << V(num_matchings);

    Base::finalize();
  }

  /*!
   * We maintain the invariant during clustering that each cluster has a unique
   * representative and all vertices also part of that cluster point to that
   * representative. Let v be the representative of a cluster C_v, then for
   * all nodes u \in C_v follows that cluster_ids[u] = v.
   * If we perform sequential clustering, we can simply set
   * cluster_ids[u] = cluster_ids[v] to maintain our invariant. However,
   * things become more complicated if we perform parallel clustering.
   * Especially, if two neighbors u and v are concurrently matched, we have
   * to guarantee that our clustering fullfils our invariant. There are mainly
   * two different cases, which needs special attention:
   *   1.) u is matched with v and v is matched with u concurrently
   *   2.) u is matched with v and v is matched an other vertex w concurrently
   * The following functions guarantees that our invariant is fullfilled, if
   * vertices are matched concurrently.
   */
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool matchVertices(const HyperGraph& hypergraph,
                                                        const HypernodeID u,
                                                        const HypernodeID v,
                                                        parallel::scalable_vector<HypernodeID>& cluster_ids,
                                                        HypernodeID& contracted_nodes,
                                                        size_t& num_conflicts) {
    ASSERT(u < hypergraph.initialNumNodes());
    ASSERT(v < hypergraph.initialNumNodes());
    uint8_t unmatched = STATE(MatchingState::UNMATCHED);
    uint8_t match_in_progress = STATE(MatchingState::MATCHING_IN_PROGRESS);

    // Indicates that u wants to join the cluster of v.
    // Will be important later for conflict resolution.
    bool success = false;
    _matching_partner[u] = v;
    const HypernodeWeight weight_u = hypergraph.nodeWeight(hypergraph.globalNodeID(u));
    HypernodeWeight weight_v = _cluster_weight[v];
    if ( weight_u + weight_v <= _max_allowed_node_weight ) {

      if ( _matching_state[u].compare_and_exchange_strong(unmatched, match_in_progress) ) {
        // Current thread gets "ownership" for vertex u. Only threads with "ownership"
        // can change the cluster id of a vertex.

        uint8_t matching_state_v = _matching_state[v].load();
        if ( matching_state_v == STATE(MatchingState::MATCHED) ) {
          // Vertex v is already matched and will not change it cluster id any more.
          // In that case, it is safe to set the cluster id of u to the cluster id of v.
          if ( v == cluster_ids[v] ) {
            // In case v is also the representative of the cluster,
            // we change the cluster id of u to v, ...
            cluster_ids[u] = v;
            _cluster_weight[v] += weight_u;
            ++contracted_nodes;
            success = true;
          } else {
            // ... otherwise, we try again to match u with the
            // representative of the cluster.
            const HypernodeID cluster_v = cluster_ids[v];
            weight_v = _cluster_weight[cluster_v];
            if ( weight_u + weight_v <= _max_allowed_node_weight ) {
              ASSERT(_matching_state[cluster_v] == STATE(MatchingState::MATCHED));
              cluster_ids[u] = cluster_v;
              _cluster_weight[cluster_v] += weight_u;
              ++contracted_nodes;
              success = true;
            }
          }
        } else if ( _matching_state[v].compare_and_exchange_strong(unmatched, match_in_progress) ) {
          // Current thread has the "ownership" for u and v and can change the cluster id
          // of both vertices thread-safe.
          cluster_ids[u] = v;
          _cluster_weight[v] += weight_u;
          ++contracted_nodes;
          _matching_state[v] = STATE(MatchingState::MATCHED);
          success = true;
        } else {
          // State of v must be either MATCHING_IN_PROGRESS or an other thread changed the state
          // in the meantime to MATCHED. We have to wait until the state of v changed to
          // MATCHED or resolve the conflict if u is matched to v and v is matched to u and both
          // are in state MATCHING_IN_PROGRESS

          // Conflict Resolution
          while ( _matching_state[v] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {
            if ( _matching_partner[v] == u && u < v) {
              cluster_ids[u] = v;
              _cluster_weight[v] += weight_u;
              ++contracted_nodes;
              _matching_state[v] = STATE(MatchingState::MATCHED);
              success = true;
            }
          }

          // If u is still in state MATCHING_IN_PROGRESS its matching partner v
          // must be matched in the meantime with an other vertex. Therefore,
          // we try to match u with the representative v's cluster.
          if ( _matching_state[u] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {
            ASSERT( _matching_state[v] == STATE(MatchingState::MATCHED) );
            const HypernodeID cluster_v = cluster_ids[v];
            const HypernodeWeight weight_v = _cluster_weight[cluster_v];
            if ( weight_u + weight_v <= _max_allowed_node_weight ){
              cluster_ids[u] = cluster_v;
              _cluster_weight[cluster_v] += weight_u;
              ++contracted_nodes;
              success = true;
            }
          }
        }
        _rater.markAsMatched(u);
        _rater.markAsMatched(v);
        _matching_state[u] = STATE(MatchingState::MATCHED);
      } else {
        ++num_conflicts;
      }
    }
    return success;
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
  Rater _rater;
  parallel::scalable_vector<AtomicMatchingState> _matching_state;
  parallel::scalable_vector<AtomicWeight> _cluster_weight;
  parallel::scalable_vector<HypernodeID> _matching_partner;
  HypernodeWeight _max_allowed_node_weight;
  tbb::enumerable_thread_specific<size_t> _num_matchings;
  tbb::enumerable_thread_specific<size_t> _num_conflicts;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
using MultilevelCoarsener = MultilevelCoarsenerT<GlobalTypeTraits, ScorePolicy,
                                                 HeavyNodePenaltyPolicy, AcceptancePolicy>;
}  // namespace mt_kahypar
