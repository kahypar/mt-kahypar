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
#include "tbb/parallel_reduce.h"

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
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class MultilevelCoarsener : public ICoarsener,
                            private MultilevelCoarsenerBase {
 private:

  using Base = MultilevelCoarsenerBase;
  using Rater = MultilevelVertexPairRater<ScorePolicy,
                                          HeavyNodePenaltyPolicy,
                                          AcceptancePolicy>;
  using Rating = typename Rater::Rating;

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
  MultilevelCoarsener(Hypergraph& hypergraph,
                      const Context& context,
                      const TaskGroupID task_group_id,
                      const bool top_level) :
    Base(hypergraph, context, task_group_id, top_level),
    _rater(hypergraph, context),
    _current_vertices(),
    _matching_state(),
    _cluster_weight(),
    _matching_partner(),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += hypergraph.numRemovedHypernodes();

    // Initialize internal data structures parallel
    tbb::parallel_invoke([&] {
      _current_vertices.resize(hypergraph.initialNumNodes());
    }, [&] {
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
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
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

  MultilevelCoarsener(const MultilevelCoarsener&) = delete;
  MultilevelCoarsener(MultilevelCoarsener&&) = delete;
  MultilevelCoarsener & operator= (const MultilevelCoarsener &) = delete;
  MultilevelCoarsener & operator= (MultilevelCoarsener &&) = delete;

  ~MultilevelCoarsener() {
    parallel::parallel_free(
      _current_vertices, _matching_state,
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
    while ( Base::currentNumNodes() > _context.coarsening.contraction_limit ) {
      HighResClockTimepoint round_start = std::chrono::high_resolution_clock::now();
      Hypergraph& current_hg = Base::currentHypergraph();
      DBG << V(pass_nr)
          << V(current_hg.initialNumNodes())
          << V(current_hg.initialNumEdges())
          << V(current_hg.initialNumPins());

      // Random shuffle vertices of current hypergraph
      utils::Timer::instance().start_timer("shuffle_vertices", "Shuffle Vertices");
      _current_vertices.resize(current_hg.initialNumNodes());
      parallel::scalable_vector<HypernodeID> cluster_ids(current_hg.initialNumNodes());
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID hn) {
        ASSERT(hn < _current_vertices.size());
        // Reset clustering
        _current_vertices[hn] = hn;
        _matching_state[hn] = STATE(MatchingState::UNMATCHED);
        _matching_partner[hn] = hn;
        cluster_ids[hn] = hn;
        if ( current_hg.nodeIsEnabled(hn) ) {
          _cluster_weight[hn] = current_hg.nodeWeight(hn);
        }
      });

      if ( _enable_randomization ) {
        utils::Randomize::instance().parallelShuffleVector(
          _current_vertices, 0UL, _current_vertices.size());
      }
      utils::Timer::instance().stop_timer("shuffle_vertices");

      // We iterate in parallel over all vertices of the hypergraph and compute its contraction partner.
      // Matched vertices are linked in a concurrent union find data structure, that also aggregates
      // weights of the resulting clusters and keep track of the number of nodes left, if we would
      // contract all matched vertices.
      utils::Timer::instance().start_timer("parallel_clustering", "Parallel Clustering");
      if ( _context.partition.show_detailed_clustering_timings ) {
        utils::Timer::instance().start_timer("clustering_level_" + std::to_string(pass_nr), "Level " + std::to_string(pass_nr));
      }
      _rater.resetMatches();
      _rater.setCurrentNumberOfNodes(current_hg.initialNumNodes());
      const HypernodeID num_hns_before_pass = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
      const HypernodeID num_pins_before_pass = current_hg.initialNumPins();
      const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
      DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
      HypernodeID current_num_nodes = num_hns_before_pass;
      tbb::enumerable_thread_specific<HypernodeID> contracted_nodes(0);
      tbb::enumerable_thread_specific<HypernodeID> num_nodes_update_threshold(0);
      tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID id) {
        ASSERT(id < _current_vertices.size());
        const HypernodeID hn = _current_vertices[id];
        if ( current_hg.nodeIsEnabled(hn) ) {
          // We perform rating if ...
          //  1.) The contraction limit of the current level is not reached
          //  2.) Vertex hn is not matched before
          const HypernodeID u = hn;
          if ( _matching_state[u] == STATE(MatchingState::UNMATCHED) ) {
            if ( current_num_nodes > hierarchy_contraction_limit ) {
              ASSERT(current_hg.nodeIsEnabled(hn));
              const Rating rating = _rater.rate(current_hg, hn,
                cluster_ids, _cluster_weight, _max_allowed_node_weight);
              if ( rating.target != kInvalidHypernode ) {
                const HypernodeID v = rating.target;
                HypernodeID& local_contracted_nodes = contracted_nodes.local();
                matchVertices(current_hg, u, v, cluster_ids, local_contracted_nodes);

                // To maintain the current number of nodes of the hypergraph each PE sums up
                // its number of contracted nodes locally. To compute the current number of
                // nodes, we have to sum up the number of contracted nodes of each PE. This
                // operation becomes more expensive the more PEs are participating in coarsening.
                // In order to prevent expensive updates of the current number of nodes, we
                // define a threshold which the local number of contracted nodes have to exceed
                // before the current PE updates the current number of nodes. This threshold is defined
                // by the distance to the current contraction limit divided by the number of PEs.
                // Once one PE exceeds this bound the first time it is not possible that the
                // contraction limit is reached, because otherwise an other PE would update
                // the global current number of nodes before. After update the threshold is
                // increased by the new difference (in number of nodes) to the contraction limit
                // divided by the number of PEs.
                if (  local_contracted_nodes >= num_nodes_update_threshold.local() ) {
                  current_num_nodes = num_hns_before_pass -
                    contracted_nodes.combine(std::plus<HypernodeID>());
                  num_nodes_update_threshold.local() +=
                    (current_num_nodes - hierarchy_contraction_limit) /
                    _context.shared_memory.num_threads;
                }
              }
            }
          }
        }
      });
      if ( _context.partition.show_detailed_clustering_timings ) {
        utils::Timer::instance().stop_timer("clustering_level_" + std::to_string(pass_nr));
      }
      utils::Timer::instance().stop_timer("parallel_clustering");
      current_num_nodes = num_hns_before_pass -
        contracted_nodes.combine(std::plus<HypernodeID>());
      DBG << V(current_num_nodes);

      HEAVY_COARSENING_ASSERT([&] {
        parallel::scalable_vector<HypernodeWeight> expected_weights(current_hg.initialNumNodes());
        // Verify that clustering is correct
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          const HypernodeID u = hn;
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
          const HypernodeID u = hn;
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
      Base::performMultilevelContraction(std::move(cluster_ids), round_start);
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
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool matchVertices(const Hypergraph& hypergraph,
                                                        const HypernodeID u,
                                                        const HypernodeID v,
                                                        parallel::scalable_vector<HypernodeID>& cluster_ids,
                                                        HypernodeID& contracted_nodes) {
    ASSERT(u < hypergraph.initialNumNodes());
    ASSERT(v < hypergraph.initialNumNodes());
    uint8_t unmatched = STATE(MatchingState::UNMATCHED);
    uint8_t match_in_progress = STATE(MatchingState::MATCHING_IN_PROGRESS);

    // Indicates that u wants to join the cluster of v.
    // Will be important later for conflict resolution.
    bool success = false;
    const HypernodeWeight weight_u = hypergraph.nodeWeight(u);
    HypernodeWeight weight_v = _cluster_weight[v];
    if ( weight_u + weight_v <= _max_allowed_node_weight ) {

      if ( _matching_state[u].compare_exchange_strong(unmatched, match_in_progress) ) {
        _matching_partner[u] = v;
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
        } else if ( _matching_state[v].compare_exchange_strong(unmatched, match_in_progress) ) {
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
          // MATCHED or resolve the conflict if u is matched within a cyclic matching dependency

          // Conflict Resolution
          while ( _matching_state[v] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {

            // Check if current vertex is in a cyclic matching dependency
            HypernodeID cur_u = u;
            HypernodeID smallest_node_id_in_cycle = cur_u;
            while ( _matching_partner[cur_u] != u && _matching_partner[cur_u] != cur_u ) {
              cur_u = _matching_partner[cur_u];
              smallest_node_id_in_cycle = std::min(smallest_node_id_in_cycle, cur_u);
            }

            // Resolve cyclic matching dependency
            // Vertex with smallest id starts to resolve conflict
            const bool is_in_cyclic_dependency = _matching_partner[cur_u] == u;
            if ( is_in_cyclic_dependency && u == smallest_node_id_in_cycle) {
              cluster_ids[u] = v;
              _cluster_weight[v] += weight_u;
              ++contracted_nodes;
              _matching_state[v] = STATE(MatchingState::MATCHED);
              _matching_state[u] = STATE(MatchingState::MATCHED);
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
        _matching_partner[u] = u;
        _matching_state[u] = STATE(MatchingState::MATCHED);
      }
    }
    return success;
  }

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) override {
    return Base::doUncoarsen(label_propagation, fm);
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

  HypernodeID hierarchyContractionLimit(const Hypergraph& hypergraph) const {
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
  parallel::scalable_vector<HypernodeID> _current_vertices;
  parallel::scalable_vector<AtomicMatchingState> _matching_state;
  parallel::scalable_vector<AtomicWeight> _cluster_weight;
  parallel::scalable_vector<HypernodeID> _matching_partner;
  HypernodeWeight _max_allowed_node_weight;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
