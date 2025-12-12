/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <string>

#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

#include "kahypar-resources/meta/mandatory.h"

#include "include/mtkahypartypes.h"

#include "mt-kahypar/partition/coarsening/multilevel/clustering_context.h"
#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"
#include "mt-kahypar/partition/coarsening/multilevel/clustering_algorithms/single_round_lp.h"
#include "mt-kahypar/partition/coarsening/multilevel/clustering_algorithms/two_hop_clustering.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

template <class TypeTraits = Mandatory,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = NoWeightPenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class ThreePhaseCoarsener : public ICoarsener,
                            private MultilevelCoarsenerBase<TypeTraits> {
 private:
  using Base = MultilevelCoarsenerBase<TypeTraits>;
  using Rating = MultilevelVertexPairRater::Rating;
  using LPClustering = SingleRoundLP<ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy>;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using Base::_hg;
  using Base::_context;
  using Base::_timer;
  using Base::_uncoarseningData;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  ThreePhaseCoarsener(mt_kahypar_hypergraph_t hypergraph,
                      const Context& context,
                      uncoarsening_data_t* uncoarseningData) :
    Base(utils::cast<Hypergraph>(hypergraph),
         context,
         uncoarsening::to_reference<TypeTraits>(uncoarseningData)),
    _lp_clustering(_hg.initialNumNodes(), context),
    _two_hop_clustering(_hg.initialNumNodes(), context),
    _rater(_hg.initialNumNodes(), _hg.maxEdgeSize(), context),
    _clustering_data(_hg.initialNumNodes(), context),
    _initial_num_nodes(_hg.initialNumNodes()),
    _num_communities(kInvalidPartition),
    _current_vertices(),
    _pass_nr(0),
    _progress_bar(_hg.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += _hg.numRemovedHypernodes();
    _current_vertices.resize(_hg.initialNumNodes());
  }

  ThreePhaseCoarsener(const ThreePhaseCoarsener&) = delete;
  ThreePhaseCoarsener(ThreePhaseCoarsener&&) = delete;
  ThreePhaseCoarsener & operator= (const ThreePhaseCoarsener &) = delete;
  ThreePhaseCoarsener & operator= (ThreePhaseCoarsener &&) = delete;

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void initializeImpl() override {
    if ( _context.partition.enable_logging && _context.partition.enable_progress_bar ) {
      _progress_bar.enable();
    }
  }

  bool shouldNotTerminateImpl() const override {
    return Base::currentNumNodes() > _context.coarsening.contraction_limit;
  }

  bool coarseningPassImpl() override {
    HighResClockTimepoint round_start = std::chrono::high_resolution_clock::now();
    Hypergraph& current_hg = Base::currentHypergraph();
    DBG << V(_pass_nr)
        << V(current_hg.initialNumNodes())
        << V(current_hg.initialNumEdges())
        << V(current_hg.initialNumPins());

    // Random shuffle vertices of current hypergraph
    _current_vertices.resize(current_hg.initialNumNodes());
    parallel::scalable_vector<HypernodeID> cluster_ids;
    tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID hn) {
      _current_vertices[hn] = hn;
    });

    // initialization of various things
    const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
    ClusteringContext<Hypergraph> cc(_context, hierarchy_contraction_limit, cluster_ids,
                                     _rater, _clustering_data);
    cc.may_ignore_communities = shouldIgnoreCommunities(hierarchy_contraction_limit);
    cc.initializeCoarseningPass(current_hg, _context);

    // TODO: degree zero nodes?!
    // Step 1: standard label propagation clustering
    HypernodeID current_num_nodes = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
    coarseningRound("lp_clustering", "LP Clustering", current_hg, _lp_clustering, cc);
    _progress_bar += (current_num_nodes - cc.finalNumNodes());
    current_num_nodes = cc.currentNumNodes();

    // Step 2: if the size didn't shrink far enough, use two-hop clustering for low degree nodes
    if (current_num_nodes > hierarchy_contraction_limit) {
      DBG << "Start Two-Hop Coarsening: " << V(current_num_nodes) << V(hierarchy_contraction_limit);
      coarseningRound("two_hop_clustering", "Two-Hop Clustering", current_hg, _two_hop_clustering, cc);
      _progress_bar += (current_num_nodes - cc.finalNumNodes());
      current_num_nodes = cc.currentNumNodes();
    }

    if (current_num_nodes > hierarchy_contraction_limit) {
      // If the size is still too large, the reason could be that there are too many communities.
      // We initialize the community count, so the next round can decide to ignore communities
      // (delayed initialization since it is not completely free)
      initializeCommunityCount(current_hg);
    }

    DBG << V(current_num_nodes) << V(hierarchy_contraction_limit);
    bool should_continue = cc.finalize(current_hg, _context);
    if (!should_continue) {
      return false;
    }

    _timer.start_timer("contraction", "Contraction");
    // Perform parallel contraction
    _uncoarseningData.performMultilevelContraction(std::move(cluster_ids), false /* deterministic */, round_start);
    _timer.stop_timer("contraction");

    ++_pass_nr;
    return true;
  }

  template<typename ClusteringAlgo>
  HypernodeID coarseningRound(const char* timer_key,
                              const char* timer_name,
                              const Hypergraph& current_hg,
                              ClusteringAlgo& algo,
                              ClusteringContext<Hypergraph>& cc) {
    _timer.start_timer(timer_key, timer_name);
    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.start_timer(timer_key + std::string("_level_") + std::to_string(_pass_nr), "Level " + std::to_string(_pass_nr));
    }

    if ( _enable_randomization ) {
      utils::Randomize::instance().parallelShuffleVector( _current_vertices, UL(0), _current_vertices.size());
    }
    algo.performClustering(current_hg, _current_vertices, cc, current_hg.hasFixedVertices());

    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.stop_timer(timer_key + std::string("_level_") + std::to_string(_pass_nr));
    }
    _timer.stop_timer(timer_key);
    return cc.finalNumNodes();
  }

  void terminateImpl() override {
    _progress_bar += (_initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();
    _uncoarseningData.finalizeCoarsening();
  }

  HypernodeID currentNumberOfNodesImpl() const override {
    return Base::currentNumNodes();
  }

  mt_kahypar_hypergraph_t coarsestHypergraphImpl() override {
    return mt_kahypar_hypergraph_t {
      reinterpret_cast<mt_kahypar_hypergraph_s*>(
        &Base::currentHypergraph()), Hypergraph::TYPE };
  }

  mt_kahypar_partitioned_hypergraph_t coarsestPartitionedHypergraphImpl() override {
    return mt_kahypar_partitioned_hypergraph_t {
      reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(
        &Base::currentPartitionedHypergraph()), PartitionedHypergraph::TYPE };
  }

  HypernodeID hierarchyContractionLimit(const Hypergraph& hypergraph) const {
    return std::max( static_cast<HypernodeID>( static_cast<double>(hypergraph.initialNumNodes() -
      hypergraph.numRemovedHypernodes()) / _context.coarsening.maximum_shrink_factor ),
      _context.coarsening.contraction_limit );
  }

  void initializeCommunityCount(const Hypergraph& hypergraph) {
    if (_num_communities == kInvalidPartition) {
      _num_communities = utils::communityCount(hypergraph);
    }
  }

  bool shouldIgnoreCommunities(HypernodeID hierarchy_contraction_limit) {
    return _num_communities != kInvalidPartition && UL(_num_communities) > hierarchy_contraction_limit;
  }

  LPClustering _lp_clustering;
  TwoHopClustering _two_hop_clustering;
  MultilevelVertexPairRater _rater;
  ConcurrentClusteringData _clustering_data;
  HypernodeID _initial_num_nodes;
  PartitionID _num_communities;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  int _pass_nr;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
