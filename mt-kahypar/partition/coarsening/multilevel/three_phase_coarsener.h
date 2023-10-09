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

#include "include/libmtkahypartypes.h"

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
#include "mt-kahypar/partition/coarsening/policies/rating_degree_similarity_policy.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

template <class TypeTraits = Mandatory,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = NoWeightPenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched,
          class SimilarityPolicy = PreserveRebalancingNodesPolicy>
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
    _similarity_policy(_hg.initialNumNodes()),
    _always_accept_policy(_hg.initialNumNodes()),
    _rater(_hg.initialNumNodes(), _hg.maxEdgeSize(), context),
    _clustering_data(_hg.initialNumNodes(), context),
    _num_nodes_tracker(),
    _initial_num_nodes(_hg.initialNumNodes()),
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
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar ) {
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

    const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
    const HypernodeID target_contraction_size = targetContractionSize(current_hg);
    ASSERT(target_contraction_size >= hierarchy_contraction_limit);
    const HypernodeID num_hns_before_pass = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
    HypernodeID current_num_nodes = num_hns_before_pass;

    // initialization of various things
    ClusteringContext<Hypergraph> cc(_context, hierarchy_contraction_limit, cluster_ids,
                                     _rater, _clustering_data, _num_nodes_tracker);
    cc.initializeCoarseningPass(current_hg, _context);
    _timer.start_timer("init_similarity", "Initialize Similarity Data");
    _similarity_policy.initialize(current_hg, _context);
    _timer.stop_timer("init_similarity");
    _always_accept_policy.initialize(current_hg, _context);

    // TODO: degree zero nodes?!
    // Phase 1: LP coarsening, but forbid contraction of low degree nodes onto high degree nodes
    coarseningRound("first_lp_round", "First LP round",
                    current_hg, _lp_clustering, _similarity_policy, cc);
    _progress_bar += (current_num_nodes - _num_nodes_tracker.finalNumNodes());
    current_num_nodes = _num_nodes_tracker.currentNumNodes();

    // Phase 2: Two-hop coarsening for low degree nodes
    if (current_num_nodes > (_context.coarsening.delayed_two_hop_coarsening ?
          target_contraction_size : hierarchy_contraction_limit)) {
      DBG << "Start Two-Hop Coarsening: " << V(_num_nodes_tracker.currentNumNodes()) << V(hierarchy_contraction_limit);
      coarseningRound("first_two_hop_round", "First two-hop round",
                     current_hg, _two_hop_clustering, _similarity_policy, cc);
      _progress_bar += (current_num_nodes - _num_nodes_tracker.finalNumNodes());
      current_num_nodes = _num_nodes_tracker.currentNumNodes();
    }

    // Phase 3: LP and two-hop coarsening with all contractions allowed (as well as contracting size 1 communities)
    cc.may_ignore_communities = true;
    cc.contract_aggressively = true;
    cc.hierarchy_contraction_limit = target_contraction_size;
    if (current_num_nodes > target_contraction_size) {
      DBG << "Start Second LP round: " << V(_num_nodes_tracker.currentNumNodes()) << V(target_contraction_size);
      coarseningRound("second_lp_round", "Second LP round",
                      current_hg, _lp_clustering, _always_accept_policy, cc);
      _progress_bar += (current_num_nodes - _num_nodes_tracker.finalNumNodes());
      current_num_nodes = _num_nodes_tracker.currentNumNodes();
    }
    if (current_num_nodes > target_contraction_size) {
      DBG << "Start Second Two-Hop Coarsening: " << V(_num_nodes_tracker.currentNumNodes()) << V(target_contraction_size);
      coarseningRound("second_two_hop_round", "Second two-hop round",
                     current_hg, _two_hop_clustering, _always_accept_policy, cc);
      _progress_bar += (current_num_nodes - _num_nodes_tracker.finalNumNodes());
      current_num_nodes = _num_nodes_tracker.currentNumNodes();
    }

    DBG << V(current_num_nodes) << V(target_contraction_size) << V(hierarchy_contraction_limit);
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

  template<typename ClusteringAlgo, typename CurrentSimilarityPolicy>
  HypernodeID coarseningRound(const char* timer_key, const char* timer_name,
                              const Hypergraph& current_hg, ClusteringAlgo& algo,
                              const CurrentSimilarityPolicy& similarity, ClusteringContext<Hypergraph>& cc) {
    _timer.start_timer(timer_key, timer_name);
    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.start_timer(timer_key + std::string("_level_") + std::to_string(_pass_nr), "Level " + std::to_string(_pass_nr));
    }
    if ( _enable_randomization ) {
      utils::Randomize::instance().parallelShuffleVector( _current_vertices, UL(0), _current_vertices.size());
    }

    auto weight_ratio_fn = [&](const HypernodeID hn) {
      return _similarity_policy.weightRatioForNode(current_hg, hn);
    };
    if ( current_hg.hasFixedVertices() ) {
      algo.template performClustering<true>(current_hg, _current_vertices, similarity, cc, weight_ratio_fn);
    } else {
      algo.template performClustering<false>(current_hg, _current_vertices, similarity, cc, weight_ratio_fn);
    }

    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.stop_timer(timer_key + std::string("_level_") + std::to_string(_pass_nr));
    }
    _timer.stop_timer(timer_key);
    return _num_nodes_tracker.finalNumNodes();
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

  HypernodeID targetContractionSize(const Hypergraph& hypergraph) const {
    return std::max( static_cast<HypernodeID>( static_cast<double>(hypergraph.initialNumNodes() -
      hypergraph.numRemovedHypernodes()) / _context.coarsening.min_accepted_shrink_factor ),
      _context.coarsening.contraction_limit );
  }

  LPClustering _lp_clustering;
  TwoHopClustering _two_hop_clustering;
  SimilarityPolicy _similarity_policy;
  AlwaysAcceptPolicy _always_accept_policy;
  MultilevelVertexPairRater _rater;
  ConcurrentClusteringData _clustering_data;
  NumNodesTracker _num_nodes_tracker;
  HypernodeID _initial_num_nodes;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  int _pass_nr;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
