/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

template <class TypeTraits = Mandatory,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = NoWeightPenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class MultilevelCoarsener : public ICoarsener,
                            private MultilevelCoarsenerBase<TypeTraits> {
 private:

  using Base = MultilevelCoarsenerBase<TypeTraits>;
  using Rating = MultilevelVertexPairRater::Rating;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  MultilevelCoarsener(mt_kahypar_hypergraph_t hypergraph,
                      const Context& context,
                      uncoarsening_data_t* uncoarseningData) :
    Base(utils::cast<Hypergraph>(hypergraph),
         context,
         uncoarsening::to_reference<TypeTraits>(uncoarseningData)),
    _rater(utils::cast<Hypergraph>(hypergraph).initialNumNodes(),
           utils::cast<Hypergraph>(hypergraph).maxEdgeSize(), context),
    _clustering_data(_hg.initialNumNodes(), _context),
    _initial_num_nodes(utils::cast<Hypergraph>(hypergraph).initialNumNodes()),
    _current_vertices(),
    _pass_nr(0),
    _progress_bar(utils::cast<Hypergraph>(hypergraph).initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += _hg.numRemovedHypernodes();
    _current_vertices.resize(_hg.initialNumNodes());
  }

  MultilevelCoarsener(const MultilevelCoarsener&) = delete;
  MultilevelCoarsener(MultilevelCoarsener&&) = delete;
  MultilevelCoarsener & operator= (const MultilevelCoarsener &) = delete;
  MultilevelCoarsener & operator= (MultilevelCoarsener &&) = delete;

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
    _clustering_data.initializeCoarseningPass(current_hg, cluster_ids);

    if ( _enable_randomization ) {
      utils::Randomize::instance().parallelShuffleVector( _current_vertices, UL(0), _current_vertices.size());
    }

    const HypernodeID num_hns_before_pass =
      current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
    HypernodeID current_num_nodes = 0;
    if ( current_hg.hasFixedVertices() ) {
      current_num_nodes = performClustering<true>(current_hg, cluster_ids);
    } else {
      current_num_nodes = performClustering<false>(current_hg, cluster_ids);
    }
    DBG << V(current_num_nodes);

    HEAVY_COARSENING_ASSERT(_clustering_data.verifyClustering(current_hg, cluster_ids),
                            "Parallel clustering computed invalid cluster ids and weights");

    const double reduction_vertices_percentage =
      static_cast<double>(num_hns_before_pass) /
      static_cast<double>(current_num_nodes);
    if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
      return false;
    }
    _progress_bar += (num_hns_before_pass - current_num_nodes);

    _timer.start_timer("contraction", "Contraction");
    // Perform parallel contraction
    _uncoarseningData.performMultilevelContraction(std::move(cluster_ids), false /* deterministic */, round_start);
    _timer.stop_timer("contraction");

    ++_pass_nr;
    return true;
  }

  template<bool has_fixed_vertices>
  HypernodeID performClustering(const Hypergraph& current_hg, vec<HypernodeID>& cluster_ids) {
    // We iterate in parallel over all vertices of the hypergraph and compute its contraction partner.
    // Matched vertices are linked in a concurrent union find data structure, that also aggregates
    // weights of the resulting clusters and keep track of the number of nodes left, if we would
    // contract all matched vertices.
    _timer.start_timer("clustering", "Clustering");
    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.start_timer("clustering_level_" + std::to_string(_pass_nr), "Level " + std::to_string(_pass_nr));
    }

    _rater.resetMatches();
    _rater.setCurrentNumberOfNodes(current_hg.initialNumNodes());
    const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
    NumNodesTracker num_nodes_tracker(current_hg.initialNumNodes() - current_hg.numRemovedHypernodes());
    ds::FixedVertexSupport<Hypergraph> fixed_vertices = current_hg.copyOfFixedVertexSupport();
    fixed_vertices.setMaxBlockWeight(_context.partition.max_part_weights);
    DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
    tbb::parallel_for(0U, current_hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < _current_vertices.size());
      const HypernodeID hn = _current_vertices[id];
      // We perform rating if ...
      //  1.) The contraction limit of the current level is not reached
      //  2.) Vertex hn is not matched before
      if (current_hg.nodeIsEnabled(hn)
          && num_nodes_tracker.currentNumNodes() > hierarchy_contraction_limit
          && _clustering_data.vertexIsUnmatched(hn)) {
        const Rating rating = _rater.template rate<ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy, has_fixed_vertices>(
                                     current_hg, hn, cluster_ids, _clustering_data.clusterWeight(),
                                     fixed_vertices, _context.coarsening.max_allowed_node_weight);
        if (rating.target != kInvalidHypernode) {
          bool success = _clustering_data.template matchVertices<has_fixed_vertices>(
            current_hg, hn, rating.target, cluster_ids, _rater, fixed_vertices);
          if (success) {
            // update the number of nodes in a way that minimizes synchronization overhead
            num_nodes_tracker.subtractNode(_context.shared_memory.original_num_threads, hierarchy_contraction_limit);
          }
        }
      }
    });

    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.stop_timer("clustering_level_" + std::to_string(_pass_nr));
    }
    _timer.stop_timer("clustering");
    if constexpr ( has_fixed_vertices ) {
      ASSERT(fixed_vertices.verifyClustering(current_hg, cluster_ids), "Fixed vertex support is corrupted");
    }

    return num_nodes_tracker.finalNumNodes();
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

  using Base::_hg;
  using Base::_context;
  using Base::_timer;
  using Base::_uncoarseningData;
  MultilevelVertexPairRater _rater;
  ConcurrentClusteringData _clustering_data;
  HypernodeID _initial_num_nodes;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  int _pass_nr;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
