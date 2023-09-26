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
    _initial_num_nodes(utils::cast<Hypergraph>(hypergraph).initialNumNodes()),
    _current_vertices(),
    _clustering_data(_hg.initialNumNodes(), _context),
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
  HypernodeID performClustering(const Hypergraph& current_hg,
                                vec<HypernodeID>& cluster_ids) {
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
    const HypernodeID num_hns_before_pass = current_hg.initialNumNodes() - current_hg.numRemovedHypernodes();
    const HypernodeID hierarchy_contraction_limit = hierarchyContractionLimit(current_hg);
    DBG << V(current_hg.initialNumNodes()) << V(hierarchy_contraction_limit);
    HypernodeID current_num_nodes = num_hns_before_pass;
    tbb::enumerable_thread_specific<HypernodeID> contracted_nodes(0);
    tbb::enumerable_thread_specific<HypernodeID> num_nodes_update_threshold(0);
    ds::FixedVertexSupport<Hypergraph> fixed_vertices = current_hg.copyOfFixedVertexSupport();
    fixed_vertices.setMaxBlockWeight(_context.partition.max_part_weights);
    tbb::parallel_for(0U, current_hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < _current_vertices.size());
      const HypernodeID hn = _current_vertices[id];
      if (current_hg.nodeIsEnabled(hn)) {
        // We perform rating if ...
        //  1.) The contraction limit of the current level is not reached
        //  2.) Vertex hn is not matched before
        const HypernodeID u = hn;
        if (_clustering_data.vertexIsUnmatched(u) && current_num_nodes > hierarchy_contraction_limit) {
          ASSERT(current_hg.nodeIsEnabled(hn));
          const Rating rating = _rater.template rate<ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy, has_fixed_vertices>(
            current_hg, hn, cluster_ids, _clustering_data.clusterWeight(), fixed_vertices, _context.coarsening.max_allowed_node_weight);
          if (rating.target != kInvalidHypernode) {
            const HypernodeID v = rating.target;
            HypernodeID& local_contracted_nodes = contracted_nodes.local();
            _clustering_data.template matchVertices<has_fixed_vertices>(current_hg,
              u, v, cluster_ids, local_contracted_nodes, _rater, fixed_vertices);

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
            if (local_contracted_nodes >= num_nodes_update_threshold.local()) {
              current_num_nodes = num_hns_before_pass -
                                  contracted_nodes.combine(std::plus<HypernodeID>());
              const HypernodeID dist_to_contraction_limit =
                current_num_nodes > hierarchy_contraction_limit ?
                current_num_nodes - hierarchy_contraction_limit : 0;
              num_nodes_update_threshold.local() +=
                dist_to_contraction_limit / _context.shared_memory.original_num_threads;
            }
          }
        }
      }
    });
    if ( _context.partition.show_detailed_clustering_timings ) {
      _timer.stop_timer("clustering_level_" + std::to_string(_pass_nr));
    }
    _timer.stop_timer("clustering");

    if constexpr ( has_fixed_vertices ) {
      // Verify fixed vertices
      ASSERT([&] {
        vec<PartitionID> fixed_vertex_blocks(current_hg.initialNumNodes(), kInvalidPartition);
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          if ( current_hg.isFixed(hn) ) {
            if ( fixed_vertex_blocks[cluster_ids[hn]] != kInvalidPartition &&
                 fixed_vertex_blocks[cluster_ids[hn]] != current_hg.fixedVertexBlock(hn)) {
              LOG << "There are two nodes assigned to same cluster that belong to different fixed vertex blocks";
              return false;
            }
            fixed_vertex_blocks[cluster_ids[hn]] = current_hg.fixedVertexBlock(hn);
          }
        }

        vec<HypernodeWeight> expected_block_weights(_context.partition.k, 0);
        for ( const HypernodeID& hn : current_hg.nodes() ) {
          if ( fixed_vertex_blocks[cluster_ids[hn]] != kInvalidPartition ) {
            if ( !fixed_vertices.isFixed(cluster_ids[hn]) ) {
              LOG << "Cluster" << cluster_ids[hn] << "should be fixed to block"
                  << fixed_vertex_blocks[cluster_ids[hn]];
              return false;
            }
            expected_block_weights[fixed_vertex_blocks[cluster_ids[hn]]] += current_hg.nodeWeight(hn);
          }
        }

        for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
          if ( fixed_vertices.fixedVertexBlockWeight(block) != expected_block_weights[block] ) {
            LOG << "Fixed vertex block" << block << "should have weight" << expected_block_weights[block]
                << ", but it is" << fixed_vertices.fixedVertexBlockWeight(block);
            return false;
          }
        }
        return true;
      }(), "Fixed vertex support is corrupted");
    }

    return num_hns_before_pass - contracted_nodes.combine(std::plus<>());
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
  HypernodeID _initial_num_nodes;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  ConcurrentClusteringData _clustering_data;
  int _pass_nr;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
