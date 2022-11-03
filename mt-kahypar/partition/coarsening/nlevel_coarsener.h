/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
#include "tbb/parallel_sort.h"
#include "tbb/parallel_scan.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/nlevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched>
class NLevelCoarsener : public ICoarsener,
                        private NLevelCoarsenerBase {
 private:

  #define HIGH_DEGREE_VERTEX_THRESHOLD ID(200000)

  using Base = NLevelCoarsenerBase;
  using Rater = NLevelVertexPairRater<ScorePolicy,
                                      HeavyNodePenaltyPolicy,
                                      AcceptancePolicy>;
  using Rating = typename Rater::Rating;

  class ContractionLimitTracker {

   public:
    explicit ContractionLimitTracker(const Context& context) :
      _context(context),
      _initial_num_nodes(0),
      _current_num_nodes(0),
      _contracted_nodes(0),
      _num_nodes_update_threshold(0) { }

    void initialize(const HypernodeID initial_num_nodes) {
      _initial_num_nodes = initial_num_nodes;
      _current_num_nodes = initial_num_nodes;
    }

    HypernodeID currentNumNodes() const {
      return _current_num_nodes;
    }

    void update(const HypernodeID num_contractions, const HypernodeID contraction_limit) {
      if ( num_contractions > 0 ) {
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
        HypernodeID& local_contracted_nodes = _contracted_nodes.local();
        local_contracted_nodes += num_contractions;
        if (  local_contracted_nodes >= _num_nodes_update_threshold.local() ) {
          _current_num_nodes = _initial_num_nodes -
            _contracted_nodes.combine(std::plus<HypernodeID>());
          const HypernodeID dist_to_cl = _current_num_nodes > contraction_limit ?
            _current_num_nodes - contraction_limit : 0;
          _num_nodes_update_threshold.local() +=
            dist_to_cl / _context.shared_memory.num_threads;
        }
      }
    }

    void updateCurrentNumNodes() {
      _current_num_nodes = _initial_num_nodes -
        _contracted_nodes.combine(std::plus<HypernodeID>());
    }

   private:
    const Context& _context;
    HypernodeID _initial_num_nodes;
    HypernodeID _current_num_nodes;
    tbb::enumerable_thread_specific<HypernodeID> _contracted_nodes;
    tbb::enumerable_thread_specific<HypernodeID> _num_nodes_update_threshold;
  };

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  NLevelCoarsener(Hypergraph& hypergraph,
                  const Context& context,
                  UncoarseningData& uncoarseningData) :
    Base(hypergraph, context, uncoarseningData),
    _rater(hypergraph, context),
    _current_vertices(),
    _tmp_current_vertices(),
    _enabled_vertex_flag_array(),
    _cl_tracker(context),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += hypergraph.numRemovedHypernodes();
    tbb::parallel_invoke([&] {
      _current_vertices.resize(hypergraph.initialNumNodes());
      tbb::parallel_for(ID(0), hypergraph.initialNumNodes(), [&](const HypernodeID hn) {
        _current_vertices[hn] = hn;
      });
      utils::Randomize::instance().parallelShuffleVector(_current_vertices, 0UL, _current_vertices.size());
    }, [&] {
      _tmp_current_vertices.resize(hypergraph.initialNumNodes());
    }, [&] {
      _enabled_vertex_flag_array.resize(hypergraph.initialNumNodes());
    });
  }

  NLevelCoarsener(const NLevelCoarsener&) = delete;
  NLevelCoarsener(NLevelCoarsener&&) = delete;
  NLevelCoarsener & operator= (const NLevelCoarsener &) = delete;
  NLevelCoarsener & operator= (NLevelCoarsener &&) = delete;

  ~NLevelCoarsener() = default;

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void coarsenImpl() override {
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar ) {
      _progress_bar.enable();
    }

    // Contraction function
    auto contraction = [&](const HypernodeID& hn) {
      HypernodeID num_contractions = 0;
      if ( _hg.nodeIsEnabled(hn) ) {
        const Rating rating = _rater.rate(_hg, hn, _max_allowed_node_weight);
        if ( rating.target != kInvalidHypernode ) {
          HypernodeID u = hn;
          HypernodeID v = rating.target;
          // In case v is a high degree vertex, we reverse contraction order to improve performance
          if ( _hg.nodeDegree(u) < _hg.nodeDegree(v) && _hg.nodeDegree(v) > HIGH_DEGREE_VERTEX_THRESHOLD ) {
            u = rating.target;
            v = hn;
          }

          if ( _hg.registerContraction(u, v) ) {
            _rater.markAsMatched(u);
            _rater.markAsMatched(v);
            num_contractions = _hg.contract(v, _max_allowed_node_weight);
            _progress_bar += num_contractions;
          }
        }
      }
      return num_contractions;
    };

    const HypernodeID initial_num_nodes = _hg.initialNumNodes() - _hg.numRemovedHypernodes();
    _cl_tracker.initialize(initial_num_nodes);
    int pass_nr = 0;
    while ( _cl_tracker.currentNumNodes() > _context.coarsening.contraction_limit ) {
      DBG << V(pass_nr) << V(_cl_tracker.currentNumNodes());
      const HypernodeID num_hns_before_pass = _cl_tracker.currentNumNodes();

      // Coarsening Pass
      _rater.resetMatches();
      coarseningPass(contraction);

      // Writes all enabled vertices to _current_vertices
      _cl_tracker.updateCurrentNumNodes();
      compactifyVertices();
      utils::Randomize::instance().parallelShuffleVector(
        _current_vertices, 0UL, _current_vertices.size());

      // Terminate contraction if the number of contracted vertices in this round
      // is smaller than a certain fraction.
      const double reduction_vertices_percentage =
        static_cast<double>(num_hns_before_pass) /
        static_cast<double>(_cl_tracker.currentNumNodes());
      if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
        break;
      }

      ++pass_nr;
    }

    _progress_bar += (initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();
    _uncoarseningData.finalizeCoarsening();
  }


  template<class F>
  void coarseningPass(const F& contraction) {
    double contraction_limit =
      std::max(static_cast<HypernodeID>(_cl_tracker.currentNumNodes() /
        _context.coarsening.maximum_shrink_factor), _context.coarsening.contraction_limit);
    if ( _context.coarsening.maximum_shrink_factor > 99.0 ) {
      contraction_limit = _context.coarsening.contraction_limit;
    }

    HighResClockTimepoint round_start = std::chrono::high_resolution_clock::now();
    _timer.start_timer("clustering", "Clustering");
    tbb::parallel_for(0UL, _current_vertices.size(), [&](const size_t i) {
      if ( _cl_tracker.currentNumNodes() > contraction_limit ) {
        const HypernodeID& hn = _current_vertices[i];
        const HypernodeID num_contractions = contraction(hn);
        _cl_tracker.update(num_contractions, contraction_limit);
      }
    });
    _timer.stop_timer("clustering");

    // Remove single-pin and parallel nets
    Base::removeSinglePinAndParallelNets(round_start);
    HEAVY_COARSENING_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::compactifiedHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::compactifiedPartitionedHypergraph();
  }

  void compactifyVertices() {
    // Mark all vertices that are still enabled
    const HypernodeID current_num_nodes = _cl_tracker.currentNumNodes();
    tbb::parallel_for(0UL, _current_vertices.size(), [&](const size_t i) {
      const HypernodeID hn = _current_vertices[i];
      _enabled_vertex_flag_array[i] = _hg.nodeIsEnabled(hn);
    });

    // Calculate prefix sum over all enabled vertices to determine their new position
    // in _current_vertices
    parallel::TBBPrefixSum<size_t> active_vertex_prefix_sum(_enabled_vertex_flag_array);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, _enabled_vertex_flag_array.size()), active_vertex_prefix_sum);
    ASSERT(active_vertex_prefix_sum.total_sum() == static_cast<size_t>(current_num_nodes));

    // Write all enabled vertices to _tmp_current_vertices
    _tmp_current_vertices.resize(current_num_nodes);
    tbb::parallel_for(0UL, _current_vertices.size(), [&](const size_t i) {
      const HypernodeID hn = _current_vertices[i];
      if ( _hg.nodeIsEnabled(hn) ) {
        const size_t pos = active_vertex_prefix_sum[i];
        ASSERT(pos < _tmp_current_vertices.size());
        _tmp_current_vertices[pos] = hn;
      }
    });
    _current_vertices.swap(_tmp_current_vertices);
    _enabled_vertex_flag_array.resize(current_num_nodes);
    ASSERT(_current_vertices.size() == static_cast<size_t>(current_num_nodes));
  }

  using Base::_hg;
  Rater _rater;
  parallel::scalable_vector<HypernodeID> _current_vertices;
  parallel::scalable_vector<HypernodeID> _tmp_current_vertices;
  parallel::scalable_vector<size_t> _enabled_vertex_flag_array;
  ContractionLimitTracker _cl_tracker;
  HypernodeWeight _max_allowed_node_weight;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
