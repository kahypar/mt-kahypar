/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/nlevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
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


  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  NLevelCoarsener(Hypergraph& hypergraph,
                  const Context& context,
                  const TaskGroupID task_group_id,
                  const bool top_level) :
    Base(hypergraph, context, task_group_id, top_level),
    _rater(hypergraph, context),
    _max_allowed_node_weight(context.coarsening.max_allowed_node_weight),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) {
    _progress_bar += hypergraph.numRemovedHypernodes();
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

    const HypernodeID initial_num_nodes = _hg.initialNumNodes() - _hg.numRemovedHypernodes();
    HypernodeID current_num_nodes = initial_num_nodes;
    tbb::enumerable_thread_specific<HypernodeID> contracted_nodes(0);
    tbb::enumerable_thread_specific<HypernodeID> num_nodes_update_threshold(0);
    parallel::scalable_vector<HypernodeID> current_hns(current_num_nodes, kInvalidHypernode);
    int pass_nr = 0;
    while ( current_num_nodes > _context.coarsening.contraction_limit ) {
      DBG << V(pass_nr) << V(current_num_nodes);

      utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
      std::atomic<size_t> idx(0);
      current_hns.resize(current_num_nodes);
      _hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const size_t i = idx++;
        ASSERT(i < current_hns.size());
        current_hns[i] = hn;
      });
      ASSERT(idx == current_num_nodes);
      utils::Randomize::instance().parallelShuffleVector(current_hns, 0UL, current_hns.size());
      utils::Timer::instance().stop_timer("random_shuffle");

      utils::Timer::instance().start_timer("n_level_coarsening", "n-Level Coarsening");
      const HypernodeID num_hns_before_pass = current_num_nodes;
      _rater.resetMatches();
      tbb::parallel_for(0UL, current_hns.size(), [&](const size_t i) {
        const HypernodeID hn = current_hns[i];
        if ( current_num_nodes > _context.coarsening.contraction_limit && _hg.nodeIsEnabled(hn) ) {
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
              // TODO(heuer): Think what should happen if a contraction failed due to the max node weight
              // It might be that other concurrent running contractions to that vertex are relinked to
              // an other vertex in the contraction tree.
              const size_t num_contractions = _hg.contract(v, _max_allowed_node_weight);
              _progress_bar += num_contractions; // TODO: should be updated outside this parallel for loop

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
              HypernodeID& local_contracted_nodes = contracted_nodes.local();
              local_contracted_nodes += num_contractions;
              if (  local_contracted_nodes >= num_nodes_update_threshold.local() ) {
                current_num_nodes = initial_num_nodes -
                  contracted_nodes.combine(std::plus<HypernodeID>());
                num_nodes_update_threshold.local() +=
                  (current_num_nodes - _context.coarsening.contraction_limit) /
                  _context.shared_memory.num_threads;
              }
            }
          }
        }
      });
      utils::Timer::instance().stop_timer("n_level_coarsening");

      // Remove single-pin and parallel nets
      Base::removeSinglePinAndParallelNets();
      current_num_nodes = initial_num_nodes -
        contracted_nodes.combine(std::plus<HypernodeID>());

      HEAVY_COARSENING_ASSERT(_hg.verifyIncidenceArrayAndIncidentNets());

      // Terminate contraction if the number of contracted vertices in this round
      // is smaller than a certain fraction.
      const double reduction_vertices_percentage =
        static_cast<double>(num_hns_before_pass) /
        static_cast<double>(current_num_nodes);
      if ( reduction_vertices_percentage <= _context.coarsening.minimum_shrink_factor ) {
        break;
      }

      ++pass_nr;
    }

    _progress_bar += (initial_num_nodes - _progress_bar.count());
    _progress_bar.disable();
    Base::finalize();
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::compactifiedHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::compactifiedPartitionedHypergraph();
  }

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) override {
    return Base::doUncoarsen(label_propagation, fm);
  }

  using Base::_hg;
  Rater _rater;
  HypernodeWeight _max_allowed_node_weight;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

}  // namespace mt_kahypar
