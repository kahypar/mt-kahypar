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
#include "mt-kahypar/partition/coarsening/community_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/community_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/hypergraph_pruner.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/policies/community_assignment_objective.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <typename TypeTraits,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched,
          class CommunityAssignmentObjective = PinObjectivePolicy>
class CommunityCoarsenerT : public ICoarsener,
                            private CommunityCoarsenerBase<TypeTraits> {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  using Base = CommunityCoarsenerBase<TypeTraits>;
  using Memento = typename StreamingHyperGraph::Memento;
  using HypergraphPruner = HypergraphPrunerT<TypeTraits>;
  using Rater = CommunityVertexPairRater<TypeTraits,
                                         ScorePolicy,
                                         HeavyNodePenaltyPolicy,
                                         AcceptancePolicy>;
  using Rating = typename Rater::Rating;

  static constexpr bool debug = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

 public:
  CommunityCoarsenerT(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    Base(hypergraph, context, task_group_id),
    _progress_bar(hypergraph.initialNumNodes(), 0, false),
    _enable_randomization(true) { }

  CommunityCoarsenerT(const CommunityCoarsenerT&) = delete;
  CommunityCoarsenerT(CommunityCoarsenerT&&) = delete;
  CommunityCoarsenerT & operator= (const CommunityCoarsenerT &) = delete;
  CommunityCoarsenerT & operator= (CommunityCoarsenerT &&) = delete;

  ~CommunityCoarsenerT() = default;

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void coarsenImpl() override {
    // Initialize community coarsener
    this->init();

    // Compute execution order
    utils::Timer::instance().start_timer("compute_execution_order", "Compute Execution Order");
    std::vector<parallel::scalable_vector<HypernodeID> > community_hns(_hg.numCommunities());
    std::vector<HypernodeWeight> community_weights(_hg.numCommunities(), 0);
    for (const HypernodeID& hn : _hg.nodes()) {
      ASSERT(_hg.communityID(hn) < _hg.numCommunities());
      community_hns[_hg.communityID(hn)].emplace_back(hn);
      community_weights[_hg.communityID(hn)] += _hg.nodeWeight(hn);
    }

    // Execute coarsening in decreasing order of their community assignment objective
    // Note, this objective is the same than in the preprocessing step where we
    // redistribute the communities to numa nodes
    std::vector<PartitionID> community_ids(_hg.numCommunities(), 0);
    std::iota(community_ids.begin(), community_ids.end(), 0);
    std::sort(community_ids.begin(), community_ids.end(), [&](const PartitionID& lhs, const PartitionID& rhs) {
        return CommunityAssignmentObjective::objective(_hg, lhs) > CommunityAssignmentObjective::objective(_hg, rhs);
      });
    int used_numa_nodes = TBB::instance().num_used_numa_nodes();
    std::vector<tbb::concurrent_queue<PartitionID> > community_queues(used_numa_nodes);
    for (const PartitionID& community_id : community_ids) {
      int node = _hg.communityNumaNode(community_id);
      ASSERT(node < used_numa_nodes);
      community_queues[node].push(community_id);
    }
    utils::Timer::instance().stop_timer("compute_execution_order");

    // Parallel Community Coarsening
    // We schedule exactly number of available threads tasks. Each task
    // polls from the tbb::concurrent_queue (responsible for the numa node
    // where the task is executed) a community and executes the contractions
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar ) {
      _progress_bar.enable();
    }
    utils::Timer::instance().start_timer("parallel_community_coarsening", "Parallel Community Coarsening");
    for (int node = 0; node < used_numa_nodes; ++node) {
      int num_threads = TBB::instance().number_of_threads_on_numa_node(node);
      for (int i = 0; i < num_threads; ++i) {
        TBB::instance().numa_task_arena(node).execute([&, node] {
            TBB::instance().numa_task_group(_task_group_id, node).run([&, node] {
              tbb::concurrent_queue<PartitionID>& queue = community_queues[node];
              while (!queue.empty()) {
                PartitionID community_id = -1;
                bool success = queue.try_pop(community_id);
                if (success) {
                  ASSERT(community_id >= 0);
                  ASSERT(community_id < _hg.numCommunities());
                  if (community_hns[community_id].size() <= 1) {
                    continue;
                  }
                  // Compute contraction limit for community relative to
                  // community size and original contraction limit
                  HypernodeID contraction_limit =
                    std::ceil((((double)community_weights[community_id]) /
                               this->_hg.totalWeight()) * _context.coarsening.contraction_limit);
                  parallelCommunityCoarsening(community_id, contraction_limit, community_hns[community_id]);
                }
              }
            });
          });
      }
    }
    TBB::instance().wait(_task_group_id);
    _progress_bar += (_hg.initialNumNodes() - _progress_bar.count());
    _progress_bar.disable();
    utils::Timer::instance().stop_timer("parallel_community_coarsening");

    // Finalize community coarsening
    this->finalize();
  }

  void parallelCommunityCoarsening(const PartitionID community_id,
                                   const HypernodeID contraction_limit,
                                   parallel::scalable_vector<HypernodeID>& community_nodes) {
    ASSERT(community_nodes.size() > 1);
    DBG << "Start coarsening of community" << community_id
        << "with" << _hg.numCommunityHypernodes(community_id) << "vertices"
        << "and community assignment objective of" << CommunityAssignmentObjective::objective(_hg, community_id)
        << "and contraction limit" << contraction_limit
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu();

    // Sort community nodes in increasing order of their local community node id
    // Note, this is required by vertex pair rater that the node ids are sorted
    // in that order
    std::sort(community_nodes.begin(), community_nodes.end(),
              [&](const HypernodeID& lhs, const HypernodeID& rhs) {
        return _hg.communityNodeId(lhs) < _hg.communityNodeId(rhs);
      });
    parallel::scalable_vector<HypernodeID> nodes;
    for (const HypernodeID& hn : community_nodes) {
      nodes.emplace_back(hn);
    }

    Rater rater(_hg, _context, community_id, community_nodes);
    HypernodeID current_num_nodes = nodes.size();
    size_t pass_nr = 0;
    while (current_num_nodes > contraction_limit) {
      DBG << V(community_id) << V(pass_nr) << V(current_num_nodes) << V(sched_getcpu());

      rater.resetMatches();
      size_t num_hns_before_pass = current_num_nodes;
      parallel::scalable_vector<HypernodeID> tmp_nodes;

      if (_enable_randomization) {
        utils::Randomize::instance().shuffleVector(nodes, nodes.size(), sched_getcpu());
      }

      for (const HypernodeID& hn : nodes) {
        if (_hg.nodeIsEnabled(hn) && !_hg.isHighDegreeVertex(hn)) {
          Rating rating = rater.rate(hn);

          if (rating.target != kInvalidHypernode) {
            rater.markAsMatched(hn);
            rater.markAsMatched(rating.target);
            if (_hg.isHighDegreeVertex(rating.target)) {
              this->performContraction(rating.target, hn);
            } else {
              this->performContraction(hn, rating.target);
              tmp_nodes.emplace_back(hn);
            }
            --current_num_nodes;
          }

          if (current_num_nodes <= contraction_limit) {
            break;
          }
        }
      }

      _progress_bar += (num_hns_before_pass - current_num_nodes);

      if (num_hns_before_pass == current_num_nodes) {
        break;
      }

      nodes = std::move(tmp_nodes);
      ++pass_nr;
    }

    utils::Stats::instance().update_stat("num_removed_single_node_hes",
                                         (int64_t)_pruner[community_id].removedSingleNodeHyperedges().size());
  }

  bool uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation) override {
    return this->doUncoarsen(label_propagation);
  }

  using Base::_hg;
  using Base::_context;
  using Base::_pruner;
  using Base::_task_group_id;
  utils::ProgressBar _progress_bar;
  bool _enable_randomization;
};

template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched,
          class CommunityAssignmentObjective = PinObjectivePolicy>
using CommunityCoarsener = CommunityCoarsenerT<GlobalTypeTraits, ScorePolicy, HeavyNodePenaltyPolicy,
                                               AcceptancePolicy, CommunityAssignmentObjective>;
}  // namespace mt_kahypar
