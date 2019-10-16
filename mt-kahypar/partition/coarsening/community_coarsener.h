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

#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/community_coarsener_base.h"
#include "mt-kahypar/partition/coarsening/hypergraph_pruner.h"
#include "mt-kahypar/partition/coarsening/community_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"

namespace mt_kahypar {

template< typename TypeTraits,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched >
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
  CommunityCoarsenerT(HyperGraph& hypergraph, const Context& context) :
    Base(hypergraph, context),
    _enable_randomization(true) { }

  CommunityCoarsenerT(const CommunityCoarsenerT&) = delete;
  CommunityCoarsenerT(CommunityCoarsenerT&&) = delete;
  CommunityCoarsenerT& operator= (const CommunityCoarsenerT&) = delete;
  CommunityCoarsenerT& operator= (CommunityCoarsenerT&&) = delete;

  ~CommunityCoarsenerT() = default;

  void disableRandomization() {
    _enable_randomization = false;
  }

 private:
  void coarsenImpl() override {
    // Initialize community coarsener
    this->init();

    // Compute execution order
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    std::vector<parallel::scalable_vector<HypernodeID>> community_hns(_hg.numCommunities());
    for ( const HypernodeID& hn : _hg.nodes() ) {
      ASSERT(_hg.communityID(hn) < _hg.numCommunities());
      community_hns[_hg.communityID(hn)].emplace_back(hn);
    }
    // Execute coarsening in decreasing order of their community sizes
    std::sort(community_hns.begin(), community_hns.end(),
      [&](const parallel::scalable_vector<HypernodeID>& lhs, const parallel::scalable_vector<HypernodeID>& rhs) {
      return lhs.size() > rhs.size();
    });
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_execution_order", "Compute Execution Order",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 2, std::chrono::duration<double>(end - start).count());

    // Parallel Community Coarsening
    start = std::chrono::high_resolution_clock::now();
    for ( size_t i = 0; i < community_hns.size(); ++i ) {
      if ( community_hns[i].size() <= 1 ) continue;

      PartitionID community_id = _hg.communityID(community_hns[i][0]);
      int node = _hg.communityNumaNode(community_id);
      TBB::instance().numa_task_arena(node).execute([&, community_id, i] {
        TBB::instance().numa_task_group(node).run([&, community_id, i] {
          // Compute contraction limit for community relative to
          // community size and original contraction limit
          HypernodeID contraction_limit =
            std::ceil((((double) this->_hg.numCommunityHypernodes(community_id)) /
              this->_hg.totalWeight()) * _context.coarsening.contraction_limit);
          parallelCommunityCoarsening(community_id, contraction_limit, community_hns[i]);
        });
      });
    }
    TBB::instance().wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("parallel_community_coarsening", "Parallel Community Coarsening",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 3, std::chrono::duration<double>(end - start).count());

    // Finalize community coarsening
    this->finalize();
  }

  void parallelCommunityCoarsening(const PartitionID community_id,
                                   const HypernodeID contraction_limit,
                                   parallel::scalable_vector<HypernodeID>& community_nodes) {
    ASSERT(community_nodes.size() > 1);
    DBG << "Start coarsening of community" << community_id
        << "with" << _hg.numCommunityHypernodes(community_id) << "vertices"
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
    for ( const HypernodeID& hn : community_nodes ) {
      nodes.emplace_back(hn);
    }

    Rater rater(_hg, _context, community_id, community_nodes);
    HypernodeID current_num_nodes = nodes.size();
    size_t pass_nr = 0;
    while ( current_num_nodes > contraction_limit ) {
      DBG << V(community_id) << V(pass_nr) << V(current_num_nodes) << V(sched_getcpu());

      rater.resetMatches();
      size_t num_hns_before_pass = current_num_nodes;
      parallel::scalable_vector<HypernodeID> tmp_nodes;

      if ( _enable_randomization ) {
        utils::Randomize::instance().shuffleVector(nodes, nodes.size(), sched_getcpu());
      }

      for ( const HypernodeID& hn : nodes ) {
        if ( _hg.nodeIsEnabled(hn) && _hg.nodeDegree(hn) <= _context.coarsening.hypernode_degree_threshold ) {
          Rating rating = rater.rate(hn);

          if ( rating.target != kInvalidHypernode ) {
            DBG << "Contract: (" << hn << "," << rating.target << ")";
            rater.markAsMatched(hn);
            rater.markAsMatched(rating.target);
            if ( _hg.nodeDegree(rating.target) < _context.coarsening.hypernode_degree_threshold ) {
              this->performContraction(hn, rating.target);
              tmp_nodes.emplace_back(hn);
            } else {
              this->performContraction(rating.target, hn);
            }
            --current_num_nodes;
          }

          if ( current_num_nodes <= contraction_limit ) {
            break;
          }
        }
      }

      if ( num_hns_before_pass == current_num_nodes ) {
        break;
      }

      nodes = std::move(tmp_nodes);
      ++pass_nr;
    }
    utils::Stats::instance().update_stat("num_removed_single_node_hes",
      (int64_t) _pruner[community_id].removedSingleNodeHyperedges().size());
  }

  bool uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation) override {
    return this->doUncoarsen(label_propagation);
  }

  using Base::_hg;
  using Base::_context;
  using Base::_pruner;
  bool _enable_randomization;
};

template< class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched >
using CommunityCoarsener = CommunityCoarsenerT<GlobalTypeTraits, ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy>;

}  // namespace mt_kahypar
