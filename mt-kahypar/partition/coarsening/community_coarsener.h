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
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/hypergraph_pruner.h"
#include "mt-kahypar/partition/coarsening/community_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"

namespace mt_kahypar {

template< typename TypeTraits,
          class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched<> >
class CommunityCoarsenerT : public ICoarsener {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

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
    _hg(hypergraph),
    _context(context),
    _pruner() {
    for ( PartitionID community_id = 0; community_id < _hg.numCommunities(); ++community_id ) {
      _pruner.emplace_back(_hg.initialNumCommunityHypernodes(community_id));
    }
  }

  CommunityCoarsenerT(const CommunityCoarsenerT&) = delete;
  CommunityCoarsenerT(CommunityCoarsenerT&&) = delete;
  CommunityCoarsenerT& operator= (const CommunityCoarsenerT&) = delete;
  CommunityCoarsenerT& operator= (CommunityCoarsenerT&&) = delete;

  ~CommunityCoarsenerT() = default;


 private:
  void coarsenImpl() override {
    // Initialize community hyperedges such that parallel contractions in each community are possible
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    _hg.initializeCommunityHyperedges();
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_community_hyperedges", "Initialize Community Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 1, std::chrono::duration<double>(end - start).count());

    // Compute execution order
    start = std::chrono::high_resolution_clock::now();
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
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("compute_execution_order", "Compute Execution Order",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 2, std::chrono::duration<double>(end - start).count());

    // Parallel Community Coarsening
    start = std::chrono::high_resolution_clock::now();
    tbb::task_group group;
    parallel::scalable_vector<parallel::scalable_vector<Memento>> mementos(_hg.numCommunities());
    for ( size_t i = 0; i < community_hns.size(); ++i ) {
      if ( community_hns[i].size() <= 1 ) continue;

      PartitionID community_id = _hg.communityID(community_hns[i][0]);
      int node = _hg.communityNumaNode(community_id);
      TBB::instance().numa_task_arena(node).execute([&, community_id, i] {
        group.run([&, community_id, i] {
          HypernodeID contraction_limit =
            std::ceil((((double) this->_hg.initialNumCommunityHypernodes(community_id)) /
              this->_hg.totalWeight()) * _context.coarsening.contraction_limit);
          parallelCommunityCoarsening(community_id, contraction_limit, community_hns[i], mementos[community_id]);
        });
      });
    }
    group.wait();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("parallel_community_coarsening", "Parallel Community Coarsening",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 3, std::chrono::duration<double>(end - start).count());

    // Merge Coarsening Mementos
    start = std::chrono::high_resolution_clock::now();
    std::vector<Memento> history = mergeContractions(mementos);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("merge_contractions", "Merge Contractions",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 4, std::chrono::duration<double>(end - start).count());

    // Reset community hyperedges
    start = std::chrono::high_resolution_clock::now();
    _hg.resetCommunityHyperedges(history);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("reset_community_hyperedges", "Reset Community Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 5, std::chrono::duration<double>(end - start).count());
  }

  void parallelCommunityCoarsening(const PartitionID community_id,
                                   const HypernodeID contraction_limit,
                                   parallel::scalable_vector<HypernodeID>& community_nodes,
                                   parallel::scalable_vector<Memento>& mementos) {
    ASSERT(community_nodes.size() > 1);
    LOG << "Start coarsening of community" << community_id
        << "with" << _hg.initialNumCommunityHypernodes(community_id) << "vertices"
        << "and contraction limit" << contraction_limit
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu();

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
      utils::Randomize::instance().shuffleVector(nodes, nodes.size(), sched_getcpu());

      for ( const HypernodeID& hn : nodes ) {
        if ( _hg.nodeIsEnabled(hn) ) {
          Rating rating = rater.rate(hn);

          if ( rating.target != kInvalidHypernode ) {
            ASSERT(_hg.communityID(hn) == community_id);
            ASSERT(_hg.communityID(rating.target) == community_id);

            DBG << "Contract: (" << hn << "," << rating.target << ")";
            rater.markAsMatched(hn);
            rater.markAsMatched(rating.target);
            mementos.emplace_back(_hg.contract(hn, rating.target, community_id));
            _pruner[community_id].removeSingleNodeHyperedges(_hg, community_id, mementos.back());
            _pruner[community_id].removeParallelHyperedges(_hg, community_id, mementos.back());
            tmp_nodes.emplace_back(hn);
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
  }

  std::vector<Memento> mergeContractions(parallel::scalable_vector<parallel::scalable_vector<Memento>>& mementos) {
    std::vector<Memento> history;
    std::vector<size_t> pos(mementos.size(), 0);
    bool found = false;
    do {
      found = false;
      for ( size_t i = 0; i < mementos.size(); ++i ) {
        if ( pos[i] < mementos[i].size() ) {
          found = true;
          history.emplace_back(std::move(mementos[i][pos[i]]));
          ++pos[i];
        }
      }
    } while ( found );

    return history;
  }

  bool uncoarsenImpl() override {
    return false;
  }

  HyperGraph& _hg;
  const Context& _context;
  parallel::scalable_vector<HypergraphPruner> _pruner;
};

template< class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = MultiplicativePenalty,
          class AcceptancePolicy = BestRatingPreferringUnmatched<> >
using CommunityCoarsener = CommunityCoarsenerT<GlobalTypeTraits, ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy>;

}  // namespace mt_kahypar
