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
#include <queue>

#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/coarsening/hypergraph_pruner.h"

namespace mt_kahypar {

template< typename TypeTraits >
class CommunityCoarsenerBase {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;

  using Memento = typename StreamingHyperGraph::Memento;
  using HypergraphPruner = HypergraphPrunerT<TypeTraits>;

  static constexpr bool debug = false;
  static HypernodeID kInvalidHyperedge;

 public:
  CommunityCoarsenerBase(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _init(false),
    _contained_hypernodes(hypergraph.initialNumNodes()),
    _parallel_he_representative(hypergraph.initialNumEdges(), kInvalidHyperedge),
    _pruner(),
    _community_history(_hg.numCommunities()),
    _history() {
    for ( PartitionID community_id = 0; community_id < _hg.numCommunities(); ++community_id ) {
      _pruner.emplace_back(_contained_hypernodes, _parallel_he_representative);
    }
  }

  CommunityCoarsenerBase(const CommunityCoarsenerBase&) = delete;
  CommunityCoarsenerBase(CommunityCoarsenerBase&&) = delete;
  CommunityCoarsenerBase& operator= (const CommunityCoarsenerBase&) = delete;
  CommunityCoarsenerBase& operator= (CommunityCoarsenerBase&&) = delete;

  virtual ~CommunityCoarsenerBase() = default;

 protected:
  void init() {
    // Initialize community hyperedges such that parallel contractions in each community are possible
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    _hg.initializeCommunityHyperedges();
    _init = true;
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("initialize_community_hyperedges", "Initialize Community Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 1, std::chrono::duration<double>(end - start).count());
  }

  void finalize() {
    // Merge Coarsening Mementos
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    mergeCommunityContractions();
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("merge_contractions", "Merge Contractions",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 4, std::chrono::duration<double>(end - start).count());

    // Reset community hyperedges
    start = std::chrono::high_resolution_clock::now();
    _hg.removeCommunityHyperedges(_history);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("reset_community_hyperedges", "Reset Community Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 5, std::chrono::duration<double>(end - start).count());

    start = std::chrono::high_resolution_clock::now();
    postprocessParallelHyperedges();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("postprocess_parallel_hyperedges", "Postprocess Parallel Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 6, std::chrono::duration<double>(end - start).count());

    _init = false;
  }

  void performContraction(const HypernodeID u, const HypernodeID v) {
    ASSERT(_init, "Community coarsener not initialized");
    ASSERT(_hg.communityID(u) == _hg.communityID(v));
    PartitionID community_id = _hg.communityID(u);
    DBG << "Contract: (" << u << "," << v << ")";
    _community_history[community_id].emplace_back(_hg.contract(u, v, community_id));
    _pruner[community_id].removeSingleNodeHyperedges(_hg, community_id, _community_history[community_id].back());
    _pruner[community_id].removeParallelHyperedges(_hg, community_id, _community_history[community_id].back());
  }

  bool doUncoarsen() {
    ASSERT(!_init, "Community coarsener must be finalized before uncoarsening");

    while ( !_history.empty() ) {
      PartitionID community_id = _hg.communityID(_history.back().u);
      DBG << "Uncontracting: (" << _history.back().u << "," << _history.back().v << ")" << V(_history.size());
      _pruner[community_id].restoreParallelHyperedges(_hg, _history.back());
      _pruner[community_id].restoreSingleNodeHyperedges(_hg, _history.back());
      _hg.uncontract(_history.back(), _parallel_he_representative);
      _hg.updateGlobalPartInfos();
      _history.pop_back();
    }

    return true;
  }

 private:
  void mergeCommunityContractions() {
    size_t size = 0;
    for ( const auto& community_mementos : _community_history ) {
      size += community_mementos.size();
    }
    std::sort(_community_history.begin(), _community_history.end(),
              [&](const auto& lhs, const auto& rhs) {
      return lhs.size() > rhs.size();
    });
    while ( _community_history.back().size() == 0 ) {
      _community_history.pop_back();
    }

    // Mementos from all communities are written in parallel
    // to one history vector in round-robin fashion.
    _history.resize(size);
    tbb::parallel_for(tbb::blocked_range<size_t>(0UL, _community_history.size()),
      [&](const tbb::blocked_range<size_t>& range) {
      for ( size_t i = range.begin(); i < range.end(); ++i ) {
        size_t current_pos = i;
        int64_t k = _community_history.size();
        for ( size_t pos = 0; pos < _community_history[i].size(); ++pos ) {
          while ( pos >= _community_history[k - 1].size() ) {
            --k;
            ASSERT(k > 0);
          }

          ASSERT(current_pos < _history.size(), V(i) << V(current_pos) << V(_history.size()));
          ASSERT(_history[current_pos].u == std::numeric_limits<HypernodeID>::max(), V(current_pos) << V(i) << V(k));
          _history[current_pos] = std::move(_community_history[i][pos]);
          current_pos += k;
        }
      }
    });

    // Release memory of community histories
    for ( auto& community_mementos : _community_history ) {
      parallel::scalable_vector<Memento> tmp;
      community_mementos = std::move(tmp);
    }

    ASSERT([&] {
      for ( size_t i = 0; i < _history.size(); ++i ) {
        if ( _history[i].u == std::numeric_limits<HypernodeID>::max() ) {
          return false;
        }
      }
      return true;
    }(), "Merging contraction hierarchies failed");
  }

  void postprocessParallelHyperedges() {
    // During parallel community coarsening we detect parallel hyperedges
    // by checking if two hyperedges become parallel in all its community
    // hyperedges. This check is performed in each community hyperegraph pruner.
    // Thus it can happen that several pruner detect that two hyperedges become
    // parallel. Once two parallel hyperedges are detected within a pruner, the
    // hyperedge is only removed in the community for which the pruner is responsible
    // for. The representative of the removed hyperedge is stored in the global
    // array _parallel_he_representative (which is shared). It can happen that several
    // pruner remove the same hyperedge, but with different representatives. In that
    // case, the representative of the hyperedge is the one which is last written
    // to the array (last write wins, no concurrent access control here).
    // The removal of parallel hyperedges together with their representatives
    // build a forest in _parallel_he_representatives. In a postprocessing step
    // we determine the weight of each hyperedge in the forest by performing
    // a buttom up traversal (from the leaves). The weight of each leaf hyperedge
    // is it original weight and the weight of an inner hyperedge is the weight
    // of its childs. All hyperedges except the roots are disabled.

    // Determine the in degree of all tree nodes
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    std::vector<HyperedgeID> in_degree(_hg.initialNumEdges(), 0);
    for ( HyperedgeID he = 0; he < _hg.initialNumEdges(); ++he ) {
      HyperedgeID representative = _parallel_he_representative[he];
      if ( representative != kInvalidHyperedge ) {
        ASSERT(representative < in_degree.size());
        ++in_degree[representative];
      }
    }

    // Push all leaves into a queue
    std::queue<HyperedgeID> q;
    for ( HyperedgeID he = 0; he < _hg.initialNumEdges(); ++he ) {
      if ( in_degree[he] == 0 && _parallel_he_representative[he] != kInvalidHyperedge ) {
        q.push(he);
      }
    }

    // Perform a bottom up traversal of the forest and compute weight
    // of each hyperedge
    std::vector<HyperedgeWeight> weights(_hg.initialNumEdges(), 0);
    while ( !q.empty() ) {
      HyperedgeID original_he = q.front();
      HyperedgeID he = _hg.globalEdgeID(original_he);

      bool was_enabled = false;
      if ( !_hg.edgeIsEnabled(he) ) {
        // Only roots are allowed to be disabled, because of single-pin hyperedge removal.
        ASSERT(_parallel_he_representative[original_he] == kInvalidHyperedge);
        // We temporary enable the hyperedge in order to set the edge weight
        was_enabled = true;
        _hg.enableHyperedge(he);
      }

      HyperedgeWeight weight = weights[original_he] + _hg.edgeWeight(he);
      _hg.setEdgeWeight(he, weight);
      q.pop();

      HyperedgeID representative = _parallel_he_representative[original_he];
      if ( representative != kInvalidHyperedge ) {
        ASSERT(representative < weights.size());
        weights[representative] += _hg.edgeWeight(he);
        --in_degree[representative];
        if ( in_degree[representative] == 0 ) {
          q.push(representative);
        }
      }

      // If the hyperedge was temporary enabled (to set edge weight) or
      // the hyperedge is an inner node of the parallel hyperedge forest
      // we disable the hyperedge
      if ( was_enabled || representative != kInvalidHyperedge ) {
        ASSERT(_hg.edgeIsEnabled(he));
        _hg.disableHyperedge(he);
      }
    }
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("determine_he_weights", "Determine HE weights",
      "postprocess_parallel_hyperedges", mt_kahypar::utils::Timer::Type::COARSENING, 1, std::chrono::duration<double>(end - start).count());


    start = std::chrono::high_resolution_clock::now();
    _hg.invalidateDisabledHyperedgesFromIncidentNets();
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("remove_disabled_hyperedges_from_incident_nets", "Remove Disabled HE from HNs",
      "postprocess_parallel_hyperedges", mt_kahypar::utils::Timer::Type::COARSENING, 2, std::chrono::duration<double>(end - start).count());
  }

 protected:
  HyperGraph& _hg;
  const Context& _context;
  bool _init;
  ThreadLocalFastResetFlagArray _contained_hypernodes;
  parallel::scalable_vector<HyperedgeID> _parallel_he_representative;
  parallel::scalable_vector<HypergraphPruner> _pruner;
  parallel::scalable_vector<parallel::scalable_vector<Memento>> _community_history;
  std::vector<Memento> _history;
};

template< typename TypeTraits >
HyperedgeID CommunityCoarsenerBase<TypeTraits>::kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();


}  // namespace mt_kahypar
