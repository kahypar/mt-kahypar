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

#include "mt-kahypar/definitions.h"
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

 public:
  CommunityCoarsenerBase(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _init(false),
    _pruner(),
    _community_history(_hg.numCommunities()),
    _history() {
    for ( PartitionID community_id = 0; community_id < _hg.numCommunities(); ++community_id ) {
      _pruner.emplace_back(_hg.initialNumCommunityHypernodes(community_id));
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
    _hg.resetCommunityHyperedges(_history);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("reset_community_hyperedges", "Reset Community Hyperedges",
      "coarsening", mt_kahypar::utils::Timer::Type::COARSENING, 5, std::chrono::duration<double>(end - start).count());

    _init = false;
  }

  void performContraction(const HypernodeID u, const HypernodeID v) {
    ASSERT(_init, "Community coarsener not initialized");
    ASSERT(_hg.communityID(u) == _hg.communityID(v));
    PartitionID community_id = _hg.communityID(u);
    _community_history[community_id].emplace_back(_hg.contract(u, v, community_id));
    _pruner[community_id].removeSingleNodeHyperedges(_hg, community_id, _community_history[community_id].back());
    _pruner[community_id].removeParallelHyperedges(_hg, community_id, _community_history[community_id].back());
  }

  bool doUncoarsen() {
    ASSERT(!_init, "Community coarsener must be finalized before uncoarsening");

    while ( !_history.empty() ) {
      DBG << "Uncontracting: (" << _history.back().u << "," << _history.back().v << ")";
      PartitionID community_id = _hg.communityID(_history.back().u);
      _pruner[community_id].restoreParallelHyperedges(_hg, _history.back());
      _pruner[community_id].restoreSingleNodeHyperedges(_hg, _history.back());
      _hg.uncontract(_history.back());
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

 protected:
  HyperGraph& _hg;
  const Context& _context;
  bool _init;
  parallel::scalable_vector<HypergraphPruner> _pruner;
  parallel::scalable_vector<parallel::scalable_vector<Memento>> _community_history;
  std::vector<Memento> _history;
};

}  // namespace mt_kahypar
