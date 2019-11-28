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

#include <queue>
#include <string>

#include "tbb/parallel_for.h"
#include "tbb/task_group.h"
#include "tbb/parallel_invoke.h"
#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_vector.h"

#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/hypergraph_pruner.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
template <typename TypeTraits>
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
    _batch_hyperedges(hypergraph.initialNumEdges()),
    _incident_nets_of_v(hypergraph.initialNumEdges()),
    _batch_arena(std::max(_context.shared_memory.num_threads - 1, 1UL), 0),
    _batch_group(),
    _parallel_he_representative(hypergraph.initialNumEdges(), kInvalidHyperedge),
    _pruner(),
    _community_history(_hg.numCommunities()),
    _history() {
    for (PartitionID community_id = 0; community_id < _hg.numCommunities(); ++community_id) {
      _pruner.emplace_back(_contained_hypernodes, _parallel_he_representative);
    }
  }

  CommunityCoarsenerBase(const CommunityCoarsenerBase&) = delete;
  CommunityCoarsenerBase(CommunityCoarsenerBase&&) = delete;
  CommunityCoarsenerBase & operator= (const CommunityCoarsenerBase &) = delete;
  CommunityCoarsenerBase & operator= (CommunityCoarsenerBase &&) = delete;

  virtual ~CommunityCoarsenerBase() throw() { }

 protected:
  void init() {
    // Initialize community hyperedges such that parallel contractions in each community are possible
    utils::Timer::instance().start_timer("initialize_community_hyperedges", "Initialize Community Hyperedges");
    _hg.initializeCommunityHyperedges();
    _init = true;
    utils::Timer::instance().stop_timer("initialize_community_hyperedges");
  }

  void finalize() {
    // Merge Coarsening Mementos
    utils::Timer::instance().start_timer("merge_contractions", "Merge Contractions");
    mergeCommunityContractions();
    utils::Timer::instance().stop_timer("merge_contractions");

    // Reset community hyperedges
    utils::Timer::instance().start_timer("reset_community_hyperedges", "Reset Community Hyperedges");
    _hg.removeCommunityHyperedges(_history);
    utils::Timer::instance().stop_timer("reset_community_hyperedges");

    utils::Timer::instance().start_timer("postprocess_parallel_hyperedges", "Postprocess Parallel Hyperedges");
    postprocessParallelHyperedges();
    utils::Timer::instance().stop_timer("postprocess_parallel_hyperedges");

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

  bool doUncoarsen(std::unique_ptr<IRefiner>& label_propagation) {
    ASSERT(!_init, "Community coarsener must be finalized before uncoarsening");

    kahypar::Metrics current_metrics = { metrics::hyperedgeCut(_hg),
                                         metrics::km1(_hg),
                                         metrics::imbalance(_hg, _context) };
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);

    std::vector<HypernodeID> refinement_nodes;
    size_t max_batch_size = _context.refinement.use_batch_uncontractions &&
                            _context.shared_memory.num_threads > 1 ? _context.refinement.batch_size : 1;
    while (!_history.empty()) {
      // utils::Timer::instance().start_timer("uncontraction", "Uncontraction");
        std::vector<Memento> batch;
      while ( !_history.empty() && batch.size() < max_batch_size ) {
        batch.emplace_back(_history.back());
        _history.pop_back();
      }
      batchUncontraction(batch, refinement_nodes);
      _hg.updateGlobalPartInfos();
      // utils::Timer::instance().stop_timer("uncontraction");

      // Call label propagation refiner
      if (label_propagation) {
        // NOTE, label propagation refiner relies on the assumption, that it is called after
        // each uncontraction. Do not move the refiner call out of this loop.
        label_propagation->refine(refinement_nodes, current_metrics);
      }

      refinement_nodes.clear();
    }

    ASSERT(metrics::objective(_hg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_hg, _context.partition.objective)));
    return true;
  }

 private:
  void batchUncontraction(const std::vector<Memento>& batch,
                          std::vector<HypernodeID>& refinement_nodes) {
    ASSERT(!batch.empty());
    if ( batch.size() == 1 ) {
      uncontract(batch[0]);
      // NOTE, label propagation refiner relies on the assumption, that on the second position
      // of the refinement_nodes vector always the contraction partner occurs. Do not change the
      // order here.
      refinement_nodes.push_back(batch[0].u);
      refinement_nodes.push_back(batch[0].v);
      return;
    }

    tbb::concurrent_vector<tbb::concurrent_queue<Memento>> batch_queues;
    int64_t valid_batch_queues = 0;
    std::atomic<int64_t> current_batch_queue(-1);
    size_t current_batch_size = 1;
    std::atomic<size_t> num_processed(0);
    bool terminated = false;
    tbb::parallel_invoke([&] {
      // Master Thread

      auto switch_batch_queue = [&] {
        size_t tmp_num_processed = valid_batch_queues == 1 && current_batch_queue == -1 ? 1 : num_processed.load();
        ASSERT( tmp_num_processed <= current_batch_size );
        if ( tmp_num_processed == current_batch_size &&
             current_batch_queue.load() + 1 < valid_batch_queues ) {
          ASSERT( current_batch_queue < 0 || batch_queues[current_batch_queue].empty() );
          num_processed = 0;
          current_batch_size = batch_queues[current_batch_queue + 1].unsafe_size();
          ++current_batch_queue;
        } else if ( tmp_num_processed == current_batch_size && terminated ) {
          ASSERT( current_batch_queue.load() + 1 == valid_batch_queues );
          ++current_batch_queue;
        }
      };

      kahypar::ds::FastResetFlagArray<>& batch_hypernodes = _contained_hypernodes.local();
      size_t current_batch_idx = 0;
      while ( current_batch_idx < batch.size() ) {
        batch_hypernodes.reset();
        _batch_hyperedges.reset();
        batch_queues.emplace_back();
        for ( ; current_batch_idx < batch.size(); ++current_batch_idx ) {
          _incident_nets_of_v.reset();
          const Memento& current_memento = batch[current_batch_idx];
          bool found_conflict = false;

          const HypernodeID original_u_id = _hg.originalNodeID(current_memento.u);
          const HypernodeID original_v_id = _hg.originalNodeID(current_memento.v);
          if ( batch_hypernodes[original_u_id] || batch_hypernodes[original_v_id] ) {
            found_conflict = true;
            break;
          }
          batch_hypernodes.set(original_u_id, true);
          batch_hypernodes.set(original_v_id, true);

          ASSERT(!_hg.nodeIsEnabled(current_memento.v));
          _hg.enableHypernode(current_memento.v);
          for ( const HyperedgeID& he : _hg.incidentEdges(current_memento.v) ) {
            const HyperedgeID original_id = _hg.originalEdgeID(he);
            if ( _batch_hyperedges[original_id] ) {
              found_conflict = true;
              break;
            }
            _incident_nets_of_v.set(original_id, true);
            _batch_hyperedges.set(original_id, true);
          }
          _hg.disableHypernode(current_memento.v);

          if ( found_conflict ) {
            break;
          }

          ASSERT(_hg.nodeIsEnabled(current_memento.u));
          for ( const HyperedgeID& he : _hg.incidentEdges(current_memento.u) ) {
            const HyperedgeID original_id = _hg.originalEdgeID(he);
            if ( !_incident_nets_of_v[original_id] && _batch_hyperedges[original_id] ) {
              found_conflict = true;
              break;
            }
            _batch_hyperedges.set(original_id, true);
          }

          if ( found_conflict ) {
            break;
          }

          batch_queues.back().push(batch[current_batch_idx]);
          refinement_nodes.push_back(current_memento.u);
          refinement_nodes.push_back(current_memento.v);
          switch_batch_queue();
        }
        ++valid_batch_queues;
        switch_batch_queue();
      }
      terminated = true;
      while ( current_batch_queue < valid_batch_queues ) {
        switch_batch_queue();
      }
    }, [&] {
      // Worker Thread
      ASSERT(_context.shared_memory.num_threads > 1);
      size_t num_threads = _context.shared_memory.num_threads - 1;
      _batch_arena.execute([&] {
        for ( size_t i = 0; i < num_threads; ++i ) {
          _batch_group.run([&] {

            Memento current_memento;
            while ( !terminated || current_batch_queue < valid_batch_queues ) {

              int64_t tmp_current_batch_queue = current_batch_queue.load();
              if ( tmp_current_batch_queue >= 0 &&
                   tmp_current_batch_queue < valid_batch_queues &&
                   batch_queues[tmp_current_batch_queue].try_pop(current_memento) ) {
                // Uncontract
                uncontract(current_memento);
                ++num_processed;
              }

            }

          });
        }
      });
      _batch_arena.execute([&] { _batch_group.wait(); });
    });
  }

  void uncontract(const Memento& memento) {
    PartitionID community_id = _hg.communityID(memento.u);
    DBG << "Uncontracting: (" << memento.u << "," << memento.v << ")";
    _pruner[community_id].restoreParallelHyperedges(_hg, memento);
    _pruner[community_id].restoreSingleNodeHyperedges(_hg, memento);
    _hg.uncontract(memento, _parallel_he_representative);
  }

  void mergeCommunityContractions() {
    size_t size = 0;
    for (const auto& community_mementos : _community_history) {
      size += community_mementos.size();
    }
    std::sort(_community_history.begin(), _community_history.end(),
              [&](const auto& lhs, const auto& rhs) {
        return lhs.size() > rhs.size();
      });
    while (_community_history.back().size() == 0) {
      _community_history.pop_back();
    }

    // Mementos from all communities are written in parallel
    // to one history vector in round-robin fashion.
    _history.resize(size);
    tbb::parallel_for(0UL, _community_history.size(), [&](const size_t& i) {
        size_t current_pos = i;
        int64_t k = _community_history.size();
        for (size_t pos = 0; pos < _community_history[i].size(); ++pos) {
          while (pos >= _community_history[k - 1].size()) {
            --k;
            ASSERT(k > 0);
          }

          ASSERT(current_pos < _history.size(), V(i) << V(current_pos) << V(_history.size()));
          ASSERT(_history[current_pos].u == std::numeric_limits<HypernodeID>::max(), V(current_pos) << V(i) << V(k));
          _history[current_pos] = std::move(_community_history[i][pos]);
          current_pos += k;
        }
      });

    // Release memory of community histories
    for (auto& community_mementos : _community_history) {
      parallel::scalable_vector<Memento> tmp;
      community_mementos = std::move(tmp);
    }

    ASSERT([&] {
        for (size_t i = 0; i < _history.size(); ++i) {
          if (_history[i].u == std::numeric_limits<HypernodeID>::max()) {
            return false;
          }
        }
        return true;
      } (), "Merging contraction hierarchies failed");
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
    utils::Timer::instance().start_timer("determine_he_weights", "Determine HE weights");
    std::vector<HyperedgeID> in_degree(_hg.initialNumEdges(), 0);
    for (HyperedgeID he = 0; he < _hg.initialNumEdges(); ++he) {
      HyperedgeID representative = _parallel_he_representative[he];
      if (representative != kInvalidHyperedge) {
        ASSERT(representative < in_degree.size());
        ++in_degree[representative];
      }
    }

    // Push all leaves into a queue
    std::queue<HyperedgeID> q;
    for (HyperedgeID he = 0; he < _hg.initialNumEdges(); ++he) {
      if (in_degree[he] == 0 && _parallel_he_representative[he] != kInvalidHyperedge) {
        q.push(he);
      }
    }

    // Perform a bottom up traversal of the forest and compute weight
    // of each hyperedge
    int64_t num_removed_parallel_hyperedges = 0;
    std::vector<HyperedgeWeight> weights(_hg.initialNumEdges(), 0);
    while (!q.empty()) {
      HyperedgeID original_he = q.front();
      HyperedgeID he = _hg.globalEdgeID(original_he);

      bool was_enabled = false;
      if (!_hg.edgeIsEnabled(he)) {
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
      if (representative != kInvalidHyperedge) {
        ++num_removed_parallel_hyperedges;
        ASSERT(representative < weights.size());
        weights[representative] += _hg.edgeWeight(he);
        --in_degree[representative];
        if (in_degree[representative] == 0) {
          q.push(representative);
        }
      }

      // If the hyperedge was temporary enabled (to set edge weight) or
      // the hyperedge is an inner node of the parallel hyperedge forest
      // we disable the hyperedge
      if (was_enabled || representative != kInvalidHyperedge) {
        ASSERT(_hg.edgeIsEnabled(he));
        _hg.disableHyperedge(he);
      }
    }
    utils::Timer::instance().stop_timer("determine_he_weights");

    utils::Timer::instance().start_timer("remove_disabled_hyperedges_from_incident_nets", "Remove Disabled HE from HNs");
    _hg.invalidateDisabledHyperedgesFromIncidentNets();
    utils::Timer::instance().stop_timer("remove_disabled_hyperedges_from_incident_nets");
  }

 protected:
  HyperGraph& _hg;
  const Context& _context;
  bool _init;
  ThreadLocalFastResetFlagArray _contained_hypernodes;
  kahypar::ds::FastResetFlagArray<> _batch_hyperedges;
  kahypar::ds::FastResetFlagArray<> _incident_nets_of_v;
  tbb::task_arena _batch_arena;
  tbb::task_group _batch_group;
  parallel::scalable_vector<HyperedgeID> _parallel_he_representative;
  parallel::scalable_vector<HypergraphPruner> _pruner;
  parallel::scalable_vector<parallel::scalable_vector<Memento> > _community_history;
  std::vector<Memento> _history;
};

template <typename TypeTraits>
HyperedgeID CommunityCoarsenerBase<TypeTraits>::kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
}  // namespace mt_kahypar
