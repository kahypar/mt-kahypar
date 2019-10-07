/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/utils/math.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"

namespace mt_kahypar {

template< typename TypeTraits >
class HypergraphPrunerT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;
  using Memento = typename StreamingHyperGraph::Memento;

 protected:
  static constexpr bool debug = false;

  static constexpr HyperedgeID kInvalidID = std::numeric_limits<HyperedgeID>::max();

  struct Fingerprint {
    HyperedgeID id;
    size_t hash;
  };

 public:
  explicit HypergraphPrunerT(ds::ThreadSafeFastResetFlagArray<>& contained_hypernodes,
                             parallel::scalable_vector<HyperedgeID>& parallel_he_representative) :
    _removed_single_node_hyperedges(),
    _removed_parallel_hyperedges(),
    _fingerprints(),
    _contained_hypernodes(contained_hypernodes),
    _parallel_he_representative(parallel_he_representative) { }

  HypergraphPrunerT(const HypergraphPrunerT&) = delete;
  HypergraphPrunerT& operator= (const HypergraphPrunerT&) = delete;

  HypergraphPrunerT(HypergraphPrunerT&& other) = default;
  HypergraphPrunerT& operator= (HypergraphPrunerT&&) = default;

  virtual ~HypergraphPrunerT() = default;

  void restoreSingleNodeHyperedges(HyperGraph& hypergraph,
                                   const Memento& memento) {
    for (int i = memento.one_pin_hes_begin + memento.one_pin_hes_size - 1;
         i >= memento.one_pin_hes_begin; --i) {
      ASSERT(i >= 0 && static_cast<size_t>(i) < _removed_single_node_hyperedges.size(),
             "Index out of bounds" << i);
      DBG << "restore single-node HE "
          << _removed_single_node_hyperedges[i];
      hypergraph.restoreSinglePinHyperedge(_removed_single_node_hyperedges[i]);
      _removed_single_node_hyperedges.pop_back();
    }
  }

  void restoreParallelHyperedges(HyperGraph& hypergraph,
                                 const Memento& memento) {
    for (int i = memento.parallel_hes_begin + memento.parallel_hes_size - 1;
         i >= memento.parallel_hes_begin; --i) {
      ASSERT(i >= 0 && static_cast<size_t>(i) < _removed_parallel_hyperedges.size(),
             "Index out of bounds:" << i);
      const HyperedgeID restore_parallel_edge = _removed_parallel_hyperedges[i];
      DBG << "restore parallel HE" << restore_parallel_edge;
      hypergraph.restoreParallelHyperedge(restore_parallel_edge, _parallel_he_representative);
      _removed_parallel_hyperedges.pop_back();
    }
  }

  HyperedgeWeight removeSingleNodeHyperedges(HyperGraph& hypergraph,
                                             const PartitionID community_id,
                                             Memento& memento) {
    memento.one_pin_hes_begin = _removed_single_node_hyperedges.size();
    auto begin_it = hypergraph.incidentEdges(memento.u, community_id).first;
    auto end_it = hypergraph.incidentEdges(memento.u, community_id).second;
    HyperedgeWeight removed_he_weight = 0;
    for (auto he_it = begin_it; he_it != end_it; ++he_it) {
      if (hypergraph.numCommunitiesOfHyperedge(*he_it) == 1 &&
          hypergraph.edgeSize(*he_it, community_id) == 1) {
        ASSERT(hypergraph.edgeIsEnabled(*he_it));
        _removed_single_node_hyperedges.push_back(*he_it);
        removed_he_weight += hypergraph.edgeWeight(*he_it);
        ++memento.one_pin_hes_size;
        DBG << "removing single-node HE" << *he_it;
        hypergraph.removeSinglePinEdge(*he_it, community_id);
        --he_it;
        --end_it;
      }
    }
    return removed_he_weight;
  }

  HyperedgeWeight removeParallelHyperedges(HyperGraph& hypergraph,
                                             const PartitionID community_id,
                                           Memento& memento) {
    memento.parallel_hes_begin = _removed_parallel_hyperedges.size();

    createFingerprints(hypergraph, memento.u, community_id);
    std::sort(_fingerprints.begin(), _fingerprints.end(),
              [](const Fingerprint& a, const Fingerprint& b) {
                return a.hash < b.hash || (a.hash == b.hash && a.id < b.id); });

    HyperedgeWeight removed_parallel_hes = 0;
    size_t i = 0;
    bool filled_probe_bitset = false;
    int cpu_id = sched_getcpu();
    while (i < _fingerprints.size()) {
      size_t j = i + 1;
      DBG << "i=" << i << ", j=" << j;
      if (_fingerprints[i].id != kInvalidID) {
        ASSERT(hypergraph.edgeIsEnabled(_fingerprints[i].id), V(_fingerprints[i].id));
        while (j < _fingerprints.size() && _fingerprints[i].hash == _fingerprints[j].hash) {
          // If we are here, then we have a hash collision for _fingerprints[i].id and
          // _fingerprints[j].id.
          DBG << _fingerprints[i].hash << "==" << _fingerprints[j].hash;
          DBG << "Size:" << hypergraph.edgeSize(_fingerprints[i].id) << "=="
              << hypergraph.edgeSize(_fingerprints[j].id);
          if (_fingerprints[j].id != kInvalidID &&
              hypergraph.numCommunitiesOfHyperedge(_fingerprints[i].id) ==
              hypergraph.numCommunitiesOfHyperedge(_fingerprints[j].id) &&
              hypergraph.edgeSize(_fingerprints[i].id) ==
              hypergraph.edgeSize(_fingerprints[j].id)) {
            ASSERT(hypergraph.edgeIsEnabled(_fingerprints[j].id), V(_fingerprints[j].id));
            if (!filled_probe_bitset) {
              fillProbeBitset(hypergraph, _fingerprints[i].id, cpu_id);
              filled_probe_bitset = true;
            }
            if (isParallelHyperedge(hypergraph, _fingerprints[j].id, cpu_id)) {
              removed_parallel_hes += 1;
              removeParallelHyperedge(hypergraph, _fingerprints[i].id, _fingerprints[j].id, community_id);
              _fingerprints[j].id = kInvalidID;
              ++memento.parallel_hes_size;
            }
          }
          ++j;
        }
      }
      // We need pairwise comparisons for all HEs with same hash.
      filled_probe_bitset = false;
      ++i;
    }
    return removed_parallel_hes;
  }

  const parallel::scalable_vector<HyperedgeID> & removedSingleNodeHyperedges() const {
    return _removed_single_node_hyperedges;
  }

  const parallel::scalable_vector<HyperedgeID> & removedParallelHyperedges() const {
    return _removed_parallel_hyperedges;
  }

 private:

  void createFingerprints(HyperGraph& hypergraph, const HypernodeID u, const PartitionID community_id) {
    _fingerprints.clear();
    for (const HyperedgeID& he : hypergraph.validIncidentEdges(u, community_id)) {
      DBG << "Fingerprint for HE" << he << "= {" << he << "," << hypergraph.edgeHash(he)
          << "," << hypergraph.edgeSize(he) << "}";
      _fingerprints.emplace_back(Fingerprint { he, hypergraph.edgeHash(he) });
    }
  }

  void fillProbeBitset(HyperGraph& hypergraph, const HyperedgeID he, const int cpu_id) {
    _contained_hypernodes.reset(cpu_id);
    DBG << "Filling Bitprobe Set for HE" << he;
    for ( const PartitionID& community_id : hypergraph.communities(he) ) {
      for (const HypernodeID& pin : hypergraph.pins(he, community_id)) {
        _contained_hypernodes.set(cpu_id, hypergraph.originalNodeID(pin), true);
      }
    }
  }


  bool isParallelHyperedge(HyperGraph& hypergraph, const HyperedgeID he, const int cpu_id) const {
    bool is_parallel = true;
    for ( const PartitionID& community_id : hypergraph.communities(he) ) {
      for (const HypernodeID& pin : hypergraph.pins(he, community_id)) {
        if (!_contained_hypernodes.get(cpu_id, hypergraph.originalNodeID(pin))) {
          is_parallel = false;
          break;
        }
      }
      if ( !is_parallel ) {
        break;
      }
    }
    DBG << "HE" << he << "is parallel HE=" << is_parallel;
    return is_parallel;
  }

  void removeParallelHyperedge(HyperGraph& hypergraph,
                               const HyperedgeID representative,
                               const HyperedgeID to_remove,
                               const PartitionID community_id) {
    hypergraph.setEdgeWeight(representative, community_id,
                             hypergraph.edgeWeight(representative, community_id)
                             + hypergraph.edgeWeight(to_remove, community_id));
    DBG << "removed HE" << to_remove << "which was parallel to" << representative;
    hypergraph.removeParallelEdge(to_remove, community_id);
    _removed_parallel_hyperedges.emplace_back(to_remove);
    ASSERT(hypergraph.originalEdgeID(to_remove) < _parallel_he_representative.size());
    _parallel_he_representative[hypergraph.originalEdgeID(to_remove)] = hypergraph.originalEdgeID(representative);
  }

  parallel::scalable_vector<HyperedgeID> _removed_single_node_hyperedges;
  parallel::scalable_vector<HyperedgeID> _removed_parallel_hyperedges;
  parallel::scalable_vector<Fingerprint> _fingerprints;
  ds::ThreadSafeFastResetFlagArray<>& _contained_hypernodes;
  parallel::scalable_vector<HyperedgeID>& _parallel_he_representative;
};

} // namespace mt_kahypar