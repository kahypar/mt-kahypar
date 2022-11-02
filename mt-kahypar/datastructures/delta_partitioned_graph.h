/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include <atomic>
#include <type_traits>

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Special graph data structure for the FM algorithm.
 * This is a variant of DeltaPartitionedHypergraph specialized for graphs.
 * See delte_partitioned_hypergraph.h for more details.
 */
template <typename PartitionedGraph = Mandatory>
class DeltaPartitionedGraph {
 private:
  static constexpr size_t MAP_SIZE_LARGE = 16384;
  static constexpr size_t MAP_SIZE_MOVE_DELTA = 8192;
  static constexpr size_t MAP_SIZE_SMALL = 128;

  using HypernodeIterator = typename PartitionedGraph::HypernodeIterator;
  using HyperedgeIterator = typename PartitionedGraph::HyperedgeIterator;
  using IncidenceIterator = typename PartitionedGraph::IncidenceIterator;
  using IncidentNetsIterator = typename PartitionedGraph::IncidentNetsIterator;

 public:
  static constexpr bool supports_connectivity_set = false;
  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = PartitionedGraph::HIGH_DEGREE_THRESHOLD;

  DeltaPartitionedGraph() :
    _k(kInvalidPartition),
    _pg(nullptr),
    _part_weights_delta(0, 0),
    _part_ids_delta(),
    _incident_weight_in_part_delta() {}

  DeltaPartitionedGraph(const Context& context) :
    _k(context.partition.k),
    _pg(nullptr),
    _part_weights_delta(context.partition.k, 0),
    _part_ids_delta(),
    _incident_weight_in_part_delta() {
      const bool top_level = context.type == ContextType::main;
      _part_ids_delta.initialize(MAP_SIZE_SMALL);
      _incident_weight_in_part_delta.initialize(top_level ? MAP_SIZE_LARGE : MAP_SIZE_MOVE_DELTA);
    }

  DeltaPartitionedGraph(const DeltaPartitionedGraph&) = delete;
  DeltaPartitionedGraph & operator= (const DeltaPartitionedGraph &) = delete;

  DeltaPartitionedGraph(DeltaPartitionedGraph&& other) = default;
  DeltaPartitionedGraph & operator= (DeltaPartitionedGraph&& other) = default;

  ~DeltaPartitionedGraph() = default;

  void setPartitionedHypergraph(PartitionedGraph* pg) {
    _pg = pg;
  }

  // ####################### Iterators #######################

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    ASSERT(_pg);
    return _pg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    ASSERT(_pg);
    return _pg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(_pg);
    return _pg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(_pg);
    return _pg->pins(e);
  }

  // ####################### Hypernode Information #######################

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(_pg);
    return _pg->nodeWeight(u);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(_pg);
    return _pg->nodeDegree(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Target of an edge
  HypernodeID edgeTarget(const HyperedgeID e) const {
    return _pg->edgeTarget(e);
  }

  // ! Source of an edge
  HypernodeID edgeSource(const HyperedgeID e) const {
    return _pg->edgeSource(e);
  }

  // ! Number of pins of an edge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(_pg);
    return _pg->edgeSize(e);
  }

  HyperedgeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(_pg);
    return _pg->edgeWeight(e);
  }

  // ####################### Partition Information #######################

  // ! Changes the block of hypernode u from 'from' to 'to'.
  // ! Move is successful, if it is not violating the balance
  // ! constraint specified by 'max_weight_to'.
  template<typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to,
                      DeltaFunc&& delta_func) {
    ASSERT(_pg);
    ASSERT(partID(u) == from);
    ASSERT(from != to);

    const HypernodeWeight weight = _pg->nodeWeight(u);
    if (partWeight(to) + weight <= max_weight_to) {
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += weight;
      _part_weights_delta[from] -= weight;

      for (const HyperedgeID edge : _pg->incidentEdges(u)) {
        const PartitionID target_part = partID(_pg->edgeTarget(edge));
        const HypernodeID pin_count_in_from_part_after = target_part == from ? 1 : 0;
        const HypernodeID pin_count_in_to_part_after = target_part == to ? 2 : 1;
        delta_func(edge, _pg->edgeWeight(edge), _pg->edgeSize(edge), pin_count_in_from_part_after, pin_count_in_to_part_after);
      }
      return true;
    } else {
      return false;
    }
  }

  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to) {
    ASSERT(_pg);
    ASSERT(partID(u) == from);
    ASSERT(from != to);

    const HypernodeWeight weight = _pg->nodeWeight(u);
    if (partWeight(to) + weight <= max_weight_to) {
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += weight;
      _part_weights_delta[from] -= weight;
      return true;
    } else {
      return false;
    }
  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         const PartitionID from,
                                         const PartitionID to,
                                         const HypernodeWeight max_weight_to) {
    auto delta_gain_func = [&](HyperedgeID he, HyperedgeWeight edge_weight,
                               HypernodeID ,HypernodeID pcip_from, HypernodeID pcip_to) {
      gainCacheUpdate(he, edge_weight, from, pcip_from, to, pcip_to);
    };
    return changeNodePart(u, from, to, max_weight_to, delta_gain_func);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID he, const HyperedgeWeight we,
                       const PartitionID from, const HypernodeID /*pin_count_in_from_part_after*/,
                       const PartitionID to, const HypernodeID /*pin_count_in_to_part_after*/) {
    const HypernodeID target = _pg->edgeTarget(he);
    const size_t index_in_from_part = incident_weight_index(target, from);
    _incident_weight_in_part_delta[index_in_from_part] -= we;
    const size_t index_in_to_part = incident_weight_index(target, to);
    _incident_weight_in_part_delta[index_in_to_part] += we;
  }

  // ! Returns the block of hypernode u
  PartitionID partID(const HypernodeID u) const {
    ASSERT(_pg);
    const PartitionID* part_id = _part_ids_delta.get_if_contained(u);
    return part_id ? *part_id : _pg->partID(u);
  }

  // ! Returns the total weight of block p
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(_pg);
    ASSERT(p != kInvalidPartition && p < _k);
    return _pg->partWeight(p) + _part_weights_delta[p];
  }

  // ! Returns the number of pins of edge e in block p
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(_pg);
    ASSERT(e < _pg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_pg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    HypernodeID count = 0;
    if (p == partID(edgeSource(e))) {
      count++;
    }
    if (!_pg->isSinglePin(e) && p == partID(edgeTarget(e))) {
      count++;
    }
    return count;
  }

  HyperedgeWeight moveFromPenalty(const HypernodeID u) const {
    ASSERT(_pg);
    const PartitionID part_id = partID(u);
    const HyperedgeWeight* incident_weight_delta_u =
      _incident_weight_in_part_delta.get_if_contained(incident_weight_index(u, part_id));
    const HyperedgeWeight incident_weight_u = _pg->incidentWeightInPart(u, part_id) +
                                              (incident_weight_delta_u ? *incident_weight_delta_u : 0);
    return incident_weight_u;
  }

  HyperedgeWeight moveToBenefit(const HypernodeID u, const PartitionID p) const {
    ASSERT(_pg);
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight* incident_weight_delta_p =
      _incident_weight_in_part_delta.get_if_contained(incident_weight_index(u, p));
    const HyperedgeWeight incident_weight_p = _pg->incidentWeightInPart(u, p) +
                                              (incident_weight_delta_p ? *incident_weight_delta_p : 0);
    return incident_weight_p;
  }

  Gain km1Gain(const HypernodeID u, const PartitionID from, const PartitionID to) const {
    unused(from);
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveToBenefit(u, to) - moveFromPenalty(u);
  }

  void initializeGainCacheEntry(const HypernodeID u, vec<Gain>& benefit_aggregator) {
    _pg->initializeGainCacheEntry(u, benefit_aggregator);
  }

  // ! Clears all deltas applied to the partitioned hypergraph
  void clear() {
    // O(k)
    _part_weights_delta.assign(_k, 0);
    // Constant Time
    _part_ids_delta.clear();
    _incident_weight_in_part_delta.clear();
  }

  void dropMemory() {
    if (!_memory_dropped) {
      _memory_dropped = true;
      _part_ids_delta.freeInternalData();
      _incident_weight_in_part_delta.freeInternalData();
    }
  }

  size_t combinedMemoryConsumption() const {
    return _part_ids_delta.size_in_bytes()
           + _incident_weight_in_part_delta.size_in_bytes();
  }

  PartitionID k() const {
    return _k;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* delta_pg_node = parent->addChild("Delta Partitioned Hypergraph");
    utils::MemoryTreeNode* part_weights_node = delta_pg_node->addChild("Delta Part Weights");
    part_weights_node->updateSize(_part_weights_delta.capacity() * sizeof(HypernodeWeight));
    utils::MemoryTreeNode* part_ids_node = delta_pg_node->addChild("Delta Part IDs");
    part_ids_node->updateSize(_part_ids_delta.size_in_bytes());
    utils::MemoryTreeNode* _incident_weight_in_part_node = delta_pg_node->addChild("Delta Incident Weight In Part");
    _incident_weight_in_part_node->updateSize(_incident_weight_in_part_delta.size_in_bytes());
  }

 private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t incident_weight_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * _k  + p;
  }

  bool _memory_dropped = false;

  // ! Number of blocks
  PartitionID _k;

  // ! Partitioned graph where all deltas are stored relative to
  PartitionedGraph* _pg;

  // ! Delta for block weights
  vec< HypernodeWeight > _part_weights_delta;

  // ! Stores for each locally moved node its new block id
  DynamicFlatMap<HypernodeID, PartitionID> _part_ids_delta;

  // ! Stores the delta of each locally touched incident weight in part entry
  // ! relative to the _incident_weight_in_part member in '_pg'
  DynamicFlatMap<size_t, HyperedgeWeight> _incident_weight_in_part_delta;
};

} // namespace ds
} // namespace mt_kahypar