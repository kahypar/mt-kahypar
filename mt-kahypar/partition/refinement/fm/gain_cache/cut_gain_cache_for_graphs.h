/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/meta/policy_registry.h"

#include "tbb/parallel_for.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/concurrent_vector.h"

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

// Forward
class DeltaGraphCutGainCache;

class GraphCutGainCache final : public kahypar::meta::PolicyBase {

 public:
  static constexpr FMGainCacheType TYPE = FMGainCacheType::cut_gain_cache_for_graphs;
  using DeltaGainCache = DeltaGraphCutGainCache;

  GraphCutGainCache() :
    _is_initialized(false),
    _k(kInvalidPartition),
    _gain_cache() { }

  GraphCutGainCache(const GraphCutGainCache&) = delete;
  GraphCutGainCache & operator= (const GraphCutGainCache &) = delete;

  GraphCutGainCache(GraphCutGainCache&& other) = default;
  GraphCutGainCache & operator= (GraphCutGainCache&& other) = default;

  // ####################### Initialization #######################

  bool isInitialized() const {
    return _is_initialized;
  }

  void reset(const bool run_parallel = true) {
    if ( _is_initialized ) {
      _gain_cache.assign(_gain_cache.size(),  CAtomic<HyperedgeWeight>(0), run_parallel);
    }
    _is_initialized = false;
  }

  size_t size() const {
    return _gain_cache.size();
  }

  // ! Initializes all gain cache entries
  template<typename PartitionedGraph>
  void initializeGainCache(const PartitionedGraph& partitioned_graph) {
    ASSERT(!_is_initialized, "Gain cache is already initialized");
    ASSERT(_k == kInvalidPartition || _k == partitioned_graph.k(), "Gain cache was already initialized for a different k");
    allocateGainTable(partitioned_graph.topLevelNumNodes(), partitioned_graph.k());

    // assert that current gain values are zero
    ASSERT(!_is_initialized &&
           std::none_of(_gain_cache.begin(), _gain_cache.end(),
             [&](const auto& weight) { return weight.load() != 0; }));

    // Initialize gain cache
    partitioned_graph.doParallelForAllEdges([&](const HyperedgeID e) {
      const HypernodeID u = partitioned_graph.edgeSource(e);
      if (partitioned_graph.nodeIsEnabled(u) && !partitioned_graph.isSinglePin(e)) {
        size_t index = incident_weight_index(u,
          partitioned_graph.partID(partitioned_graph.edgeTarget(e)));
        _gain_cache[index].fetch_add(partitioned_graph.edgeWeight(e), std::memory_order_relaxed);
      }
    });

    _is_initialized = true;
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := w(u, partID(u))
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID from) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[incident_weight_index(u, from)].load(std::memory_order_relaxed);
  }

  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputePenaltyTermEntry(const PartitionedGraph&,
                                 const HypernodeID) {
    // Do nothing here (only relevant for hypergraph gain cache)
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w(u, V_j)
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[incident_weight_index(u, to)].load(std::memory_order_relaxed);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u, const PartitionID from, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

  // ####################### Delta Gain Update #######################

  // ! This functions implements the delta gain updates for the cut metric on plain graphs.
  // ! When moving a node from its current block from to a target block to, we iterate
  // ! over its incident edges and syncronize the move on each edge. After syncronization,
  // ! we call this function to update the gain cache to changes associated with
  // ! corresponding edge.
  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedGraph& partitioned_graph,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID /* only relevant for hypergraphs */,
                       const PartitionID to,
                       const HypernodeID /* only relevant for hypergraphs */) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    const HypernodeID target = partitioned_graph.edgeTarget(he);
    const size_t index_in_from_part = incident_weight_index(target, from);
    _gain_cache[index_in_from_part].fetch_sub(we, std::memory_order_relaxed);
    const size_t index_in_to_part = incident_weight_index(target, to);
    _gain_cache[index_in_to_part].fetch_add(we, std::memory_order_relaxed);
  }

  // ####################### Uncontraction #######################

  // ! This function implements the gain cache update after an uncontraction that restores node v in
  // ! an edge he. After the uncontraction the corresponding edge turns from a selfloop to a regular edge.
  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterRestore(const PartitionedGraph& partitioned_graph,
                                    const HypernodeID u,
                                    const HypernodeID v,
                                    const HyperedgeID he,
                                    const HypernodeID) {
    if ( _is_initialized ) {
      // the edge weight is added to u and v
      const PartitionID block = partitioned_graph.partID(u);
      const HyperedgeWeight we = partitioned_graph.edgeWeight(he);
      _gain_cache[incident_weight_index(u, block)].fetch_add(we, std::memory_order_relaxed);
      _gain_cache[incident_weight_index(v, block)].fetch_add(we, std::memory_order_relaxed);
    }
  }

  // ! This function implements the gain cache update after an uncontraction that replaces u with v in
  // ! an edge he. After the uncontraction only node v is part of edge he.
  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterReplacement(const PartitionedGraph& partitioned_graph,
                                        const HypernodeID u,
                                        const HypernodeID v,
                                        const HyperedgeID he) {
    if ( _is_initialized ) {
      // the edge weight shifts from u to v
      const HypernodeID w = partitioned_graph.edgeTarget(he);
      const PartitionID block_of_w = partitioned_graph.partID(w);
      const HyperedgeWeight we = partitioned_graph.edgeWeight(he);
      _gain_cache[incident_weight_index(u, block_of_w)].fetch_sub(we, std::memory_order_relaxed);
      _gain_cache[incident_weight_index(v, block_of_w)].fetch_add(we, std::memory_order_relaxed);
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void restoreSinglePinHyperedge(const HypernodeID,
                                 const PartitionID,
                                 const HyperedgeWeight) {
    // Do nothing here (only relevant for hypergraph gain cache)
  }

  // ####################### Only for Testing #######################

  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedGraph& partitioned_graph,
                                       const HypernodeID u) const {
    PartitionID block_of_u = partitioned_graph.partID(u);
    HyperedgeWeight penalty = 0;
    for (HyperedgeID e : partitioned_graph.incidentEdges(u)) {
      if (!partitioned_graph.isSinglePin(e) &&
          partitioned_graph.partID(partitioned_graph.edgeTarget(e)) == block_of_u) {
        penalty += partitioned_graph.edgeWeight(e);
      }
    }
    return penalty;
  }

  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedGraph& partitioned_graph,
                                       const HypernodeID u,
                                       const PartitionID to) const {
    HyperedgeWeight benefit = 0;
    for (HyperedgeID e : partitioned_graph.incidentEdges(u)) {
      if (!partitioned_graph.isSinglePin(e) &&
          partitioned_graph.partID(partitioned_graph.edgeTarget(e)) == to) {
        benefit += partitioned_graph.edgeWeight(e);
      }
    }
    return benefit;
  }

 private:
  friend class DeltaGraphCutGainCache;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t incident_weight_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * _k  + p;
  }

  // ! Allocates the memory required to store the gain cache
  void allocateGainTable(const HypernodeID num_nodes,
                         const PartitionID k) {
    if (_gain_cache.size() == 0) {
      _k = k;
      _gain_cache.resize("Refinement", "incident_weight_in_part", num_nodes * size_t(_k), true);
    }
  }

  // ! Indicate whether or not the gain cache is initialized
  bool _is_initialized;

  // ! Number of blocks
  PartitionID _k;

  // ! The gain for moving a node u to from its current block V_i to a target block V_j
  // ! can be expressed as follows for the cut metric of plain graphs
  // ! g(u, V_j) := w(u, V_j) - w(u, V_i) = b(u, V_j) - p(u)
  // ! Here, w(u, V') is the weight of all edges connecting u to block V'.
  // ! We call b(u, V_j) the benefit term and p(u) the penalty term. Our gain cache stores and maintains these
  // ! entries for each node and block. Thus, the gain cache stores k entries per node (since b(u, V_i) = p(u)).
  ds::Array< CAtomic<HyperedgeWeight> > _gain_cache;
};

class DeltaGraphCutGainCache {

 public:
  DeltaGraphCutGainCache(const GraphCutGainCache& gain_cache) :
    _gain_cache(gain_cache),
    _incident_weight_in_part_delta() { }

  // ####################### Initialize & Reset #######################

  void initialize(const size_t size) {
    _incident_weight_in_part_delta.initialize(size);
  }

  void clear() {
    _incident_weight_in_part_delta.clear();
  }

  void dropMemory() {
    _incident_weight_in_part_delta.freeInternalData();
  }

  size_t size_in_bytes() const {
    return _incident_weight_in_part_delta.size_in_bytes();
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID from) const {
    const HyperedgeWeight* penalty_delta =
      _incident_weight_in_part_delta.get_if_contained(
        _gain_cache.incident_weight_index(u, from));
    return _gain_cache.penaltyTerm(u, from) + (penalty_delta ? *penalty_delta : 0);
  }

  // ! Returns the benefit term for moving node u to block to.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    const HyperedgeWeight* benefit_delta =
      _incident_weight_in_part_delta.get_if_contained(
        _gain_cache.incident_weight_index(u, to));
    return _gain_cache.benefitTerm(u, to) + (benefit_delta ? *benefit_delta : 0);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID from,
                       const PartitionID to ) const {
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedGraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedGraph& partitioned_graph,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID,
                       const PartitionID to,
                       const HypernodeID) {
    const HypernodeID target = partitioned_graph.edgeTarget(he);
    const size_t index_in_from_part = _gain_cache.incident_weight_index(target, from);
    _incident_weight_in_part_delta[index_in_from_part] -= we;
    const size_t index_in_to_part = _gain_cache.incident_weight_index(target, to);
    _incident_weight_in_part_delta[index_in_to_part] += we;
  }


 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    utils::MemoryTreeNode* gain_cache_delta_node = parent->addChild("Delta Gain Cache");
    gain_cache_delta_node->updateSize(size_in_bytes());
  }

 private:
  const GraphCutGainCache& _gain_cache;

  // ! Stores the delta of each locally touched gain cache entry
  // ! relative to the gain cache in '_phg'
  ds::DynamicFlatMap<size_t, HyperedgeWeight> _incident_weight_in_part_delta;
};

}  // namespace mt_kahypar
