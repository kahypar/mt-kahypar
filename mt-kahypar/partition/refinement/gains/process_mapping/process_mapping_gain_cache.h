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

#include <algorithm>

#include "kahypar/meta/policy_registry.h"

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/datastructures/static_bitset.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {

/**
 * The gain cache stores the gain values for all possible node moves for the process mapping metric.
 *
 * The process mapping problem asks for a mapping Π: V -> V_p of the node set V of a weighted hypergraph H = (V,E,c,w)
 * onto a process graph P = (V_P, E_P) such that the following objective function is minimized:
 * process_mapping(H, P, Π) := sum_{e \in E} dist_P(Λ(e)) * w(e)
 * Here, dist_P(Λ(e)) is shortest connections between all blocks Λ(e) contained in a hyperedge e using only edges
 * of the process graph. Computing dist_P(Λ(e)) reduces to the steiner tree problem which is an NP-hard problem.
 * However, we precompute all steiner trees up to a certain size and for larger connectivity sets Λ(e), we compute
 * a 2-approximation.
 *
 * The gain of moving a node u from its current block V_i to a target block V_j can be expressed as follows:
 * g(u,V_j) := sum_{e \in I(u): Φ(e,V_i) = 1 and Φ(e, V_j) > 0} Δdist_P(e, Λ(e)\{V_i}) * w(e) +
 *             sum_{e \in I(u): Φ(e,V_i) = 1 and Φ(e, V_j) = 0} Δdist_P(e, Λ(e)\{V_i} u {V_j}) * w(e) +
 *             sum_{e \in I(u): Φ(e,V_i) > 1 and Φ(e, V_j) = 0} Δdist_P(e, Λ(e) u {V_j}) * w(e)
 * For a set of blocks A, we define Δdist_P(e, A) := (dist_P(Λ(e)) - dist_P(A)). Moreover, Φ(e,V') is the number
 * of pins contained in hyperedge e which are also part of block V'. More formally, Φ(e,V') := |e n V'|.
 *
 * This gain cache implementation maintains the gain values g(u,V_j) for all nodes and their adjacent blocks.
 * Thus, the gain cache stores and maintains at most k entries per node where k := |V_P|.
*/
class ProcessMappingGainCache {

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

  using AdjacentBlocksIterator = IteratorRange<typename ds::ConnectivitySets::Iterator>;

 public:
  struct HyperedgeState {
    HyperedgeState() :
      version(0),
      update_version(0) { }

    CAtomic<uint32_t> version;
    CAtomic<uint32_t> update_version;
  };

  static constexpr GainPolicy TYPE = GainPolicy::process_mapping;
  static constexpr bool requires_notification_before_update = true;
  static constexpr bool initializes_gain_cache_entry_after_batch_uncontractions = true;

  ProcessMappingGainCache() :
    _is_initialized(false),
    _k(kInvalidPartition),
    _gain_cache(),
    _ets_benefit_aggregator([&] { return initializeBenefitAggregator(); }),
    _num_incident_edges_of_block(),
    _adjacent_blocks(),
    _version(),
    _ets_version() { }

  ProcessMappingGainCache(const ProcessMappingGainCache&) = delete;
  ProcessMappingGainCache & operator= (const ProcessMappingGainCache &) = delete;

  ProcessMappingGainCache(ProcessMappingGainCache&& other) = default;
  ProcessMappingGainCache & operator= (ProcessMappingGainCache&& other) = default;

  // ####################### Initialization #######################

  bool isInitialized() const {
    return _is_initialized;
  }

  void reset(const bool run_parallel = true) {
    unused(run_parallel);
    _is_initialized = false;
  }

  size_t size() const {
    return _gain_cache.size();
  }

  // ! Initializes all gain cache entries
  template<typename PartitionedHypergraph>
  void initializeGainCache(const PartitionedHypergraph& partitioned_hg);

  // ! Initializes the gain cache entry for a node
  template<typename PartitionedHypergraph>
  void initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID hn);

  // ! Returns an iterator over the adjacent blocks of a node
  AdjacentBlocksIterator adjacentBlocks(const HypernodeID hn) {
    return _adjacent_blocks.connectivitySet(hn);
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! Note that the process mapping gain cache does not maintain a
  // ! penalty term and returns zero in this case.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID,
                              const PartitionID) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return 0;
  }

  // ! Recomputes all gain cache entries for node u
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputeInvalidTerms(const PartitionedHypergraph& partitioned_hg,
                             const HypernodeID u) {
    vec<HyperedgeWeight>& benefit_aggregator = _ets_benefit_aggregator.local();
    initializeGainCacheEntryForNode(partitioned_hg, u, benefit_aggregator);
  }

  // ! Returns the gain value for moving node u to block to.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[benefit_index(u, to)].load(std::memory_order_relaxed);
  }

  // ! Returns the gain value for moving node u to block to.
  // ! (same as benefitTerm(...))
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID, /* only relevant for graphs */
                       const PartitionID to ) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return benefitTerm(u, to);
  }

  // ####################### Delta Gain Update #######################

  // ! This function returns true if the corresponding pin count values triggers
  // ! a gain cache update.
  static bool triggersDeltaGainUpdate(const SyncronizedEdgeUpdate& sync_update);

  // ! The partitioned (hyper)graph call this function when updating the pin count values
  // ! and connectivity set. Updating these data structures and calling this function are done
  // ! within one transaction (protected via a spin-lock).
  void updateVersionOfHyperedge(const SyncronizedEdgeUpdate& sync_update);

  // ! This functions implements the delta gain updates for the connecitivity metric.
  // ! When moving a node from its current block from to a target block to, we iterate
  // ! over its incident hyperedges and update their pin count values. After each pin count
  // ! update, we call this function to update the gain cache to changes associated with
  // ! corresponding hyperedge.
  template<typename PartitionedHypergraph>
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const SyncronizedEdgeUpdate& sync_update);

  // ####################### Uncontraction #######################

  // ! This function implements the gain cache update after an uncontraction that restores node v in
  // ! hyperedge he. After the uncontraction operation, node u and v are contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                    const HypernodeID u,
                                    const HypernodeID v,
                                    const HyperedgeID he,
                                    const HypernodeID pin_count_in_part_after);

  // ! This function implements the gain cache update after an uncontraction that replaces u with v in
  // ! hyperedge he. After the uncontraction only node v is contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID u,
                                        const HypernodeID v,
                                        const HyperedgeID he);

  // ! This function is called after restoring a single-pin hyperedge. The function assumes that
  // ! u is the only pin of the corresponding hyperedge, while block_of_u is its corresponding block ID.
  void restoreSinglePinHyperedge(const HypernodeID u,
                                 const PartitionID block_of_u,
                                 const HyperedgeWeight weight_of_he);

  // ! This function is called after restoring a net that became identical to another due to a contraction.
  template<typename PartitionedHypergraph>
  void restoreIdenticalHyperedge(const PartitionedHypergraph&,
                                 const HyperedgeID);

  // ####################### Only for Testing #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedHypergraph&,
                                       const HypernodeID) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return 0;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       const PartitionID to) const {
    ASSERT(partitioned_hg.hasProcessGraph());
    HyperedgeWeight gain = 0;
    const ProcessGraph& process_graph = *partitioned_hg.processGraph();
    const PartitionID from = partitioned_hg.partID(u);
    for (HyperedgeID e : partitioned_hg.incidentEdges(u)) {
      ds::Bitset& connectivity_set = partitioned_hg.deepCopyOfConnectivitySet(e);
      const HyperedgeWeight current_distance = process_graph.distance(connectivity_set);
      if ( partitioned_hg.pinCountInPart(e, from) == 1 ) {
        // Moving the node out of its current block removes
        // its block from the connectivity set
        connectivity_set.unset(from);
      }
      const HyperedgeWeight distance_with_to = process_graph.distanceWithBlock(connectivity_set, to);
      gain += (current_distance - distance_with_to) * partitioned_hg.edgeWeight(e);
    }
    return gain;
  }

 private:
  friend class DeltaProcessMappingGainCache;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t benefit_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * _k + p;
  }

  // ! Allocates the memory required to store the gain cache
  void allocateGainTable(const HypernodeID num_nodes,
                         const HyperedgeID num_edges,
                         const PartitionID k) {
    if (_gain_cache.size() == 0) {
      ASSERT(_k == kInvalidPartition);
      _k = k;
      tbb::parallel_invoke([&] {
        _gain_cache.resize(
          "Refinement", "gain_cache", num_nodes * _k, true);
      }, [&] {
        _num_incident_edges_of_block.resize(
          "Refinement", "num_incident_edges_of_block", num_nodes * _k, true);
      }, [&] {
        _adjacent_blocks = ds::ConnectivitySets(num_nodes, k, true);
      }, [&] {
        _version.assign(num_edges, HyperedgeState());
      });
    }
  }

  // ! Initializes the adjacent blocks of all nodes
  template<typename PartitionedHypergraph>
  void initializeAdjacentBlocks(const PartitionedHypergraph& partitioned_hg);

    // ! Initializes the adjacent blocks of for a node
  template<typename PartitionedHypergraph>
  void initializeAdjacentBlocksOfNode(const PartitionedHypergraph& partitioned_hg,
                                      const HypernodeID hn);

  // ! Updates the adjacent blocks of a node based on a synronized hyperedge update
  template<typename PartitionedHypergraph>
  void updateAdjacentBlocks(const PartitionedHypergraph& partitioned_hg,
                            const SyncronizedEdgeUpdate& sync_update);

  // ! Increments the number of incident edges of node u that contains pins of block to.
  // ! If the value increases to one, we add the block to the connectivity set of the node
  // ! u and initialize the gain cache entry for moving u to that block.
  HyperedgeID incrementIncidentEdges(const HypernodeID u, const PartitionID to);

  // ! Decrements the number of incident edges of node u that contains pins of block to
  // ! If the value decreases to zero, we remove the block from the connectivity set of the node.
  HyperedgeID decrementIncidentEdges(const HypernodeID u, const PartitionID to);

  // ! Initializes the benefit and penalty terms for a node u
  template<typename PartitionedHypergraph>
  void initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       vec<Gain>& benefit_aggregator);

  // ! Initializes the gain cache entry of moving u to block 'to'. The function is
  // ! thread-safe, meaning that it supports correct initialization while simultanously
  // ! performing gain cache updates.
  template<typename PartitionedHypergraph>
  void initializeGainCacheEntry(const PartitionedHypergraph& partitioned_hg,
                                const HypernodeID hn,
                                const PartitionID to,
                                ds::Array<SpinLock>& edge_locks);

  bool nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    if ( p == kInvalidPartition || p >= _k ) {
      LOG << "Invalid block ID (Node" << u << "is part of block" << p
          << ", but valid block IDs must be in the range [ 0," << _k << "])";
      return false;
    }
    if ( benefit_index(u, p) >= _gain_cache.size() ) {
      LOG << "Access to gain cache would result in an out-of-bounds access ("
          << "Benefit Index =" << benefit_index(u, p)
          << ", Gain Cache Size =" << _gain_cache.size() << ")";
      return false;
    }
    return true;
  }

  vec<Gain> initializeBenefitAggregator() const {
    return vec<Gain>(_k, std::numeric_limits<Gain>::min());
  }

  // ! Indicate whether or not the gain cache is initialized
  bool _is_initialized;

  // ! Number of blocks
  PartitionID _k;

  // ! Array of size |V| * (k + 1), which stores the benefit and penalty terms of each node.
  ds::Array< CAtomic<HyperedgeWeight> > _gain_cache;

  // ! Thread-local for initializing gain cache entries
  tbb::enumerable_thread_specific<vec<Gain>> _ets_benefit_aggregator;

  // ! This array stores the number of incident hyperedges that contains
  // ! pins of a particular block for each node.
  ds::Array< CAtomic<HyperedgeID> > _num_incident_edges_of_block;

  // ! Stores the adjacent blocks of a node
  ds::ConnectivitySets _adjacent_blocks;

  // ! This array stores a version ID for each hyperedge. The partitioned hypergraph
  // ! increments the version for a hyperedge before it updates it internal data structure
  // ! (see updateVersionOfHyperedge(...)). This can be use when initialize a new gain cache entries, while
  // ! other threads perform concurrent moves on the data structure.
  vec<HyperedgeState> _version;

  // ! Array to store version IDs when we lazily initialize a gain cache entry
  tbb::enumerable_thread_specific<vec<uint32_t>> _ets_version;
};

/**
 * In our FM algorithm, the different local searches perform nodes moves locally not visible for other
 * threads. The delta gain cache stores these local changes relative to the shared
 * gain cache. For example, the penalty term can be computed as follows
 * p'(u) := p(u) + Δp(u)
 * where p(u) is the penalty term stored in the shared gain cache and Δp(u) is the penalty term stored in
 * the delta gain cache after performing some moves locally. To maintain Δp(u) and Δb(u,V_j), we use a hash
 * table that only stores entries affected by a gain cache update.
*/
class DeltaProcessMappingGainCache {

 public:
  DeltaProcessMappingGainCache(const ProcessMappingGainCache& gain_cache) :
    _gain_cache(gain_cache),
    _gain_cache_delta() { }

  // ####################### Initialize & Reset #######################

  void initialize(const size_t size) {
    _gain_cache_delta.initialize(size);
  }

  void clear() {
    _gain_cache_delta.clear();
  }

  void dropMemory() {
    _gain_cache_delta.freeInternalData();
  }

  size_t size_in_bytes() const {
    return _gain_cache_delta.size_in_bytes();
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := w({ e \in I(u) | pin_count(e, V_i) > 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID,
                              const PartitionID) const {
    return 0;
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) >= 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(to != kInvalidPartition && to < _gain_cache._k);
    const HyperedgeWeight* benefit_delta =
      _gain_cache_delta.get_if_contained(_gain_cache.benefit_index(u, to));
    return _gain_cache.benefitTerm(u, to) + ( benefit_delta ? *benefit_delta : 0 );
  }

  // ! Returns the gain of moving node u from its current block to a target block V_j.
  // ! More formally, g(u, V_j) := b(u, V_j) - p(u).
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID from,
                       const PartitionID to ) const {
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph&,
                       const SyncronizedEdgeUpdate&) {

  }

 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    utils::MemoryTreeNode* gain_cache_delta_node = parent->addChild("Delta Gain Cache");
    gain_cache_delta_node->updateSize(size_in_bytes());
  }

 private:
  const ProcessMappingGainCache& _gain_cache;

  // ! Stores the delta of each locally touched gain cache entry
  // ! relative to the gain cache in '_phg'
  ds::DynamicFlatMap<size_t, HyperedgeWeight> _gain_cache_delta;
};

}  // namespace mt_kahypar
