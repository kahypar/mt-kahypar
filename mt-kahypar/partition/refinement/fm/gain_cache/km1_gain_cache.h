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
class DeltaKm1GainCache;

class Km1GainCache final : public kahypar::meta::PolicyBase {

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

 public:
  static constexpr FMGainCacheType TYPE = FMGainCacheType::km1_gain_cache;
  using DeltaGainCache = DeltaKm1GainCache;

  Km1GainCache() :
    _is_initialized(false),
    _k(kInvalidPartition),
    _gain_cache() { }

  Km1GainCache(const Km1GainCache&) = delete;
  Km1GainCache & operator= (const Km1GainCache &) = delete;

  Km1GainCache(Km1GainCache&& other) = default;
  Km1GainCache & operator= (Km1GainCache&& other) = default;

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
  void initializeGainCache(const PartitionedHypergraph& partitioned_hg) {
    ASSERT(!_is_initialized, "Gain cache is already initialized");
    ASSERT(_k == kInvalidPartition || _k == partitioned_hg.k(), "Gain cache was already initialized for a different k");
    allocateGainTable(partitioned_hg.topLevelNumNodes(), partitioned_hg.k());


    // Gain calculation consist of two stages
    //  1. Compute gain of all low degree vertices
    //  2. Compute gain of all high degree vertices
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mtb(_k, 0);
    tbb::concurrent_vector<HypernodeID> high_degree_vertices;
    // Compute gain of all low degree vertices
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), partitioned_hg.initialNumNodes()),
      [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& benefit_aggregator = ets_mtb.local();
        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
          if ( partitioned_hg.nodeIsEnabled(u)) {
            if ( partitioned_hg.nodeDegree(u) <= HIGH_DEGREE_THRESHOLD) {
              initializeGainCacheEntryForNode(partitioned_hg, u, benefit_aggregator);
            } else {
              // Collect high degree vertices
              high_degree_vertices.push_back(u);
            }
          }
        }
      });

    auto aggregate_contribution_of_he_for_node =
      [&](const PartitionID block_of_u,
          const HyperedgeID he,
          HyperedgeWeight& penalty_aggregator,
          vec<HyperedgeWeight>& benefit_aggregator) {
      HyperedgeWeight edge_weight = partitioned_hg.edgeWeight(he);
      if (partitioned_hg.pinCountInPart(he, block_of_u) > 1) {
        penalty_aggregator += edge_weight;
      }

      for (const PartitionID block : partitioned_hg.connectivitySet(he)) {
        benefit_aggregator[block] += edge_weight;
      }
    };

    // Compute gain of all high degree vertices
    for ( const HypernodeID& u : high_degree_vertices ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> ets_mfp(0);
      const PartitionID from = partitioned_hg.partID(u);
      const HypernodeID degree_of_u = partitioned_hg.nodeDegree(u);
      tbb::parallel_for(tbb::blocked_range<HypernodeID>(ID(0), degree_of_u),
        [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& benefit_aggregator = ets_mtb.local();
        HyperedgeWeight& penalty_aggregator = ets_mfp.local();
        size_t current_pos = r.begin();
        for ( const HyperedgeID& he : partitioned_hg.incidentEdges(u, r.begin()) ) {
          aggregate_contribution_of_he_for_node(from, he,
            penalty_aggregator, benefit_aggregator);
          ++current_pos;
          if ( current_pos == r.end() ) {
            break;
          }
        }
      });

      // Aggregate thread locals to compute overall gain of the high degree vertex
      const HyperedgeWeight penalty_term = ets_mfp.combine(std::plus<HyperedgeWeight>());
      _gain_cache[penalty_index(u)].store(penalty_term, std::memory_order_relaxed);
      for (PartitionID p = 0; p < _k; ++p) {
        HyperedgeWeight move_to_benefit = 0;
        for ( auto& l_move_to_benefit : ets_mtb ) {
          move_to_benefit += l_move_to_benefit[p];
          l_move_to_benefit[p] = 0;
        }
        _gain_cache[benefit_index(u, p)].store(move_to_benefit, std::memory_order_relaxed);
      }
    }

    _is_initialized = true;
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, partID(u)) = 1 }))
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID /* only relevant for graphs */) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[penalty_index(u)].load(std::memory_order_relaxed);
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputePenaltyTermEntry(const PartitionedHypergraph& partitioned_hg,
                                 const HypernodeID u) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    _gain_cache[penalty_index(u)].store(recomputePenaltyTerm(
      partitioned_hg, u), std::memory_order_relaxed);
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) >= 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[benefit_index(u, to)].load(std::memory_order_relaxed);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID, /* only relevant for graphs */
                       const PartitionID to ) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return benefitTerm(u, to) - penaltyTerm(u, kInvalidPartition);
  }

  // ####################### Delta Gain Update #######################

  // ! This functions implements the delta gain updates for the connecitivity metric.
  // ! When moving a node from its current block from to a target block to, we iterate
  // ! over its incident hyperedges and update their pin count values. After each pin count
  // ! update, we call this function to update the gain cache to changes associated with
  // ! corresponding hyperedge.
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to,
                       const HypernodeID pin_count_in_to_part_after) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    if ( pin_count_in_from_part_after == 1 ) {
      for (const HypernodeID& u : partitioned_hg.pins(he)) {
        ASSERT(nodeGainAssertions(u, from));
        if (partitioned_hg.partID(u) == from) {
          _gain_cache[penalty_index(u)].fetch_sub(we, std::memory_order_relaxed);
        }
      }
    } else if (pin_count_in_from_part_after == 0) {
      for (const HypernodeID& u : partitioned_hg.pins(he)) {
        ASSERT(nodeGainAssertions(u, from));
        _gain_cache[benefit_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
      }
    }

    if (pin_count_in_to_part_after == 1) {
      for (const HypernodeID& u : partitioned_hg.pins(he)) {
        ASSERT(nodeGainAssertions(u, to));
        _gain_cache[benefit_index(u, to)].fetch_add(we, std::memory_order_relaxed);
      }
    } else if (pin_count_in_to_part_after == 2) {
      for (const HypernodeID& u : partitioned_hg.pins(he)) {
        ASSERT(nodeGainAssertions(u, to));
        if (partitioned_hg.partID(u) == to) {
          _gain_cache[penalty_index(u)].fetch_add(we, std::memory_order_relaxed);
        }
      }
    }
  }

  // ####################### Uncontraction #######################

  // ! This function implements the gain cache update after an uncontraction that restores node v in
  // ! hyperedge he. After the uncontraction node u and v are contained in hyperedge he.
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                    const HypernodeID u,
                                    const HypernodeID v,
                                    const HyperedgeID he,
                                    const HypernodeID pin_count_in_part_after) {
    if ( _is_initialized ) {
      // If u was the only pin of hyperedge he in its block before then moving out vertex u
      // of hyperedge he does not decrease the connectivity any more after the
      // uncontraction => p(u) += w(he)
      const PartitionID block = partitioned_hg.partID(u);
      const HyperedgeWeight edge_weight = partitioned_hg.edgeWeight(he);
      if ( pin_count_in_part_after == 2 ) {
        // u might be replaced by an other vertex in the batch
        // => search for other pin of the corresponding block and
        // add edge weight.
        for ( const HypernodeID& pin : partitioned_hg.pins(he) ) {
          if ( pin != v && partitioned_hg.partID(pin) == block ) {
            _gain_cache[penalty_index(pin)].add_fetch(edge_weight, std::memory_order_relaxed);
            break;
          }
        }
      }

      _gain_cache[penalty_index(v)].add_fetch(edge_weight, std::memory_order_relaxed);
      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the b(u, block) for vertex v by w(e)
      for ( const PartitionID block : partitioned_hg.connectivitySet(he) ) {
        _gain_cache[benefit_index(v, block)].add_fetch(
          edge_weight, std::memory_order_relaxed);
      }
    }
  }

  // ! This function implements the gain cache update after an uncontraction that replaces u with v in
  // ! hyperedge he. After the uncontraction only node v is contained in hyperedge he.
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID u,
                                        const HypernodeID v,
                                        const HyperedgeID he) {
    // In this case, u is replaced by v in hyperedge he
    // => Pin counts of hyperedge he does not change
    if ( _is_initialized ) {
      const PartitionID block = partitioned_hg.partID(u);
      const HyperedgeWeight edge_weight = partitioned_hg.edgeWeight(he);
      // Since u is no longer incident to hyperedge he its contribution for decreasing
      // the connectivity of he is shifted to vertex v
      if ( partitioned_hg.pinCountInPart(he, block) == 1 ) {
        _gain_cache[penalty_index(u)].add_fetch(edge_weight, std::memory_order_relaxed);
        _gain_cache[penalty_index(v)].sub_fetch(edge_weight, std::memory_order_relaxed);
      }

      _gain_cache[penalty_index(u)].sub_fetch(
        edge_weight, std::memory_order_relaxed);
      _gain_cache[penalty_index(v)].add_fetch(
        edge_weight, std::memory_order_relaxed);
      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the move_to_benefit for vertex v by w(e) and decrease
      // it for vertex u by w(e)
      for ( const PartitionID block : partitioned_hg.connectivitySet(he) ) {
        _gain_cache[benefit_index(u, block)].sub_fetch(
          edge_weight, std::memory_order_relaxed);
        _gain_cache[benefit_index(v, block)].add_fetch(
          edge_weight, std::memory_order_relaxed);
      }
    }
  }

  // ! This function is called after restoring a single-pin hyperedge. The function assumes that
  // ! u is the only pin of the corresponding hyperedge, while block_of_u is its corresponding block ID.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void restoreSinglePinHyperedge(const HypernodeID u,
                                 const PartitionID block_of_u,
                                 const HyperedgeWeight weight_of_he) {
    if ( _is_initialized ) {
      _gain_cache[benefit_index(u, block_of_u)].add_fetch(
        weight_of_he, std::memory_order_relaxed);
    }
  }

  // ####################### Only for Testing #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    const PartitionID block_of_u = partitioned_hg.partID(u);
    HyperedgeWeight penalty = 0;
    for (HyperedgeID e : partitioned_hg.incidentEdges(u)) {
      if ( partitioned_hg.pinCountInPart(e, block_of_u) > 1 ) {
        penalty += partitioned_hg.edgeWeight(e);
      }
    }
    return penalty;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       const PartitionID to) const {
    HyperedgeWeight benefit = 0;
    for (HyperedgeID e : partitioned_hg.incidentEdges(u)) {
      if (partitioned_hg.pinCountInPart(e, to) >= 1) {
        benefit += partitioned_hg.edgeWeight(e);
      }
    }
    return benefit;
  }

 private:
  friend class DeltaKm1GainCache;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t penalty_index(const HypernodeID u) const {
    return size_t(u) * ( _k + 1 );
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t benefit_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * ( _k + 1 )  + p + 1;
  }

  // ! Allocates the memory required to store the gain cache
  void allocateGainTable(const HypernodeID num_nodes,
                         const PartitionID k) {
    if (_gain_cache.size() == 0) {
      ASSERT(_k == kInvalidPartition);
      _k = k;
      _gain_cache.resize(
        "Refinement", "gain_cache", num_nodes * size_t(_k + 1), true);
    }
  }

  // ! Initializes the benefit and penalty terms for a node u
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       vec<Gain>& benefit_aggregator) {
    PartitionID from = partitioned_hg.partID(u);
    Gain penalty = 0;
    for (const HyperedgeID& e : partitioned_hg.incidentEdges(u)) {
      HyperedgeWeight ew = partitioned_hg.edgeWeight(e);
      if ( partitioned_hg.pinCountInPart(e, from) > 1 ) {
        penalty += ew;
      }
      for (const PartitionID& i : partitioned_hg.connectivitySet(e)) {
        benefit_aggregator[i] += ew;
      }
    }

    _gain_cache[penalty_index(u)].store(penalty, std::memory_order_relaxed);
    for (PartitionID i = 0; i < _k; ++i) {
      _gain_cache[benefit_index(u, i)].store(benefit_aggregator[i], std::memory_order_relaxed);
      benefit_aggregator[i] = 0;
    }
  }

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


  // ! Indicate whether or not the gain cache is initialized
  bool _is_initialized;

  // ! Number of blocks
  PartitionID _k;

  // ! The gain for moving a node u to from its current block V_i to a target block V_j
  // ! can be expressed as follows for the connectivity metric
  // ! g(u, V_j) := w({ e \in I(u) | pin_count(e, V_i) = 1 }) - w({ e \in I(u) | pin_count(e, V_j) = 0 })
  // !            = w({ e \in I(u) | pin_count(e, V_i) = 1 }) - w(I(u)) + w({ e \in I(u) | pin_count(e, V_j) >= 1 }) (=: b(u, V_j))
  // !            = b(u, V_j) - (w(I(u)) - w({ e \in I(u) | pin_count(e, V_i) = 1 }))
  // !            = b(u, V_j) - w({ e \in I(u) | pin_count(e, V_i) > 1 })
  // !            = b(u, V_j) - p(u)
  // ! We call b(u, V_j) the benefit term and p(u) the penalty term. Our gain cache stores and maintains these
  // ! entries for each node and block. Thus, the gain cache stores k + 1 entries per node.
  ds::Array< CAtomic<HyperedgeWeight> > _gain_cache;
};

class DeltaKm1GainCache {

 public:
  DeltaKm1GainCache(const Km1GainCache& gain_cache) :
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
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, partID(u)) = 1 }))
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID from) const {
    const HyperedgeWeight* penalty_delta =
      _gain_cache_delta.get_if_contained(_gain_cache.penalty_index(u));
    return _gain_cache.penaltyTerm(u, from) + ( penalty_delta ? *penalty_delta : 0 );
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

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID from,
                       const PartitionID to ) const {
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to,
                       const HypernodeID pin_count_in_to_part_after) {
    if (pin_count_in_from_part_after == 1) {
      for (HypernodeID u : partitioned_hg.pins(he)) {
        if (partitioned_hg.partID(u) == from) {
          _gain_cache_delta[_gain_cache.penalty_index(u)] -= we;
        }
      }
    } else if (pin_count_in_from_part_after == 0) {
      for (HypernodeID u : partitioned_hg.pins(he)) {
        _gain_cache_delta[_gain_cache.benefit_index(u, from)] -= we;
      }
    }

    if (pin_count_in_to_part_after == 1) {
      for (HypernodeID u : partitioned_hg.pins(he)) {
        _gain_cache_delta[_gain_cache.benefit_index(u, to)] += we;
      }
    } else if (pin_count_in_to_part_after == 2) {
      for (HypernodeID u : partitioned_hg.pins(he)) {
        if (partitioned_hg.partID(u) == to) {
          _gain_cache_delta[_gain_cache.penalty_index(u)] += we;
        }
      }
    }
  }


 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    utils::MemoryTreeNode* gain_cache_delta_node = parent->addChild("Delta Gain Cache");
    gain_cache_delta_node->updateSize(size_in_bytes());
  }

 private:
  const Km1GainCache& _gain_cache;

  // ! Stores the delta of each locally touched gain cache entry
  // ! relative to the gain cache in '_phg'
  ds::DynamicFlatMap<size_t, HyperedgeWeight> _gain_cache_delta;
};

}  // namespace mt_kahypar
