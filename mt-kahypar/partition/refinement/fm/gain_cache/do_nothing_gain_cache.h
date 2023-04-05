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

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

// Forward
class DeltaDoNothingGainCache;

class DoNothingGainCache final : public kahypar::meta::PolicyBase {

 public:
  static constexpr FMGainCacheType TYPE = FMGainCacheType::none;
  using DeltaGainCache = DeltaDoNothingGainCache;

  DoNothingGainCache() { }

  DoNothingGainCache(const DoNothingGainCache&) = delete;
  DoNothingGainCache & operator= (const DoNothingGainCache &) = delete;

  DoNothingGainCache(DoNothingGainCache&& other) = default;
  DoNothingGainCache & operator= (DoNothingGainCache&& other) = default;

  // ####################### Initialization #######################

  bool isInitialized() const {
    return false;
  }

  void reset(const bool run_parallel = true) {
    unused(run_parallel);
  }

  size_t size() const {
    return 0;
  }

  template<typename PartitionedHypergraph>
  void initializeGainCache(const PartitionedHypergraph&) {
    // do nothing
  }

  // ####################### Gain Computation #######################

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID,
                              const PartitionID) const {
    return 0;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputePenaltyTermEntry(const PartitionedHypergraph&,
                                 const HypernodeID) {
    // Do nothing
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID, const PartitionID) const {
    return 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID,
                       const PartitionID,
                       const PartitionID ) const {
    return 0;
  }

  // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph&, const HyperedgeID, const HyperedgeWeight,
                       const PartitionID, const HypernodeID, const PartitionID, const HypernodeID) {
    // Do nothing
  }

  // ####################### Uncontraction #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterRestore(const PartitionedHypergraph&, const HypernodeID, const HypernodeID,
                                    const HyperedgeID, const HypernodeID) {
    // Do nothing
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void uncontractUpdateAfterReplacement(const PartitionedHypergraph&, const HypernodeID,
                                        const HypernodeID, const HyperedgeID) {
    // Do nothing
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void restoreSinglePinHyperedge(const HypernodeID, const PartitionID, const HyperedgeWeight) {
    // Do nothing
  }

  // ####################### Only for Testing #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedHypergraph&, const HypernodeID) const {
    return 0;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedHypergraph&, const HypernodeID, const PartitionID) const {
    return 0;
  }

};

class DeltaDoNothingGainCache {

 public:
  DeltaDoNothingGainCache(const DoNothingGainCache&) { }

  // ####################### Initialize & Reset #######################

  void initialize(const size_t) {
    // Do nothing
  }

  void clear() {
    // Do nothing
  }

  void dropMemory() {
    // Do nothing
  }

  size_t size_in_bytes() const {
    return 0;
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, partID(u)) = 1 }))
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID, const PartitionID) const {
    return 0;
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) >= 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID, const PartitionID) const {
    return 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID, const PartitionID, const PartitionID ) const {
    return 0;
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph&, const HyperedgeID, const HyperedgeWeight,
                       const PartitionID, const HypernodeID, const PartitionID, const HypernodeID) {
    // Do nothing
  }

 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode*) const {
    // Do nothing
  }
};

}  // namespace mt_kahypar
