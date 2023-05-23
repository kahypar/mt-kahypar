/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_gain_computation.h"

namespace mt_kahypar {
template <typename TypeTraits, typename GainTypes>
class JetRebalancer final : public IRefiner {
 private:
  using BucketMap = ds::ConcurrentBucketMap<HypernodeID>;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using GainCalculator = typename GainTypes::GainComputation;
  using RatingMap = typename GainCalculator::RatingMap;
  using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

  static constexpr size_t NUM_BUCKETS = 12;
  static constexpr size_t BUCKET_FACTOR = 32;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

public:
  // TODO: add repairEmptyBlocks functionality?

  explicit JetRebalancer(const Context& context,
                         GainCache& gain_cache) :
    _context(context),
    _gain_cache(gain_cache),
    _current_k(_context.partition.k),
    _gain(context),
    _part_weights(_context.partition.k),
    _buckets(),
    _bucket_weights(_context.partition.k * NUM_BUCKETS, 0),
    _local_bucket_weights([&] {
      return constructBucketWeightVector();
    }) {
      for (size_t i = 0; i < NUM_BUCKETS; ++i) {
        _buckets.emplace_back(BUCKET_FACTOR);
      }
    }

  JetRebalancer(const JetRebalancer&) = delete;
  JetRebalancer(JetRebalancer&&) = delete;

  JetRebalancer & operator= (const JetRebalancer &) = delete;
  JetRebalancer & operator= (JetRebalancer &&) = delete;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>&,
                  Metrics& best_metrics,
                  double);

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final { }

private:
  template<bool ensure_balanced_moves>
  void rebalancingRound(PartitionedHypergraph& phg);

  std::pair<Gain, PartitionID> computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
                                                        const HypernodeID hn,
                                                        bool non_adjacent_blocks,
                                                        bool use_precise_part_weights = false);

  PartitionID updateImbalance(const PartitionedHypergraph& hypergraph);

  void initializeDataStructures(const PartitionedHypergraph& hypergraph);

  HypernodeWeight imbalance(PartitionID block) const {
    return _part_weights[block].load(std::memory_order_relaxed)
           - _context.partition.max_part_weights[block];
  }

  bool isAcceptableTarget(const PartitionedHypergraph& hypergraph,
                          PartitionID block,
                          HypernodeWeight hn_weight,
                          bool use_precise_part_weights) const {
    const HypernodeWeight block_weight = use_precise_part_weights ?
        hypergraph.partWeight(block) : _part_weights[block].load(std::memory_order_relaxed);
    // TODO: use variable deadzone?
    return block_weight < _context.partition.perfect_balance_part_weights[block] &&
           block_weight + hn_weight <= _context.partition.max_part_weights[block];
  }

  bool changeNodePart(PartitionedHypergraph& phg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to,
                      bool ensure_balanced) {
    if (to == kInvalidPartition) {
      return false;
    }

    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      _gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // TODO: parameter for unconstrained move
    HypernodeWeight max_weight = ensure_balanced ? _context.partition.max_part_weights[to]
                                  : std::numeric_limits<HypernodeWeight>::max();
    bool success = false;
    if ( _context.forceGainCacheUpdates() && _gain_cache.isInitialized() ) {
      success = phg.changeNodePart(_gain_cache, hn, from, to, max_weight, []{}, objective_delta);
    } else {
      success = phg.changeNodePart(hn, from, to, max_weight, []{}, objective_delta);
    }
    ASSERT(success || ensure_balanced);
    return success;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t log2(Gain gain) const {
    size_t value = 0;
    while (gain >>= 1) {
      ++value;
    }
    return value;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t getBucketID(Gain gain) const {
    if (gain > 0) {
      return std::min(2 + log2(gain), NUM_BUCKETS - 1);
    } else if (gain == 0) {
      return 1;
    }
    return 0;
  }

  void resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if ( _current_k != _context.partition.k ) {
      _current_k = _context.partition.k;
      _gain.changeNumberOfBlocks(_current_k);
      _part_weights = parallel::scalable_vector<AtomicWeight>(_current_k);
    }
  }

  parallel::scalable_vector<HypernodeWeight> constructBucketWeightVector() const {
    return parallel::scalable_vector<HypernodeWeight>(_current_k * NUM_BUCKETS, 0);
  }

  const Context& _context;
  GainCache& _gain_cache;
  PartitionID _current_k;
  GainCalculator _gain;
  parallel::scalable_vector<AtomicWeight> _part_weights;
  parallel::scalable_vector<BucketMap> _buckets;
  parallel::scalable_vector<HypernodeWeight> _bucket_weights;
  tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> _local_bucket_weights;
};

}  // namespace kahypar
