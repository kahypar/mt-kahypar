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
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {
template <typename TypeTraits, typename GainTypes>
class JetRebalancer final : public IRebalancer {
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
    _max_part_weights(nullptr),
    _gain_cache(gain_cache),
    _current_k(_context.partition.k),
    _gain(context),
    _num_imbalanced_blocks(0),
    _num_valid_targets(0),
    _part_weights(_context.partition.k),
    _buckets(),
    _bucket_weights(_context.partition.k * NUM_BUCKETS, 0),
    _local_bucket_weights([&] {
      return constructBucketWeightVector();
    }),
    _node_was_moved() {
      for (size_t i = 0; i < NUM_BUCKETS; ++i) {
        _buckets.emplace_back(BUCKET_FACTOR);
      }
    }

  explicit JetRebalancer(const Context& context,
                         gain_cache_t gain_cache) :
    JetRebalancer(context, GainCachePtr::cast<GainCache>(gain_cache)) {}

  JetRebalancer(const JetRebalancer&) = delete;
  JetRebalancer(JetRebalancer&&) = delete;

  JetRebalancer & operator= (const JetRebalancer &) = delete;
  JetRebalancer & operator= (JetRebalancer &&) = delete;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  double) {
    ASSERT(refinement_nodes.empty());
    unused(refinement_nodes);
    return refineInternal(hypergraph, nullptr, best_metrics);
  }

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) final {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    _node_was_moved.resize(phg.initialNumNodes(), uint8_t(false));  // TODO
  }

  void setMaxPartWeightsForRoundImpl(const std::vector<HypernodeWeight>& max_part_weights) final {
    _max_part_weights = &max_part_weights[0];
  }

  bool refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                const vec<HypernodeID>& refinement_nodes,
                                vec<vec<Move>>& moves_by_part,
                                Metrics& best_metrics,
                                const double) {
    ASSERT(refinement_nodes.empty());
    unused(refinement_nodes);
    return refineInternal(hypergraph, &moves_by_part, best_metrics);
  }

private:
  bool refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                      vec<vec<Move>>* moves_by_part,
                      Metrics& best_metric);

  template<bool ensure_balanced_moves>
  void weakRebalancingRound(PartitionedHypergraph& phg);

  void strongRebalancingRound(PartitionedHypergraph& phg);

  template<typename F>
  void insertNodesIntoBuckets(PartitionedHypergraph& phg, F compute_gain_fn);

  template<typename F>
  void processBuckets(PartitionedHypergraph& phg, F move_node_fn, bool retry_invalid_moves, bool update_local_part_weights);

  std::pair<Gain, PartitionID> computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
                                                        const HypernodeID hn,
                                                        bool non_adjacent_blocks,
                                                        bool use_precise_part_weights = false,
                                                        bool use_deadzone = true);

  // used for Jetrs (strong rebalancing), rounded down
  Gain computeAverageGain(const PartitionedHypergraph& hypergraph, const HypernodeID hn);

  void updateImbalance(const PartitionedHypergraph& hypergraph, bool read_weights_from_graph = true);

  void initializeDataStructures(const PartitionedHypergraph& hypergraph);

  bool isBalanced(const PartitionedHypergraph& phg) {
    ASSERT(_max_part_weights != nullptr);
    for (PartitionID i = 0; i < _context.partition.k; ++i) {
      if (phg.partWeight(i) > _max_part_weights[i]) {
        return false;
      }
    }
    return true;
  }

  bool mayMoveNode(PartitionID block, HypernodeWeight hn_weight) const {
    double allowed_weight = _part_weights[block].load(std::memory_order_relaxed)
                            - _context.partition.perfect_balance_part_weights[block];
    allowed_weight *= _context.refinement.jet_rebalancing.heavy_vertex_exclusion_factor;
    return hn_weight <= allowed_weight;
  }

  HypernodeWeight imbalance(PartitionID block) const {
    return _part_weights[block].load(std::memory_order_relaxed) - _max_part_weights[block];
  }

  HypernodeWeight deadzoneForPart(PartitionID block) const {
    const HypernodeWeight balanced = _context.partition.perfect_balance_part_weights[block];
    const HypernodeWeight max = _max_part_weights[block];
    return max - _context.refinement.jet_rebalancing.relative_deadzone_size * (max - balanced);
  }

  bool isValidTarget(const PartitionedHypergraph& hypergraph,
                     PartitionID block,
                     HypernodeWeight hn_weight,
                     bool use_precise_part_weights,
                     bool use_deadzone = true) const {
    const HypernodeWeight block_weight = use_precise_part_weights ?
        hypergraph.partWeight(block) : _part_weights[block].load(std::memory_order_relaxed);
    return (!use_deadzone || block_weight < deadzoneForPart(block)) &&
           block_weight + hn_weight <= _max_part_weights[block];
  }

  bool changeNodePart(PartitionedHypergraph& phg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to,
                      bool ensure_balanced) {
    // it happens spuriously that from == to, not entirely sure why (possibly due to moving heavy nodes)
    if (from == to || to == kInvalidPartition) {
      return false;
    }

    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SyncronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };

    HypernodeWeight max_weight = ensure_balanced ? _max_part_weights[to] : std::numeric_limits<HypernodeWeight>::max();
    bool success = false;
    if ( _gain_cache.isInitialized() ) {
      success = phg.changeNodePart(_gain_cache, hn, from, to, max_weight, []{}, objective_delta);
    } else {
      success = phg.changeNodePart(hn, from, to, max_weight, []{}, objective_delta);
    }
    ASSERT(success || ensure_balanced);
    if (success) {
      _node_was_moved[hn] = uint8_t(true);
    }
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
  const HypernodeWeight* _max_part_weights;
  GainCache& _gain_cache;
  PartitionID _current_k;
  GainCalculator _gain;
  PartitionID _num_imbalanced_blocks;
  PartitionID _num_valid_targets;
  parallel::scalable_vector<AtomicWeight> _part_weights;
  parallel::scalable_vector<BucketMap> _buckets;
  parallel::scalable_vector<HypernodeWeight> _bucket_weights;
  tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> _local_bucket_weights;
  parallel::scalable_vector<uint8_t> _node_was_moved;
};

}  // namespace kahypar
