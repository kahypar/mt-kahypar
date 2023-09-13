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

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/strategies/i_fm_strategy.h"


namespace mt_kahypar {

  /*
   * LocalFMStrategy interface
   * static constexpr bool uses_gain_cache
   * static constexpr bool maintain_gain_cache_between_rounds
   * static constexpr bool is_unconstrained
   *
   * Constructor(context, sharedData, blockPQ, vertexPQs, runStats)
   * insertIntoPQ(phg, gain_cache, node)
   * updateGain(phg, gain_cache, node, move)
   * findNextMove(phg, gain_cache, move)
   * applyMove(phg, gain_cache, move, global)
   * reset()
   * deltaGainUpdates(phg, gain_cache, sync_update)
   *
   */

class LocalUnconstrainedStrategy {
  using VirtualWeightMap = ds::SparseMap<PartitionID, HypernodeWeight>;

 public:
  using BlockPriorityQueue = ds::ExclusiveHandleHeap< ds::MaxHeap<Gain, PartitionID> >;
  using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;    // these need external handles

  static constexpr bool uses_gain_cache = true;
  static constexpr bool maintain_gain_cache_between_rounds = true;
  static constexpr bool is_unconstrained = true;

  LocalUnconstrainedStrategy(const Context& context,
                             FMSharedData& sharedData,
                             BlockPriorityQueue& blockPQ,
                             vec<VertexPriorityQueue>& vertexPQs,
                             FMStats& runStats) :
      context(context),
      runStats(runStats),
      sharedData(sharedData),
      blockPQ(blockPQ),
      vertexPQs(vertexPQs),
      localVirtualWeightDelta(context.partition.k),
      penaltyFactor(context.refinement.fm.imbalance_penalty_max),
      upperBound(context.refinement.fm.unconstrained_upper_bound) { }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void insertIntoPQ(const PartitionedHypergraph& phg,
                    const GainCache& gain_cache,
                    const HypernodeID v) {
    const PartitionID pv = phg.partID(v);
    ASSERT(pv < context.partition.k);
    auto [target, gain] = computeBestTargetBlock(phg, gain_cache, v, pv);
    ASSERT(target < context.partition.k);
    sharedData.targetPart[v] = target;
    vertexPQs[pv].insert(v, gain);  // blockPQ updates are done later, collectively.
    runStats.pushes++;
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void updateGain(const PartitionedHypergraph& phg,
                  const GainCache& gain_cache,
                  const HypernodeID v,
                  const Move& move) {
    const PartitionID pv = phg.partID(v);
    ASSERT(vertexPQs[pv].contains(v));
    const PartitionID designatedTargetV = sharedData.targetPart[v];
    Gain gain = 0;
    PartitionID newTarget = kInvalidPartition;

    if (context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
      // penalty term of designatedTargetV is affected.
      // and may now be greater than that of other blocks --> recompute full
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, gain_cache, v, pv);
    } else {
      // penalty term of designatedTargetV is not affected.
      // only move.from and move.to may be better
      std::tie(newTarget, gain) = bestOfThree(phg, gain_cache,
        v, pv, { designatedTargetV, move.from, move.to });
    }

    sharedData.targetPart[v] = newTarget;
    vertexPQs[pv].adjustKey(v, gain);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool findNextMove(const PartitionedHypergraph& phg,
                    const GainCache& gain_cache,
                    Move& m) {
    updatePQs();

    if (blockPQ.empty()) {
      return false;
    }

    while (true) {
      const PartitionID from = blockPQ.top();
      const HypernodeID u = vertexPQs[from].top();
      const Gain estimated_gain = vertexPQs[from].topKey();
      ASSERT(estimated_gain == blockPQ.topKey());
      auto [to, gain] = computeBestTargetBlock(phg, gain_cache, u, phg.partID(u));

      bool apply_move = (gain >= estimated_gain); // accept any gain that is at least as good
      if (apply_move && to != kInvalidPartition && penaltyFactor > 0) {
        const HypernodeWeight wu = phg.nodeWeight(u);
        const HypernodeWeight to_weight = phg.partWeight(to);
        if (upperBound >= 1 && to_weight + wu > upperBound * context.partition.max_part_weights[to]) {
          apply_move = false;
        } else if (to_weight + wu > context.partition.max_part_weights[to]) {
          const Gain imbalance_penalty = estimatePenalty(to, to_weight, wu);
          if (imbalance_penalty != std::numeric_limits<Gain>::max()) {
            Gain new_gain = gain_cache.gain(u, from, to) - std::ceil(penaltyFactor * imbalance_penalty);
            gain = new_gain;
          } else {
            apply_move = false;
          }
        }
      }

      if (apply_move) {
        m.node = u; m.to = to; m.from = from;
        m.gain = gain;
        runStats.extractions++;
        vertexPQs[from].deleteTop();  // blockPQ updates are done later, collectively.
        return true;
      } else {
        runStats.retries++;
        vertexPQs[from].adjustKey(u, gain);
        sharedData.targetPart[u] = to;
        if (vertexPQs[from].topKey() != blockPQ.keyOf(from)) {
          blockPQ.adjustKey(from, vertexPQs[from].topKey());
        }
      }
    }
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void applyMove(const PartitionedHypergraph& phg, const GainCache&, Move m, bool global) {
    if (sharedData.unconstrained.isRebalancingNode(m.node)) {
      // If a node is moved which is already in use for penalty estimation, we need to make
      // an adjustment so future estimations are not overly optimistic (since in reality, the
      // node is not available anymore). This is achieved by increasing the "virtual" weight of
      // the origin block, thus pessimizing future estimations
      if (global) {
        sharedData.unconstrained.virtual_weight_delta[m.from].fetch_add(
            phg.nodeWeight(m.node), std::memory_order_relaxed);
      } else {
        localVirtualWeightDelta[m.from] += phg.nodeWeight(m.node);
      }
    }
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void revertMove(const PartitionedHypergraph& phg, const GainCache&, Move m, bool global) {
    if (sharedData.unconstrained.isRebalancingNode(m.node)) {
      if (global) {
        sharedData.unconstrained.virtual_weight_delta[m.from].fetch_sub(
            phg.nodeWeight(m.node), std::memory_order_relaxed);
      } else {
        localVirtualWeightDelta[m.from] -= phg.nodeWeight(m.node);
      }
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void flushLocalChanges() {
    for (auto [block, delta]: localVirtualWeightDelta) {
      ASSERT(delta >= 0);
      sharedData.unconstrained.virtual_weight_delta[block].fetch_add(delta, std::memory_order_relaxed);
    }
    localVirtualWeightDelta.clear();
  }

  void reset() {
    // release all nodes that were not moved
    const bool release = sharedData.release_nodes
                         && runStats.moves > 0;

    if (release) {
      // Release all nodes contained in PQ
      for (PartitionID i = 0; i < context.partition.k; ++i) {
        for (PosT j = 0; j < vertexPQs[i].size(); ++j) {
          const HypernodeID v = vertexPQs[i].at(j);
          sharedData.nodeTracker.releaseNode(v);
        }
      }
    }

    for (PartitionID i = 0; i < context.partition.k; ++i) {
      vertexPQs[i].clear();
    }
    blockPQ.clear();
    localVirtualWeightDelta.clear();
  }


  // We're letting the FM details implementation decide what happens here, since some may not want to do gain cache updates,
  // but rather update gains in their PQs or something

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdates(PartitionedHypergraph& phg,
                        GainCache& gain_cache,
                        const SynchronizedEdgeUpdate& sync_update) {
    gain_cache.deltaGainUpdate(phg, sync_update);
  }

  void setPenaltyFactor(double penalty) {
    ASSERT(penalty >= 0 && penalty <= 1);
    penaltyFactor = penalty;
  }

  void setUpperBound(double upper_bound) {
    upperBound = upper_bound;
  }

private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void updatePQs() {
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      if (!vertexPQs[i].empty()) {
        blockPQ.insertOrAdjustKey(i, vertexPQs[i].topKey());
      } else if (blockPQ.contains(i)) {
        blockPQ.remove(i);
      }
    }
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PartitionedHypergraph& phg,
                                                                 const GainCache& gain_cache,
                                                                 const HypernodeID u,
                                                                 const PartitionID from) const {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      if (i != from) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HypernodeWeight max_weight = context.partition.max_part_weights[i];
        HyperedgeWeight benefit = gain_cache.benefitTerm(u, i);
        if (upperBound >= 1 && to_weight + wu > upperBound * max_weight) {
          continue;
        } else if (to_weight + wu > max_weight && benefit <= to_benefit) {
          // don't take imbalanced move without improved gain
          continue;
        } else if (to_weight + wu > max_weight && penaltyFactor > 0) {
          const Gain imbalance_penalty = estimatePenalty(i, to_weight, wu);
          if (imbalance_penalty == std::numeric_limits<Gain>::max()) {
            continue;
          }
          benefit -= std::ceil(penaltyFactor * imbalance_penalty);
        }
        if ( benefit > to_benefit || ( benefit == to_benefit && to_weight < best_to_weight ) ) {
          to_benefit = benefit;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    ASSERT(from == phg.partID(u));
    const Gain gain = to != kInvalidPartition ? to_benefit - gain_cache.penaltyTerm(u, from)
                                              : std::numeric_limits<HyperedgeWeight>::min();
    return std::make_pair(to, gain);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, HyperedgeWeight> bestOfThree(const PartitionedHypergraph& phg,
                                                      const GainCache& gain_cache,
                                                      HypernodeID u,
                                                      PartitionID from,
                                                      std::array<PartitionID, 3> parts) const {

    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i : parts) {
      if (i != from && i != kInvalidPartition) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        HyperedgeWeight benefit = gain_cache.benefitTerm(u, i);
        if (upperBound >= 1 && to_weight + wu > upperBound * context.partition.max_part_weights[i]) {
          continue;
        } else if (to_weight + wu > context.partition.max_part_weights[i] && penaltyFactor > 0) {
          const Gain imbalance_penalty = estimatePenalty(i, to_weight, wu);
          if (imbalance_penalty == std::numeric_limits<Gain>::max()) {
            continue;
          }
          benefit -= std::ceil(penaltyFactor * imbalance_penalty);
        }
        if ( benefit > to_benefit || ( benefit == to_benefit && to_weight < best_to_weight ) ) {
          to_benefit = benefit;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }
    ASSERT(from == phg.partID(u));
    const Gain gain = to != kInvalidPartition ? to_benefit - gain_cache.penaltyTerm(u, from)
                                              : std::numeric_limits<HyperedgeWeight>::min();
    return std::make_pair(to, gain);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  Gain estimatePenalty(PartitionID to, HypernodeWeight to_weight, HypernodeWeight wu) const {
    HypernodeWeight virtual_delta = sharedData.unconstrained.virtual_weight_delta[to].load(std::memory_order_relaxed)
                                    + localVirtualWeightDelta[to];
    HypernodeWeight initial_imbalance = to_weight + virtual_delta - context.partition.max_part_weights[to];
    return sharedData.unconstrained.estimatePenaltyForImbalancedMove(to, initial_imbalance, wu);
  }

  const Context& context;

  FMStats& runStats;

  FMSharedData& sharedData;

  // ! Priority Queue that contains for each block of the partition
  // ! the vertex with the best gain value
  BlockPriorityQueue& blockPQ;

  // ! From PQs -> For each block it contains the vertices (contained
  // ! in that block) touched by the current local search associated
  // ! with their gain values
  vec<VertexPriorityQueue>& vertexPQs;

  // ! Virtual block weights are saved as delta to the actual block weight. They
  // ! are necessary to ensure a reasonable penalty estimation in some edge cases.
  VirtualWeightMap localVirtualWeightDelta;

  double penaltyFactor;
  double upperBound;
};


template<typename TypeTraits, typename GainTypes>
class UnconstrainedStrategy: public IFMStrategy {
  using Base = IFMStrategy;
  static constexpr bool debug = false;

 public:
  using LocalFM = LocalizedKWayFM<TypeTraits, GainTypes>;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

  UnconstrainedStrategy(const Context& context, FMSharedData& sharedData):
      Base(context, sharedData),
      current_penalty(context.refinement.fm.imbalance_penalty_min),
      current_upper_bound(context.refinement.fm.unconstrained_upper_bound),
      absolute_improvement_first_round(kInvalidGain),
      unconstrained_is_enabled(true),
      stats(utils::Utilities::instance().getStats(context.utility_id)) {
        ASSERT(!context.refinement.fm.activate_unconstrained_dynamically
                || context.refinement.fm.multitry_rounds > 2);
  }

  bool dispatchedFindMoves(LocalFM& local_fm, PartitionedHypergraph& phg, size_t task_id, size_t num_seeds, size_t round) {
    if (isUnconstrainedRound(round)) {
      LocalUnconstrainedStrategy local_strategy = local_fm.template initializeDispatchedStrategy<LocalUnconstrainedStrategy>();
      local_strategy.setPenaltyFactor(current_penalty);
      local_strategy.setUpperBound(current_upper_bound);
      return local_fm.findMoves(local_strategy, phg, task_id, num_seeds);
    } else {
      LocalGainCacheStrategy local_strategy = local_fm.template initializeDispatchedStrategy<LocalGainCacheStrategy>();
      return local_fm.findMoves(local_strategy, phg, task_id, num_seeds);
    }
  }

 private:
  virtual void findMovesImpl(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& phg,
                             size_t num_tasks, size_t num_seeds, size_t round) final {
    Base::findMovesWithConcreteStrategy<UnconstrainedStrategy>(
              local_fm, phg, num_tasks, num_seeds, round);
  }

  virtual bool isUnconstrainedRoundImpl(size_t round) const final {
    if (round > 0 && !unconstrained_is_enabled) {
      return false;
    }
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      return round == 1 || (round > 1 && round - 2 < context.refinement.fm.unconstrained_rounds);
    } else {
      return round < context.refinement.fm.unconstrained_rounds;
    }
  }

  virtual bool includesUnconstrainedImpl() const final {
    return true;
  }

  virtual void reportImprovementImpl(size_t round, Gain absolute_improvement, double relative_improvement) final {
    if (round == 0) {
      absolute_improvement_first_round = absolute_improvement;
    } else if (round == 1
               && context.refinement.fm.activate_unconstrained_dynamically
               && absolute_improvement < absolute_improvement_first_round) {
        // this is the decision point whether unconstrained or constrained FM is used
        unconstrained_is_enabled = false;
        DBG << "Disabling unconstrained FM after test round: " << V(absolute_improvement) << V(absolute_improvement_first_round);
    } else if (relative_improvement < context.refinement.fm.unconstrained_min_improvement) {
      unconstrained_is_enabled = false;
      DBG << "Disabling unconstrained FM due to too little improvement:" << V(relative_improvement);
    }
    if (round == 1) {
      stats.update_stat("top-level-ufm-active", unconstrained_is_enabled);
      if (unconstrained_is_enabled) {
        stats.update_stat("ufm-active-levels", 1);
      } else {
        stats.update_stat("ufm-inactive-levels", 1);
      }
    }
  }

  void initRound(size_t /*num_tasks*/, size_t /*num_seeds*/, size_t round) final {
    if (round == 0) {
      unconstrained_is_enabled = true;
    }
    if (context.refinement.fm.activate_unconstrained_dynamically) {
      if (round == 1) {
        current_penalty = context.refinement.fm.penalty_for_activation_test;
        current_upper_bound = context.refinement.fm.unconstrained_upper_bound;
      } else if (round > 1 && isUnconstrainedRound(round)) {
        size_t n_rounds = std::min(context.refinement.fm.unconstrained_rounds, context.refinement.fm.multitry_rounds - 2);
        calculateInterpolation(round - 2, n_rounds);
      }
    } else if (isUnconstrainedRound(round)) {
      calculateInterpolation(round, context.refinement.fm.unconstrained_rounds);
    }
    DBG << V(round) << V(isUnconstrainedRound(round)) << V(current_penalty) << V(current_upper_bound);
  }

  void calculateInterpolation(size_t round, size_t n_rounds) {
    ASSERT(unconstrained_is_enabled && round < context.refinement.fm.multitry_rounds);
    auto interpolate = [&](double start, double end) {
      if (round == 0) {
        return start;
      }
      double summed = (n_rounds - round - 1) * start + round * end;
      return summed / static_cast<double>(n_rounds - 1);
    };

    if (round < n_rounds) {
      // interpolate values for current penalty and upper bound
      current_penalty = interpolate(context.refinement.fm.imbalance_penalty_min,
                                    context.refinement.fm.imbalance_penalty_max);
      if (context.refinement.fm.unconstrained_upper_bound >= 1) {
        if (context.refinement.fm.unconstrained_upper_bound_min >= 1) {
          current_upper_bound = interpolate(context.refinement.fm.unconstrained_upper_bound,
                                            context.refinement.fm.unconstrained_upper_bound_min);
        } else {
          current_upper_bound = context.refinement.fm.unconstrained_upper_bound;
        }
      }
    }
  }

  double current_penalty;
  double current_upper_bound;
  Gain absolute_improvement_first_round;
  bool unconstrained_is_enabled;
  utils::Stats& stats;
};

}
