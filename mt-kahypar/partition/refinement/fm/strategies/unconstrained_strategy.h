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


// TODO: HIGH_DEGREE_THRESHOLD in PartitionedHypergraph/PartitionedGraph might be problematic
// for unconstrained refinement


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
   * skipMove(phg, gain_cache, move)
   * clearPQs()
   * deltaGainUpdates(phg, gain_cache, he, edge_weight, from,
   *                  pin_count_in_from_part_after, to,  pin_count_in_to_part_after)
   *
   */

class LocalUnconstrainedStrategy {

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
          const HypernodeWeight imbalance = std::min(wu, to_weight + wu - context.partition.max_part_weights[to]);
          // The following will update the imbalance globally, which also affects the imbalance penalty for other threads.
          // If the move is not applied, we need to undo this in skipMove
          const Gain imbalance_penalty = sharedData.unconstrained.applyEstimatedPenaltyForImbalancedMove(to, imbalance);
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

        if (context.refinement.fm.penalty_for_moved_rebalancing_nodes && sharedData.unconstrained.isRebalancingNode(m.node)) {
          // edge case: moving a rebalancing node can throw the estimation off if we don't apply a correction
          const HypernodeWeight from_weight = phg.partWeight(m.from);
          const HypernodeWeight wu = phg.nodeWeight(u);
          if (from_weight - wu < context.partition.max_part_weights[from]) {
            const HypernodeWeight reduction = std::min(wu, context.partition.max_part_weights[from] - from_weight + wu);
            sharedData.unconstrained.applyEstimatedPenaltyForImbalancedMove(from, reduction);
          }
        }
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
  void skipMove(const PartitionedHypergraph& phg, const GainCache&, Move m) {
    if (penaltyFactor > 0) {
      const HypernodeWeight to_weight = phg.partWeight(m.to);
      if (to_weight > context.partition.max_part_weights[m.to]) {
        // we need to undo the imbalance which was added to the shared data
        const HypernodeWeight hn_weight = phg.nodeWeight(m.node);
        const HypernodeWeight imbalance = std::min(hn_weight, to_weight - context.partition.max_part_weights[m.to]);
        sharedData.unconstrained.revertImbalancedMove(m.to, imbalance);

        // if (sharedData.unconstrained.isRebalancingNode(m.node)) {
            // edge case: undo moving a rebalancing node
            // Probably nothing to do here, since this is extremely unlikely and pessimizations are unproblematic
        // }
      }
    }
  }

  void clearPQs(const size_t /* bestImprovementIndex */ ) {
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
  }


  // We're letting the FM details implementation decide what happens here, since some may not want to do gain cache updates,
  // but rather update gains in their PQs or something

  template<typename PartitionedHypergraph, typename GainCache>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdates(PartitionedHypergraph& phg,
                        GainCache& gain_cache,
                        const HyperedgeID he,
                        const HyperedgeWeight edge_weight,
                        const PartitionID from,
                        const HypernodeID pin_count_in_from_part_after,
                        const PartitionID to,
                        const HypernodeID pin_count_in_to_part_after) {
    gain_cache.deltaGainUpdate(phg, he, edge_weight, from,
      pin_count_in_from_part_after, to, pin_count_in_to_part_after);
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
    // TODO(maas): bonus for balancing moves?!
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
          const HypernodeWeight imbalance = std::min(wu, to_weight + wu - max_weight);
          const Gain imbalance_penalty = sharedData.unconstrained.estimatedPenaltyForImbalancedMove(i, imbalance);
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
          const HypernodeWeight imbalance = std::min(wu, to_weight + wu - context.partition.max_part_weights[i]);
          const Gain imbalance_penalty = sharedData.unconstrained.estimatedPenaltyForImbalancedMove(i, imbalance);
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

  double penaltyFactor;
  double upperBound;
};


template<typename TypeTraits, typename GainTypes>
class UnconstrainedStrategy: public IFMStrategy {
  using Base = IFMStrategy;

 public:
  using LocalFM = LocalizedKWayFM<TypeTraits, GainTypes>;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

  UnconstrainedStrategy(const Context& context, FMSharedData& sharedData):
      Base(context, sharedData) { }

  bool dispatchedFindMoves(LocalFM& local_fm, PartitionedHypergraph& phg, size_t task_id, size_t num_seeds, size_t) {
    LocalUnconstrainedStrategy local_strategy = local_fm.template initializeDispatchedStrategy<LocalUnconstrainedStrategy>();
    return local_fm.findMoves(local_strategy, phg, task_id, num_seeds);
  }

 private:
  virtual void findMovesImpl(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& phg,
                             size_t num_tasks, size_t num_seeds, size_t round,
                             ds::StreamingVector<HypernodeID>& locally_locked_vertices) final {
    Base::findMovesWithConcreteStrategy<UnconstrainedStrategy>(
              local_fm, phg, num_tasks, num_seeds, round, locally_locked_vertices);
  }

  virtual bool isUnconstrainedRoundImpl(size_t) const final {
    return true;
  }

  virtual bool includesUnconstrainedImpl() const final {
    return true;
  }
};

}
