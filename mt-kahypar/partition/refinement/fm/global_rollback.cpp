/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/fm/global_rollback.h"

#include "tbb/parallel_scan.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"

namespace mt_kahypar {
  struct Prefix {
    Gain gain = 0;                           /** gain when using valid moves up to best_index */
    Gain penalty = 0;                        /** penalty for imbalanced moved in unconstrained case */
    MoveID best_index = 0;                   /** local ID of first move to revert */
    HypernodeWeight heaviest_weight =
            std::numeric_limits<HypernodeWeight>::max();   /** weight of the heaviest part */

    bool operator<(const Prefix& o) const {
      return gain - penalty > o.gain - o.penalty ||
              (gain - penalty == o.gain - o.penalty && std::tie(heaviest_weight, best_index) < std::tie(o.heaviest_weight, o.best_index));
    }
  };

  template<typename PartitionedHypergraph, typename Derived>
  struct BalanceAndBestIndexScanBase {
    const PartitionedHypergraph& phg;
    const vec<Move>& moves;
    std::shared_ptr< tbb::enumerable_thread_specific<Prefix> > local_best;

    Gain gain_sum = 0;

    vec<HypernodeWeight> part_weights;
    const std::vector<HypernodeWeight>& max_part_weights;

    BalanceAndBestIndexScanBase(Derived& b, tbb::split) :
            phg(b.phg),
            moves(b.moves),
            local_best(b.local_best),
            gain_sum(0),
            part_weights(b.part_weights.size(), 0),
            max_part_weights(b.max_part_weights) { }


    BalanceAndBestIndexScanBase(const PartitionedHypergraph& phg,
                                const vec<Move>& moves,
                                const vec<HypernodeWeight>& part_weights,
                                const std::vector<HypernodeWeight>& max_part_weights) :
            phg(phg),
            moves(moves),
            local_best(std::make_shared< tbb::enumerable_thread_specific<Prefix> >()),
            gain_sum(0),
            part_weights(part_weights),
            max_part_weights(max_part_weights)
    {
    }

    // subranges a | b | c | d . assuming this ran pre_scan on c,
    // then lhs ran pre_scan on b and final_scan of this will be on d
    void reverse_join(Derived& lhs) {
      for (size_t i = 0; i < part_weights.size(); ++i) {
        part_weights[i] += lhs.part_weights[i];
      }
      gain_sum += lhs.gain_sum;
    }

    void assign(Derived& b) {
      gain_sum = b.gain_sum;
    }
  };

  template<typename PartitionedHypergraph>
  struct BalanceAndBestIndexScan: public BalanceAndBestIndexScanBase<PartitionedHypergraph, BalanceAndBestIndexScan<PartitionedHypergraph>>  {
    using Base = BalanceAndBestIndexScanBase<PartitionedHypergraph, BalanceAndBestIndexScan<PartitionedHypergraph>>;
    using Base::phg, Base::moves, Base::local_best, Base::gain_sum, Base::part_weights, Base::max_part_weights;

    BalanceAndBestIndexScan(BalanceAndBestIndexScan& b, tbb::split) :
            Base(b, tbb::split()) { }

    BalanceAndBestIndexScan(const PartitionedHypergraph& phg,
                            const vec<Move>& moves,
                            const vec<HypernodeWeight>& part_weights,
                            const std::vector<HypernodeWeight>& max_part_weights) :
            Base(phg, moves, part_weights, max_part_weights) { }

    void operator()(const tbb::blocked_range<MoveID>& r, tbb::pre_scan_tag ) {
      for (MoveID i = r.begin(); i < r.end(); ++i) {
        const Move& m = moves[i];
        if (m.isValid()) {  // skip locally reverted moves
          gain_sum += m.gain;
          part_weights[m.from] -= phg.nodeWeight(m.node);
          part_weights[m.to] += phg.nodeWeight(m.node);
        }
      }
    }

    void operator()(const tbb::blocked_range<MoveID>& r, tbb::final_scan_tag ) {
      size_t overloaded = 0;
      for (size_t i = 0; i < part_weights.size(); ++i) {
        if (part_weights[i] > max_part_weights[i]) {
          overloaded++;
        }
      }

      Prefix current;
      for (MoveID i = r.begin(); i < r.end(); ++i) {
        const Move& m = moves[i];

        if (m.isValid()) {  // skip locally reverted moves
          gain_sum += m.gain;

          const bool from_overloaded = part_weights[m.from] > max_part_weights[m.from];
          part_weights[m.from] -= phg.nodeWeight(m.node);
          if (from_overloaded && part_weights[m.from] <= max_part_weights[m.from]) {
            overloaded--;
          }
          const bool to_overloaded = part_weights[m.to] > max_part_weights[m.to];
          part_weights[m.to] += phg.nodeWeight(m.node);
          if (!to_overloaded && part_weights[m.to] > max_part_weights[m.to]) {
            overloaded++;
          }

          if (overloaded == 0 && gain_sum >= current.gain) {
            Prefix new_prefix = { gain_sum, 0, i + 1, *std::max_element(part_weights.begin(), part_weights.end()) };
            current = std::min(current, new_prefix);
          }
        }
      }

      if (current.best_index != 0) {
        Prefix& lb = local_best->local();
        lb = std::min(lb, current);
      }
    }

    Prefix finalize(const vec<HypernodeWeight>& initial_part_weights) {
      Prefix res { 0, 0, 0, *std::max_element(initial_part_weights.begin(), initial_part_weights.end()) };
      for (const Prefix& x : *local_best) {
        res = std::min(res, x);
      }
      return res;
    }
  };


  template<typename PartitionedHypergraph>
  struct UnconstrainedBalanceAndBestIndexScan: public BalanceAndBestIndexScanBase<PartitionedHypergraph,
                                                          UnconstrainedBalanceAndBestIndexScan<PartitionedHypergraph>>  {
    using Base = BalanceAndBestIndexScanBase<PartitionedHypergraph, UnconstrainedBalanceAndBestIndexScan<PartitionedHypergraph>>;
    using Base::phg, Base::moves, Base::local_best, Base::gain_sum, Base::part_weights, Base::max_part_weights;

    const UnconstrainedFMData& unconstrainedData;

    UnconstrainedBalanceAndBestIndexScan(UnconstrainedBalanceAndBestIndexScan& b, tbb::split) :
            Base(b, tbb::split()), unconstrainedData(b.unconstrainedData) { }

    UnconstrainedBalanceAndBestIndexScan(const PartitionedHypergraph& phg,
                                         const vec<Move>& moves,
                                         const vec<HypernodeWeight>& part_weights,
                                         const std::vector<HypernodeWeight>& max_part_weights,
                                         const UnconstrainedFMData& unconstrainedData) :
            Base(phg, moves, part_weights, max_part_weights), unconstrainedData(unconstrainedData) { }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isRebalancingNode(const Move& m) const {
      return unconstrainedData.isRebalancingNode(m.node); // TODO: use threshold here!!!
    }

    void operator()(const tbb::blocked_range<MoveID>& r, tbb::pre_scan_tag ) {
      for (MoveID i = r.begin(); i < r.end(); ++i) {
        const Move& m = moves[i];
        if (m.isValid()) {  // skip locally reverted moves
          gain_sum += m.gain;
          // we need to be careful with nodes that are already "reserved" for rebalancing, otherwise
          // the penalty estimation for imbalanced blocks might be wrong
          // => pessimistic estimation: don't reduce part weight if rebalancing node is moved
          if (!isRebalancingNode(m)) {
            part_weights[m.from] -= phg.nodeWeight(m.node);
          }
          part_weights[m.to] += phg.nodeWeight(m.node);
        }
      }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HypernodeWeight weightLimitForPart(PartitionID block) const {
      return max_part_weights[block] + unconstrainedData.maximumImbalance(block);
    }

    Gain determinePenalty(const vec<HypernodeWeight>& input_part_weights) {
      Gain penalty = 0;
      for (size_t i = 0; i < input_part_weights.size(); ++i) {
        HypernodeWeight imbalance = input_part_weights[i] - max_part_weights[i];
        if (imbalance > 0) {
          penalty += unconstrainedData.estimatedPenaltyForImbalance(i, imbalance);
        }
      }
      return penalty;
    }

    void operator()(const tbb::blocked_range<MoveID>& r, tbb::final_scan_tag ) {
      Gain penalty = 0;
      size_t overloaded = 0;
      for (size_t i = 0; i < part_weights.size(); ++i) {
        if (part_weights[i] > weightLimitForPart(i)) {
          overloaded++;
        }
      }
      if (overloaded == 0) {
        penalty = determinePenalty(part_weights);
      }

      Prefix current;
      for (MoveID i = r.begin(); i < r.end(); ++i) {
        const Move& m = moves[i];

        if (m.isValid()) {  // skip locally reverted moves
          const HypernodeWeight node_weight = phg.nodeWeight(m.node);
          gain_sum += m.gain;

          const HypernodeWeight to_imbalance = part_weights[m.to] - max_part_weights[m.to];
          part_weights[m.to] += node_weight;
          if (to_imbalance + node_weight > 0) {
            const bool to_overloaded = to_imbalance > unconstrainedData.maximumImbalance(m.to);
            if (!to_overloaded && part_weights[m.to] > weightLimitForPart(m.to)) {
              overloaded++;
            } else if (overloaded == 0) {
              // this is the expected case: we update the penalty with the delta for the current block
              penalty += unconstrainedData.estimatedPenaltyForDelta(m.to,
                            std::max(to_imbalance, 0), to_imbalance + node_weight);
            }
          }

          // we need to be careful with nodes that are already "reserved" for rebalancing, otherwise
          // the penalty estimation for imbalanced blocks might be wrong
          // => pessimistic estimation: don't reduce part weight if rebalancing node is moved
          if (!isRebalancingNode(m)) {
            const HypernodeWeight from_imbalance = part_weights[m.from] - max_part_weights[m.from];
            part_weights[m.from] -= node_weight;
            if (from_imbalance > 0) {
              const bool from_overloaded = from_imbalance > unconstrainedData.maximumImbalance(m.from);
              if (from_overloaded && part_weights[m.from] <= weightLimitForPart(m.from)) {
                overloaded--;
              }
              if (!from_overloaded && overloaded == 0) {
                // this is the expected case: we update the penalty with the delta for the current block
                penalty -= unconstrainedData.estimatedPenaltyForDelta(m.from,
                              std::max(from_imbalance - node_weight, 0), from_imbalance);
              } else if (overloaded == 0) {
                // we need to determine the penalty anew
                penalty = determinePenalty(part_weights);
              }
            }
          }

          if (overloaded == 0 && gain_sum - penalty >= current.gain - current.penalty) {
            Prefix new_prefix = { gain_sum, penalty, i + 1, *std::max_element(part_weights.begin(), part_weights.end()) };
            current = std::min(current, new_prefix);
          }
        }
      }

      if (current.best_index != 0) {
        Prefix& lb = local_best->local();
        lb = std::min(lb, current);
      }
    }

    Prefix finalize(const vec<HypernodeWeight>& initial_part_weights) {
      Prefix res { 0, 0, 0, *std::max_element(initial_part_weights.begin(), initial_part_weights.end()) };
      res.penalty = determinePenalty(initial_part_weights);
      for (const Prefix& x : *local_best) {
        res = std::min(res, x);
      }
      return res;
    }
  };

  template<typename TypeTraits, typename GainTypes>
  HyperedgeWeight GlobalRollback<TypeTraits, GainTypes>::revertToBestPrefixParallel(
          PartitionedHypergraph& phg, FMSharedData& sharedData, const vec<HypernodeWeight>& partWeights,
          const std::vector<HypernodeWeight>& maxPartWeights, bool unconstrained) {
    const MoveID numMoves = sharedData.moveTracker.numPerformedMoves();
    if (numMoves == 0) return 0;

    const vec<Move>& move_order = sharedData.moveTracker.moveOrder;

    recalculateGains(phg, sharedData);
    HEAVY_REFINEMENT_ASSERT(verifyGains(phg, sharedData));

    // TODO set grain size in blocked_range? to avoid too many copies of part weights array. experiment with different values
    Prefix b;
    if (unconstrained) {
      UnconstrainedBalanceAndBestIndexScan<PartitionedHypergraph> s(phg, move_order,
          partWeights, maxPartWeights, sharedData.unconstrained);
      tbb::parallel_scan(tbb::blocked_range<MoveID>(0, numMoves), s);
      b = s.finalize(partWeights);
    } else {
      BalanceAndBestIndexScan<PartitionedHypergraph> s(phg, move_order, partWeights, maxPartWeights);
      tbb::parallel_scan(tbb::blocked_range<MoveID>(0, numMoves), s);
      b = s.finalize(partWeights);
    }

    tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
      const Move& m = move_order[moveID];
      if (m.isValid()) {
        moveVertex(phg, m.node, m.to, m.from);
      }
    });

    // recompute penalty term values since they are potentially invalid
    tbb::parallel_for(MoveID(0), numMoves, [&](const MoveID i) {
      gain_cache.recomputePenaltyTermEntry(phg, move_order[i].node);
    });

    // apply vertex locking
    if (context.refinement.fm.vertex_locking > 0) {
      const MoveID start = context.refinement.fm.lock_moved_nodes ? 0 : b.best_index;
      tbb::parallel_for(start, numMoves, [&](const MoveID i) {
        if (!sharedData.moveTracker.isRebalancingMove(i + sharedData.moveTracker.firstMoveID)) {
          sharedData.lockVertexForNextRound(move_order[i].node, context);
        }
      });
    }

    sharedData.moveTracker.reset();

    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(gain_cache));
    return b.gain;
  }

  template<typename TypeTraits, typename GainTypes>
  void GlobalRollback<TypeTraits, GainTypes>::recalculateGains(PartitionedHypergraph& phg,
                                                               FMSharedData& sharedData) {
    GlobalMoveTracker& tracker = sharedData.moveTracker;

    auto recalculate_and_distribute_for_hyperedge = [&](const HyperedgeID e) {
      // TODO reduce .local() calls using the blocked_range interface
      auto& r = ets_recalc_data.local();

      // compute auxiliary data
      for (HypernodeID v : phg.pins(e)) {
        if (tracker.wasNodeMovedInThisRound(v)) {
          const MoveID m_id = tracker.moveOfNode[v];
          const Move& m = tracker.getMove(m_id);
          Rollback::updateMove(m_id, m, r);
          // no change for remaining pins!
        } else {
          Rollback::updateNonMovedPinInBlock(phg.partID(v), r);
        }
      }

      // distribute gains to pins
      for (HypernodeID v : phg.pins(e)) {
        if (tracker.wasNodeMovedInThisRound(v)) {
          const MoveID m_id = tracker.moveOfNode[v];
          Move& m = tracker.getMove(m_id);

          const HyperedgeWeight benefit = Rollback::benefit(phg, e, m_id, m, r);;
          const HyperedgeWeight penalty = Rollback::penalty(phg, e, m_id, m, r);

          if ( benefit > 0 ) {
            // increase gain of v by benefit
            __atomic_fetch_add(&m.gain, benefit, __ATOMIC_RELAXED);
          }

          if ( penalty > 0 ) {
            // decrease gain of v by penalty
            __atomic_fetch_sub(&m.gain, penalty, __ATOMIC_RELAXED);
          }
        }
      }

      if (context.partition.k <= static_cast<int>(2 * phg.edgeSize(e))) {
        // this branch is an optimization. in case it is cheaper to iterate over the parts, do that
        for (PartitionID i = 0; i < context.partition.k; ++i) {
          r[i].reset();
        }
      } else {
        for (HypernodeID v : phg.pins(e)) {
          if (tracker.wasNodeMovedInThisRound(v)) {
            const Move& m = tracker.getMove(tracker.moveOfNode[v]);
            r[m.from].reset();
            r[m.to].reset();
          } else {
            r[phg.partID(v)].reset();
          }
        }
      }
    };

    if (context.refinement.fm.iter_moves_on_recalc) {
      tbb::parallel_for(0U, sharedData.moveTracker.numPerformedMoves(), [&](const MoveID local_move_id) {
        const HypernodeID u = sharedData.moveTracker.moveOrder[local_move_id].node;
        if (tracker.wasNodeMovedInThisRound(u)) {
          for (HyperedgeID e : phg.incidentEdges(u)) {
            // test-and-set whether this is the first time this hyperedge is encountered
            uint32_t expected = last_recalc_round[phg.uniqueEdgeID(e)].load(std::memory_order_relaxed);
            if (expected < round && last_recalc_round[phg.uniqueEdgeID(e)].exchange(round, std::memory_order_acquire) == expected) {
              recalculate_and_distribute_for_hyperedge(e);
            }
          }
        }
      });

      // reset bits
      if (++round == std::numeric_limits<uint32_t>::max()) {
        // should never happen on practical inputs.
        last_recalc_round.assign(phg.initialNumEdges(), CAtomic<uint32_t>(0));
      }
    } else{
      tbb::parallel_for(0U, phg.initialNumEdges(), recalculate_and_distribute_for_hyperedge);
    }
  }

  template<typename TypeTraits, typename GainTypes>
  HyperedgeWeight GlobalRollback<TypeTraits, GainTypes>::revertToBestPrefixSequential(
    PartitionedHypergraph& phg,
    FMSharedData& sharedData,
    const vec<HypernodeWeight>&,
    const std::vector<HypernodeWeight>& maxPartWeights) {

    GlobalMoveTracker& tracker = sharedData.moveTracker;
    const MoveID numMoves = tracker.numPerformedMoves();
    const vec<Move>& move_order = tracker.moveOrder;

    // revert all moves
    tbb::parallel_for(0U, numMoves, [&](const MoveID localMoveID) {
      const Move& m = move_order[localMoveID];
      if (m.isValid()) {
        moveVertex(phg, m.node, m.to, m.from);
      }
    });


    size_t num_unbalanced_slots = 0;

    size_t overloaded = 0;
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      if (phg.partWeight(i) > maxPartWeights[i]) {
        overloaded++;
      }
    }

    // roll forward sequentially
    Gain best_gain = 0, gain_sum = 0;
    MoveID best_index = 0;
    for (MoveID localMoveID = 0; localMoveID < numMoves; ++localMoveID) {
      const Move& m = move_order[localMoveID];
      if (!m.isValid()) continue;

      Gain gain = 0;
      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        const HypernodeID pin_count_in_from_part_after = phg.pinCountInPart(e, m.from) - 1;
        const HypernodeID pin_count_in_to_part_after = phg.pinCountInPart(e, m.to) + 1;
        gain -= AttributedGains::gain(e, phg.edgeWeight(e), phg.edgeSize(e),
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      }
      gain_sum += gain;

      const bool from_overloaded = phg.partWeight(m.from) > maxPartWeights[m.from];
      const bool to_overloaded = phg.partWeight(m.to) > maxPartWeights[m.to];
      moveVertex(phg, m.node, m.from, m.to);
      if (from_overloaded && phg.partWeight(m.from) <= maxPartWeights[m.from]) {
        overloaded--;
      }
      if (!to_overloaded && phg.partWeight(m.to) > maxPartWeights[m.to]) {
        overloaded++;
      }

      if (overloaded > 0) {
        num_unbalanced_slots++;
      }

      if (overloaded == 0 && gain_sum > best_gain) {
        best_index = localMoveID + 1;
        best_gain = gain_sum;
      }
    }

    // revert rejected moves again
    tbb::parallel_for(best_index, numMoves, [&](const MoveID i) {
      const Move& m = move_order[i];
      if (m.isValid()) {
        moveVertex(phg, m.node, m.to, m.from);
      }
    });

    tbb::parallel_for(0U, numMoves, [&](const MoveID i) {
      gain_cache.recomputePenaltyTermEntry(phg, move_order[i].node);
    });

    tracker.reset();

    return best_gain;
  }


  template<typename TypeTraits, typename GainTypes>
  bool GlobalRollback<TypeTraits, GainTypes>::verifyGains(PartitionedHypergraph& phg,
                                                          FMSharedData& sharedData) {
    vec<Move>& move_order = sharedData.moveTracker.moveOrder;

    auto recompute_penalty_terms = [&] {
      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        gain_cache.recomputePenaltyTermEntry(phg, move_order[localMoveID].node);
      }
    };

    recompute_penalty_terms();
    phg.checkTrackedPartitionInformation(gain_cache);

    // revert all moves
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (m.isValid()) {
        moveVertex(phg, m.node, m.to, m.from);
      }
    }

    recompute_penalty_terms();

    // roll forward sequentially and check gains
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (!m.isValid())
        continue;

      Gain gain = 0;
      for (HyperedgeID e: phg.incidentEdges(m.node)) {
        const HypernodeID pin_count_in_from_part_after = phg.pinCountInPart(e, m.from) - 1;
        const HypernodeID pin_count_in_to_part_after = phg.pinCountInPart(e, m.to) + 1;
        gain -= AttributedGains::gain(e, phg.edgeWeight(e), phg.edgeSize(e),
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      }

      ASSERT(gain_cache.penaltyTerm(m.node, phg.partID(m.node)) == gain_cache.recomputePenaltyTerm(phg, m.node));
      ASSERT(gain_cache.benefitTerm(m.node, m.to) == gain_cache.recomputeBenefitTerm(phg, m.node, m.to));
      ASSERT(gain == gain_cache.gain(m.node, m.from, m.to));

      // const HyperedgeWeight objective_before_move =
      //   metrics::quality(phg, context, false);
      moveVertex(phg, m.node, m.from, m.to);
      // const HyperedgeWeight objective_after_move =
      //   metrics::quality(phg, context, false);

      // ASSERT(objective_after_move + gain == objective_before_move,
      //   V(gain) << V(m.gain) << V(objective_after_move) << V(objective_before_move));
      // ASSERT(objective_after_move + m.gain == objective_before_move,
      //   V(gain) << V(m.gain) << V(objective_after_move) << V(objective_before_move));
      ASSERT(gain == m.gain, V(gain) << V(m.gain));
      unused(gain); // unused(objective_before_move); unused(objective_after_move);  // for release mode
    }

    recompute_penalty_terms();
    return true;
  }

  namespace {
  #define GLOBAL_ROLLBACK(X, Y) GlobalRollback<X, Y>
  }

  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(GLOBAL_ROLLBACK)
}
