/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/fm/global_rollback.h"

#include "tbb/parallel_scan.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

  struct BalanceAndBestIndexScan {
    const PartitionedHypergraph& phg;
    const vec<Move>& moves;

    struct Prefix {
      Gain gain = 0;                           /** gain when using valid moves up to best_index */
      MoveID best_index = 0;                   /** local ID of first move to revert */
      HypernodeWeight heaviest_weight =
              std::numeric_limits<HypernodeWeight>::max();   /** weight of the heaviest part */

      bool operator<(const Prefix& o) const {
        return gain > o.gain ||
               (gain == o.gain && std::tie(heaviest_weight, best_index) < std::tie(o.heaviest_weight, o.best_index));
      }
    };
    std::shared_ptr< tbb::enumerable_thread_specific<Prefix> > local_best;

    Gain gain_sum = 0;

    vec<HypernodeWeight> part_weights;
    const std::vector<HypernodeWeight>& max_part_weights;

    BalanceAndBestIndexScan(BalanceAndBestIndexScan& b, tbb::split) :
            phg(b.phg),
            moves(b.moves),
            local_best(b.local_best),
            gain_sum(0),
            part_weights(b.part_weights.size(), 0),
            max_part_weights(b.max_part_weights) { }


    BalanceAndBestIndexScan(const PartitionedHypergraph& phg,
                            const vec<Move>& moves,
                            const vec<HypernodeWeight>& part_weights,
                            const std::vector<HypernodeWeight>& max_part_weights) :
            phg(phg),
            moves(moves),
            local_best(std::make_shared< tbb::enumerable_thread_specific<Prefix> >()),
            part_weights(part_weights),
            max_part_weights(max_part_weights)
    {
    }


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

    // subranges a | b | c | d . assuming this ran pre_scan on c,
    // then lhs ran pre_scan on b and final_scan of this will be on d
    void reverse_join(BalanceAndBestIndexScan& lhs) {
      for (size_t i = 0; i < part_weights.size(); ++i) {
        part_weights[i] += lhs.part_weights[i];
      }
      gain_sum += lhs.gain_sum;
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
            Prefix new_prefix = { gain_sum, i + 1, *std::max_element(part_weights.begin(), part_weights.end()) };
            current = std::min(current, new_prefix);
          }
        }
      }

      if (current.best_index != 0) {
        Prefix& lb = local_best->local();
        lb = std::min(lb, current);
      }
    }

    void assign(BalanceAndBestIndexScan& b) {
      gain_sum = b.gain_sum;
    }

    Prefix finalize(const vec<HypernodeWeight>& initial_part_weights) {
      Prefix res { 0, 0, *std::max_element(initial_part_weights.begin(), initial_part_weights.end()) };
      for (const Prefix& x : *local_best) {
        res = std::min(res, x);
      }
      return res;
    }
  };

  template<bool update_gain_cache>
  HyperedgeWeight GlobalRollback::revertToBestPrefix(
          PartitionedHypergraph& phg, FMSharedData& sharedData,
          const vec<HypernodeWeight>& partWeights) {
    std::vector<HypernodeWeight> maxPartWeights = context.partition.perfect_balance_part_weights;
    if (max_part_weight_scaling == 0.0) {
      for (PartitionID i = 0; i < num_parts; ++i) {
        maxPartWeights[i] = std::numeric_limits<HypernodeWeight>::max();
      }
    } else {
      for (PartitionID i = 0; i < num_parts; ++i) {
        maxPartWeights[i] *= ( 1.0 + context.partition.epsilon * max_part_weight_scaling );
      }
    }

    if (context.refinement.fm.rollback_parallel) {
      return revertToBestPrefixParallel<update_gain_cache>(phg, sharedData, partWeights, maxPartWeights);
    } else {
      return revertToBestPrefixSequential<update_gain_cache>(phg, sharedData, partWeights, maxPartWeights);
    }

  }


  template<bool update_gain_cache>
  HyperedgeWeight GlobalRollback::revertToBestPrefixParallel(
          PartitionedHypergraph& phg, FMSharedData& sharedData,
          const vec<HypernodeWeight>& partWeights, const std::vector<HypernodeWeight>& maxPartWeights) {
    const MoveID numMoves = sharedData.moveTracker.numPerformedMoves();
    if (numMoves == 0) return 0;

    const vec<Move>& move_order = sharedData.moveTracker.moveOrder;
    utils::Timer& timer = utils::Timer::instance();

    timer.start_timer("recalculate_gains", "Recalculate Gains");
    recalculateGains(phg, sharedData);
    timer.stop_timer("recalculate_gains");
    HEAVY_REFINEMENT_ASSERT(verifyGains<update_gain_cache>(phg, sharedData));

    timer.start_timer("find_best_prefix_and_balance", "Find Best Balanced Prefix");
    BalanceAndBestIndexScan s(phg, move_order, partWeights, maxPartWeights);
    // TODO set grain size in blocked_range? to avoid too many copies of part weights array. experiment with different values
    tbb::parallel_scan(tbb::blocked_range<MoveID>(0, numMoves), s);
    BalanceAndBestIndexScan::Prefix b = s.finalize(partWeights);
    timer.stop_timer("find_best_prefix_and_balance");

    timer.start_timer("revert", "Revert Moves");
    tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
      const Move& m = move_order[moveID];
      if (m.isValid()) {
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      }
    });
    timer.stop_timer("revert");

    timer.start_timer("recompute_move_from_benefits", "Recompute Move-From Benefits");
    // recompute moveFromBenefit values since they are potentially invalid
    tbb::parallel_for(MoveID(0), numMoves, [&](MoveID localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    });
    timer.stop_timer("recompute_move_from_benefits");

    sharedData.moveTracker.reset();

    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
    return b.gain;
  }


  void GlobalRollback::recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    GlobalMoveTracker& tracker = sharedData.moveTracker;

    auto recalculate_and_distribute_for_hyperedge = [&](const HyperedgeID e) {
      // TODO reduce .local() calls using the blocked_range interface
      auto& r = ets_recalc_data.local();

      // compute auxiliary data
      for (HypernodeID v : phg.pins(e)) {
        if (tracker.wasNodeMovedInThisRound(v)) {
          const MoveID m_id = tracker.moveOfNode[v];
          const Move& m = tracker.getMove(m_id);
          r[m.to].first_in = std::min(r[m.to].first_in, m_id);
          r[m.from].last_out = std::max(r[m.from].last_out, m_id);
          // no change for remaining pins!
        } else {
          r[phg.partID(v)].remaining_pins++;
        }
      }

      // distribute gains to pins
      const HyperedgeWeight we = phg.edgeWeight(e);
      for (HypernodeID v : phg.pins(e)) {
        if (tracker.wasNodeMovedInThisRound(v)) {
          const MoveID m_id = tracker.moveOfNode[v];
          Move& m = tracker.getMove(m_id);

          const bool benefit = r[m.from].last_out == m_id && r[m.from].first_in > m_id && r[m.from].remaining_pins == 0;
          const bool penalty = r[m.to].first_in == m_id && r[m.to].last_out < m_id && r[m.to].remaining_pins == 0;

          if (benefit && !penalty) {    // only apply update if they're mutually exclusive
            // increase gain of v by w(e)
            __atomic_fetch_add(&m.gain, we, __ATOMIC_RELAXED);
          }

          if (!benefit && penalty) {
            // decrease gain of v by w(e)
            __atomic_fetch_sub(&m.gain, we, __ATOMIC_RELAXED);
          }
        }
      }

      if (num_parts <= static_cast<int>(2 * phg.edgeSize(e))) {
        // this branch is an optimization. in case it is cheaper to iterate over the parts, do that
        for (PartitionID i = 0; i < num_parts; ++i) {
          r[i] = RecalculationData();
        }
      } else {
        for (HypernodeID v : phg.pins(e)) {
          if (tracker.wasNodeMovedInThisRound(v)) {
            const Move& m = tracker.getMove(tracker.moveOfNode[v]);
            r[m.from] = RecalculationData();
            r[m.to] = RecalculationData();
          } else {
            r[phg.partID(v)].remaining_pins = 0;
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
            uint32_t expected = last_recalc_round[e].load(std::memory_order_relaxed);
            if (expected < round && last_recalc_round[e].exchange(round, std::memory_order_acquire) == expected) {
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

  template<bool update_gain_cache>
  HyperedgeWeight GlobalRollback::revertToBestPrefixSequential(PartitionedHypergraph& phg,
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
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      }
    });


    size_t num_unbalanced_slots = 0;

    size_t overloaded = 0;
    for (PartitionID i = 0; i < num_parts; ++i) {
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
        if (phg.pinCountInPart(e, m.from) == 1) {
          gain += phg.edgeWeight(e);
        }
        if (phg.pinCountInPart(e, m.to) == 0) {
          gain -= phg.edgeWeight(e);
        }
      }
      gain_sum += gain;

      const bool from_overloaded = phg.partWeight(m.from) > maxPartWeights[m.from];
      const bool to_overloaded = phg.partWeight(m.to) > maxPartWeights[m.to];
      moveVertex<update_gain_cache>(phg, m.node, m.from, m.to);
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
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      }
    });

    if constexpr (update_gain_cache) {
      tbb::parallel_for(0U, numMoves, [&](const MoveID i) {
        phg.recomputeMoveFromBenefit(move_order[i].node);
      });
    }

    tracker.reset();

    return best_gain;
  }



  template<bool update_gain_cache>
  bool GlobalRollback::verifyGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    vec<Move>& move_order = sharedData.moveTracker.moveOrder;

    auto recompute_move_from_benefits = [&] {
      if constexpr (update_gain_cache) {
        for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
          phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
        }
      }
    };

    recompute_move_from_benefits();
    phg.checkTrackedPartitionInformation();

    // revert all moves
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (m.isValid()) {
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      }
    }

    recompute_move_from_benefits();

    // roll forward sequentially and check gains
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (!m.isValid())
        continue;

      Gain gain = 0;
      for (HyperedgeID e: phg.incidentEdges(m.node)) {
        if (phg.pinCountInPart(e, m.from) == 1) gain += phg.edgeWeight(e);
        if (phg.pinCountInPart(e, m.to) == 0) gain -= phg.edgeWeight(e);
      }

      if constexpr (update_gain_cache) {
        const Gain gain_from_cache = phg.km1Gain(m.node, m.from, m.to); unused(gain_from_cache);
        ASSERT(phg.moveFromBenefit(m.node) == phg.moveFromBenefitRecomputed(m.node));
        ASSERT(phg.moveToPenalty(m.node, m.to) == phg.moveToPenaltyRecomputed(m.node, m.to));
        ASSERT(gain == gain_from_cache);
      }

      const HyperedgeWeight km1_before_move = metrics::km1(phg, false);
      moveVertex<update_gain_cache>(phg, m.node, m.from, m.to);
      const HyperedgeWeight km1_after_move = metrics::km1(phg, false);

      ASSERT(km1_after_move + gain == km1_before_move);
      ASSERT(km1_after_move + m.gain == km1_before_move);
      ASSERT(gain == m.gain);
      unused(gain); unused(km1_before_move); unused(km1_after_move);  // for release mode
    }

    recompute_move_from_benefits();
    return true;
  }


  // template instantiations
  template HyperedgeWeight GlobalRollback::revertToBestPrefix<false>
          (PartitionedHypergraph& , FMSharedData& , const vec<HypernodeWeight>& );

  template HyperedgeWeight GlobalRollback::revertToBestPrefix<true>
          (PartitionedHypergraph& , FMSharedData& , const vec<HypernodeWeight>& );

  template bool GlobalRollback::verifyGains<false>
          (PartitionedHypergraph& , FMSharedData& );

  template bool GlobalRollback::verifyGains<true>
          (PartitionedHypergraph& , FMSharedData& );

}