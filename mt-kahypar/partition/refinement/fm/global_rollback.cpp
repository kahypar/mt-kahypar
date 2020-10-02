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
        if (m.gain != invalidGain) {  // skip locally reverted moves
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

        if (m.gain != invalidGain) {  // skip locally reverted moves
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
    if (maxPartWeightScaling == 0.0) {
      for (PartitionID i = 0; i < numParts; ++i) {
        maxPartWeights[i] = std::numeric_limits<HypernodeWeight>::max();
      }
    } else {
      for (PartitionID i = 0; i < numParts; ++i) {
        maxPartWeights[i] *= ( 1.0 + context.partition.epsilon * maxPartWeightScaling );
      }
    }

    if (context.refinement.fm.revert_parallel) {
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
    tbb::parallel_scan(tbb::blocked_range<MoveID>(0, numMoves), s);
    BalanceAndBestIndexScan::Prefix b = s.finalize(partWeights);
    timer.stop_timer("find_best_prefix_and_balance");

    timer.start_timer("revert_and_rem_orig_pin_updates", "Revert Moves and apply updates");

    if (phg.hypergraph().maxEdgeSize() > 2) {
      // revert and apply updates (full version for hypergraphs)
      tbb::parallel_invoke([&] {
        // revert rejected moves
        tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
          const Move& m = move_order[moveID];
          if (sharedData.moveTracker.isMoveStillValid(m)) {
            moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
            for (HyperedgeID e : phg.incidentEdges(m.node)) {
              if (phg.edgeSize(e) > 2) {
                remaining_original_pins[size_t(phg.nonGraphEdgeID(e)) * numParts + m.from].fetch_add(1, std::memory_order_relaxed);
              }
            }
          }
        });
      }, [&] {
        // apply updates to remaining original pins
        tbb::parallel_for(MoveID(0), b.best_index, [&](const MoveID moveID) {
          const Move& m = move_order[moveID];
          if (sharedData.moveTracker.isMoveStillValid(m)) {
            for (HyperedgeID e : phg.incidentEdges(move_order[moveID].node)) {
              if (phg.edgeSize(e) > 2) {
                remaining_original_pins[size_t(phg.nonGraphEdgeID(e)) * numParts + m.to].fetch_add(1, std::memory_order_relaxed);
              }
            }
          }
        });
      } );

    } else {
      // faster special case for graphs
      tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
        const Move& m = move_order[moveID];
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      });
    }

    timer.stop_timer("revert_and_rem_orig_pin_updates");

    timer.start_timer("recompute_move_from_benefits", "Recompute Move-From Benefits");
    // recompute moveFromBenefit values since they are potentially invalid
    tbb::parallel_for(MoveID(0), numMoves, [&](MoveID localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    });
    timer.stop_timer("recompute_move_from_benefits");

    if (sharedData.moveTracker.reset()) {
      resetStoredMoveIDs();
    }

    HEAVY_REFINEMENT_ASSERT(verifyPinCountMatches(phg));
    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
    return b.gain;
  }


  void GlobalRollback::recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    GlobalMoveTracker& tracker = sharedData.moveTracker;
    vec<Move>& move_order = tracker.moveOrder;
    const MoveID numMoves = tracker.numPerformedMoves();
    const MoveID firstMoveID = tracker.firstMoveID;

    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("move_id_flagging", "Move Flagging");

    if (phg.hypergraph().maxEdgeSize() > 2) {
      tbb::parallel_for(0U, numMoves, [&](MoveID localMoveID) {
        const Move& m = move_order[localMoveID];
        if (!sharedData.moveTracker.isMoveStillValid(m)) return;

        const MoveID globalMoveID = localMoveID + firstMoveID;
        for (HyperedgeID e_global : phg.incidentEdges(m.node)) {
          if (phg.edgeSize(e_global) > 2) {
            const HyperedgeID e = phg.nonGraphEdgeID(e_global);
            CAtomic<MoveID>& fmi = first_move_in[size_t(e) * numParts + m.to];
            MoveID expected = fmi.load(std::memory_order_acq_rel);
            // first_move_in = min(first_move_in, this_move)
            while ((tracker.isMoveStale(expected) || expected > globalMoveID)
                   && !fmi.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

            CAtomic<MoveID>& lmo = last_move_out[size_t(e) * numParts + m.from];
            expected = lmo.load(std::memory_order_acq_rel);
            // last_move_out = max(last_move_out, this_move)
            while (expected < globalMoveID && !lmo.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

            remaining_original_pins[size_t(e) * numParts + m.from].fetch_sub(1, std::memory_order_relaxed);
          }
        }
      });
    }

    timer.stop_timer("move_id_flagging");
    timer.start_timer("gain_recalculation", "Recalculate Gains");

    tbb::parallel_for(0U, numMoves, [&](MoveID localMoveID) {
      Move& m = move_order[localMoveID];
      if (!tracker.isMoveStillValid(m)) {
        return;
      }

      MoveID globalMoveID = firstMoveID + localMoveID;
      Gain gain = 0;
      for (HyperedgeID e_global : phg.incidentEdges(m.node)) {
        const HyperedgeWeight edgeWeight = phg.edgeWeight(e_global);

        if (phg.edgeSize(e_global) > 2) {
          const HyperedgeID e = phg.nonGraphEdgeID(e_global);

          const MoveID firstInFrom = firstMoveIn(e, m.from);
          if (remainingPinsFromBeginningOfMovePhase(e, m.from) == 0 && lastMoveOut(e, m.from) == globalMoveID
              && (firstInFrom > globalMoveID || firstInFrom < firstMoveID)) {
            gain += edgeWeight;
          }

          const MoveID lastOutTo = lastMoveOut(e, m.to);
          if (remainingPinsFromBeginningOfMovePhase(e, m.to) == 0
              && firstMoveIn(e, m.to) == globalMoveID && lastOutTo < globalMoveID) {
            gain -= edgeWeight;
          }

        } else {

          const HypernodeID u = m.node;
          const HypernodeID v = phg.graphEdgeHead(e_global, u);
          const MoveID move_id_of_v = tracker.moveOfNode[v];
          const bool v_moved = !tracker.isMoveStale(move_id_of_v) && tracker.isMoveStillValid(move_id_of_v);

          PartitionID pv;
          if (v_moved) {
            const MoveID local_move_id_of_v = move_id_of_v - firstMoveID;
            const Move& move_of_v = move_order[local_move_id_of_v];
            if (local_move_id_of_v < localMoveID) {
              pv = move_of_v.to;
            } else {
              pv = move_of_v.from;
            }
          } else {
            pv = phg.partID(v);
          }

          if (pv == m.to) {
            gain += edgeWeight;
          } else if (pv == m.from) {
            gain -= edgeWeight;
          }

        }

      }

      m.gain = gain;
    });

    timer.stop_timer("gain_recalculation");

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
      if (tracker.isMoveStillValid(m)) {
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
      }
    });


    size_t num_unbalanced_slots = 0;

    size_t overloaded = 0;
    for (PartitionID i = 0; i < numParts; ++i) {
      if (phg.partWeight(i) > maxPartWeights[i]) {
        overloaded++;
      }
    }

    // roll forward sequentially
    Gain best_gain = 0, gain_sum = 0;
    MoveID best_index = 0;
    for (MoveID localMoveID = 0; localMoveID < numMoves; ++localMoveID) {
      const Move& m = move_order[localMoveID];
      if (!tracker.isMoveStillValid(m)) continue;

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
      if (tracker.isMoveStillValid(m)) {
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


  bool GlobalRollback::verifyPinCountMatches(const PartitionedHypergraph& phg) const {
    for ( const HyperedgeID& he : phg.edges() ) {
      if (phg.edgeSize(he) == 2) {
        continue;
      }
      for ( PartitionID block = 0; block < numParts; ++block ) {
        if ( phg.pinCountInPart(he, block) != remaining_original_pins[phg.nonGraphEdgeID(he) * numParts + block] ) {
          LOG << "Pin Count In Part does not match remaining original pins" << V(he) << V(block)
              << V(phg.pinCountInPart(he, block)) << V(remaining_original_pins[phg.nonGraphEdgeID(he) * numParts + block]);
          return false;
        }
      }
    }
    return true;
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
      if (sharedData.moveTracker.isMoveStillValid(m))
        moveVertex<update_gain_cache>(phg, m.node, m.to, m.from);
    }

    recompute_move_from_benefits();

    // roll forward sequentially and check gains
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (!sharedData.moveTracker.isMoveStillValid(m))
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

  void GlobalRollback::resetStoredMoveIDs() {
    if (context.refinement.fm.revert_parallel) {
      tbb::parallel_invoke(
              [&] {
                tbb::parallel_for(0UL, last_move_out.size(),
                                  [&](size_t i) { last_move_out[i].store(0, std::memory_order_relaxed); });
              },
              [&] {
                tbb::parallel_for(0UL, first_move_in.size(),
                                  [&](size_t i) { first_move_in[i].store(0, std::memory_order_relaxed); });
              }
      );
    }
  }

  void GlobalRollback::setRemainingOriginalPins(PartitionedHypergraph& phg) {
    if (phg.hypergraph().maxEdgeSize() > 2 && context.refinement.fm.revert_parallel) {
      phg.doParallelForAllEdges([&](const HyperedgeID& he) {
        if (phg.edgeSize(he) > 2) {
          const HyperedgeID he_local = phg.nonGraphEdgeID(he);
          for ( PartitionID block = 0; block < numParts; ++block ) {
            remaining_original_pins[size_t(he_local) * numParts + block]
                    .store(phg.pinCountInPart(he, block), std::memory_order_relaxed);
          }
        }
      });
    }
  }

  void GlobalRollback::memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* global_rollback_node = parent->addChild("Global Rollback");
    utils::MemoryTreeNode* remaining_original_pins_node = global_rollback_node->addChild("Remaining Original Pins");
    remaining_original_pins_node->updateSize(remaining_original_pins.size() * sizeof(HypernodeID));
    utils::MemoryTreeNode* first_move_in_node = global_rollback_node->addChild("First Move In");
    first_move_in_node->updateSize(first_move_in.size() * sizeof(MoveID));
    utils::MemoryTreeNode* last_move_out_node = global_rollback_node->addChild("Last Move Out");
    last_move_out_node->updateSize(last_move_out.size() * sizeof(MoveID));
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