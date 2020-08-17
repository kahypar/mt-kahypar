#include "global_rollback.h"

namespace mt_kahypar {

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
    return revertToBestPrefixParallel(phg, sharedData, partWeights, maxPartWeights);
  }

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
            phg.changeNodePartFullUpdate(m.node, m.to, m.from);
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
        phg.changeNodePartFullUpdate(m.node, m.to, m.from);
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

    HEAVY_REFINEMENT_ASSERT(verifyGains(phg, sharedData));
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

  bool GlobalRollback::verifyGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    vec<Move>& move_order = sharedData.moveTracker.moveOrder;
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }
    phg.checkTrackedPartitionInformation();

    // revert all moves
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (sharedData.moveTracker.isMoveStillValid(m))
        phg.changeNodePartFullUpdate(m.node, m.to, m.from);
    }

    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      if (sharedData.moveTracker.isMoveStillValid(move_order[localMoveID]))
        phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }

    // roll forward sequentially and check gains
    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
      if (!sharedData.moveTracker.isMoveStillValid(m))
        continue;

      const Gain estimated_gain = phg.km1Gain(m.node, m.from, m.to);
      ASSERT(phg.moveFromBenefit(m.node) == phg.moveFromBenefitRecomputed(m.node));
      ASSERT(phg.moveToPenalty(m.node, m.to) == phg.moveToPenaltyRecomputed(m.node, m.to));
      const HyperedgeWeight km1_before_move = metrics::km1(phg, false);
      phg.changeNodePartFullUpdate(m.node, m.from, m.to);
      const HyperedgeWeight km1_after_move = metrics::km1(phg, false);
      ASSERT(km1_after_move + estimated_gain == km1_before_move);
      ASSERT(km1_after_move + m.gain == km1_before_move);
      ASSERT(estimated_gain == m.gain);
    }

    for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    }
    return true;
  }

  void GlobalRollback::resetStoredMoveIDs() {
    tbb::parallel_invoke([&] {
                           tbb::parallel_for_each(last_move_out, [](auto& x) { x.store(0, std::memory_order_relaxed); });
                         }, [&] {
                           tbb::parallel_for_each(first_move_in, [](auto& x) { x.store(0, std::memory_order_relaxed); });
                         }
    );
  }

  void GlobalRollback::setRemainingOriginalPins(PartitionedHypergraph& phg) {
    if (phg.hypergraph().maxEdgeSize() > 2) {
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

}