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

#include "mt-kahypar/partition/refinement/fm/async_kway_fm_core.h"

namespace mt_kahypar {


  template<typename FMStrategy>
  Gain AsyncKWayFM<FMStrategy>::findMoves(PartitionedHypergraph& phg, const vec<HypernodeID>& refinement_nodes) {
    ASSERT(contraction_group_id != ds::invalidGroupID, "ContractionGroupID (Owner-ID) for locking is invalid.");

    _km1_delta = 0;
    localMoves.clear();
    attempted_to_move.reset();
    runStats.clear();
    num_nodes_in_pq = 0;


    for (const auto& seedNode : refinement_nodes ) {
      if (sharedData.nodeTracker.tryAcquireNode(seedNode, contraction_group_id, true)) {
        // Third parameter is only needed by gain cache on demand strategy which currently does not work with async so
        // just give any value
        fm_strategy.insertIntoPQ(phg, seedNode, ds::invalidGroupID);
        ++num_nodes_in_pq;
      }
    }

    if (runStats.pushes > 0) {
      deltaPhg.clear();
      deltaPhg.setPartitionedHypergraph(&phg);
      deltaPhg.markActive();
      internalFindMoves(phg);
//      deltaPhg.clear();
      deltaPhg.markInactive();
    }

    contraction_group_id = ds::invalidGroupID;
    return _km1_delta;
  }

  template<typename Partition>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight>
  heaviestPartAndWeight(const Partition& partition) {
    PartitionID p = kInvalidPartition;
    HypernodeWeight w = std::numeric_limits<HypernodeWeight>::min();
    for (PartitionID i = 0; i < partition.k(); ++i) {
      if (partition.partWeight(i) > w) {
        w = partition.partWeight(i);
        p = i;
      }
    }
    return std::make_pair(p, w);
  }

  template<typename FMStrategy>
  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void AsyncKWayFM<FMStrategy>::acquireOrUpdateNeighbors(PHG& phg, const Move& move) {
    // Note: In theory we should acquire/update all neighbors. It just turned out that this works fine
    // Actually: only vertices incident to edges with gain changes can become new boundary vertices.
    // Vertices that already were boundary vertices, can still be considered later since they are in the task queue
    // --> actually not that bad
    for (HyperedgeID e : edgesWithGainChanges) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
        auto pin_range = phg.pins(e);
        for (auto it = pin_range.begin(); it != pin_range.end(); ++it) {
          const HypernodeID v = *it;
          // Skip nodes that are not enabled to prevent inserting nodes that are in the intermediate state of
          // uncontraction where they are in the incidence array but not yet enabled
          if (!phg.nodeIsEnabled(v)) continue;
          if (!attempted_to_move[v] && neighborDeduplicator[v] != deduplicationTime) {
            SearchID searchOfV = sharedData.nodeTracker.owner(v);
            if (searchOfV == contraction_group_id) {
              fm_strategy.updateGain(phg, v, move);
            } else if (num_nodes_in_pq <= max_num_nodes_in_pq && sharedData.nodeTracker.tryAcquireNode(v,
                                                                                                       contraction_group_id,
                                                                                                       true)) {
              fm_strategy.insertIntoPQ(phg, v, searchOfV);
              ++num_nodes_in_pq;
            }
            neighborDeduplicator[v] = deduplicationTime;
          }
        }
      }
    }
    edgesWithGainChanges.clear();

    if (++deduplicationTime == 0) {
      neighborDeduplicator.assign(neighborDeduplicator.size(), 0);
      deduplicationTime = 1;
    }
  }

//  template<typename FMStrategy, typename PHG>
//  struct DeltaGainUpdatesFuncWrapper {
//
//        explicit DeltaGainUpdatesFuncWrapper(FMStrategy& strategy, PHG& phg) : _strategy(strategy), _phg(phg) {}
//
//        FMStrategy& _strategy;
//        PHG& _phg;
//
//        template<typename PinIteratorT>
//        void operator() (const HyperedgeID he, const HyperedgeWeight edge_weight, IteratorRange<PinIteratorT> pins,
//                         const PartitionID from, const HypernodeID pin_count_in_from_part_after,
//                         const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
//          _strategy.deltaGainUpdates(_phg, he, edge_weight, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
//        }
//      };

  template<typename FMStrategy>
  void AsyncKWayFM<FMStrategy>::internalFindMoves(PartitionedHypergraph& phg) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    bool foundAtLeastOneGoodPrefix = false;

    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        edgesWithGainChanges.push_back(he);
      }

//      _km1_delta += km1Delta(he, edge_weight, edge_size,pin_count_in_from_part_after, pin_count_in_to_part_after);

      fm_strategy.deltaGainUpdates(deltaPhg, he, edge_weight, deltaPhg.pins(he), move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
    };

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop() && runStats.moves < max_num_moves) {

      //Attempt to find the next move. If none can be found then stop this search.
      bool move_found = fm_strategy.findNextMove(deltaPhg, move);
//      bool move_found = fm_strategy.findNextMove(phg, move);
      if (!move_found) break;
      --num_nodes_in_pq;

      ASSERT(move.isValid());

      // Mark this node to prevent attempting to move it again
      attempted_to_move.set(move.node, true);

      bool moved = false;
      if (move.to != kInvalidPartition) {
        heaviestPartWeight = heaviestPartAndWeight(deltaPhg).second;
        fromWeight = deltaPhg.partWeight(move.from);
        toWeight = deltaPhg.partWeight(move.to);

        lockUncontractionLockWithWaiting(move.node);

        moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
                                          context.partition.max_part_weights[move.to], delta_func, true);

        if (!moved) uncontraction_locks->strongReleaseLock(move.node, contraction_group_id);

//        heaviestPartWeight = heaviestPartAndWeight(phg).second;
//        fromWeight = phg.partWeight(move.from);
//        toWeight = phg.partWeight(move.to);
//
//        lockUncontractionLockWithWaiting(move.node);
//
//        moved = phg.changeNodePart(move.node, move.from, move.to,
//                                          context.partition.max_part_weights[move.to], []{}, delta_func, DeltaGainUpdatesFuncWrapper<FMStrategy, PartitionedHypergraph>(fm_strategy, phg), true);
//
//        if (!moved) uncontraction_locks->strongReleaseLock(move.node, contraction_group_id);
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += move.gain;
        localMoves.emplace_back(move);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          foundAtLeastOneGoodPrefix = true;
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();

          applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, true);
          bestImprovementIndex = 0;
          releaseMoveLocksForLocalMovedNodes();
          localMoves.clear();
          deltaPhg.clear();   // clear hashtables, save memory :)
        }

        acquireOrUpdateNeighbors(deltaPhg, move);
//        acquireOrUpdateNeighbors(phg, move);

      } else {
        // Release move lock on node of move if the move does not go through
        sharedData.nodeTracker.releaseNode(move.node, contraction_group_id);
      }
    }

//    std::cout << "Best improvement: " << bestImprovement << ", Index: " << bestImprovementIndex << std::endl;

    std::tie(bestImprovement, bestImprovementIndex) =
        applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, false);

//    revertToBestLocalPrefix(phg, bestImprovementIndex);

    runStats.estimated_improvement = bestImprovement;
    releaseMoveLocksForLocalMovedNodes();
    fm_strategy.clearPQs(bestImprovementIndex, true);
    num_nodes_in_pq = 0;

    sharedData.total_moves.fetch_add(runStats.moves, std::memory_order_relaxed);
    sharedData.total_reverts.fetch_add(runStats.local_reverts, std::memory_order_relaxed);
    sharedData.total_find_moves_calls.fetch_add(1, std::memory_order_relaxed);
    if (foundAtLeastOneGoodPrefix) sharedData.find_moves_calls_with_good_prefix.fetch_add(1, std::memory_order_relaxed);
    sharedData.find_move_retries.fetch_add(runStats.retries, std::memory_order_relaxed);
    sharedData.total_pushes_pos_gain.fetch_add(runStats.pushes_with_pos_gain, std::memory_order_relaxed);
    sharedData.total_pushes_non_pos_gain.fetch_add(runStats.pushes_with_non_pos_gain, std::memory_order_relaxed);

    runStats.merge(stats);
  }


  template<typename FMStrategy>
  std::pair<Gain, size_t> AsyncKWayFM<FMStrategy>::applyBestLocalPrefixToSharedPartition(
          PartitionedHypergraph& phg,
          const size_t best_index_locally_observed,
          const Gain best_improvement_locally_observed,
          bool apply_all_moves) {

    Gain improvement_from_attributed_gains = 0;
    Gain attributed_gain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      attributed_gain += km1Delta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain best_improvement_from_attributed_gains = 0;
    size_t best_index_from_attributed_gains = 0;
    for (size_t i = 0; i < best_index_locally_observed; ++i) {
      ASSERT(i < localMoves.size());
      Move& local_move = localMoves[i];
      attributed_gain = 0;

      // Prevent concurrent uncontractions for call to changeNodePart() with uncontraction lock
//      lockUncontractionLockWithWaiting(local_move.node);

      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(local_move.node, local_move.from, local_move.to,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [&] {},
                                              delta_gain_func, true);
      } else {
        phg.changeNodePart(local_move.node, local_move.from, local_move.to,
                           std::numeric_limits<HypernodeWeight>::max(),
                           [&] {},
                           delta_gain_func, true);
      }

//      uncontraction_locks->strongReleaseLock(local_move.node, contraction_group_id);

      _km1_delta += attributed_gain;
      attributed_gain = -attributed_gain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      improvement_from_attributed_gains += attributed_gain;
      if (improvement_from_attributed_gains >= best_improvement_from_attributed_gains) {
        best_improvement_from_attributed_gains = improvement_from_attributed_gains;
        best_index_from_attributed_gains = i;
      }
    }

    runStats.local_reverts += localMoves.size() - best_index_locally_observed;
    if (!apply_all_moves && best_index_from_attributed_gains != best_index_locally_observed) {
      runStats.best_prefix_mismatch++;
    }

    std::pair<Gain, size_t> result;

    // kind of double rollback, if attributed gains say we overall made things worse
    if (!apply_all_moves && improvement_from_attributed_gains < 0) {
      // always using the if-branch gave similar results
      runStats.local_reverts += best_index_locally_observed - best_index_from_attributed_gains + 1;
      for (size_t i = best_index_from_attributed_gains + 1; i < best_index_locally_observed; ++i) {
        Move& m = localMoves[i];
        attributed_gain = 0;

        // Prevent concurrent uncontractions for call to changeNodePart() with uncontraction lock
//        lockUncontractionLockWithWaiting(m.node);

        if constexpr (FMStrategy::uses_gain_cache) {
          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from,
                                                std::numeric_limits<HypernodeWeight>::max(),
                                                [&] {},
                                                delta_gain_func, true);
        } else {
          phg.changeNodePart(m.node, m.to, m.from,
                             std::numeric_limits<HypernodeWeight>::max(),
                             [&] {},
                             delta_gain_func, true);
        }

//        uncontraction_locks->strongReleaseLock(m.node, contraction_group_id);

        _km1_delta += attributed_gain;

        m.invalidate();
      }
      result = std::make_pair(best_improvement_from_attributed_gains, best_index_from_attributed_gains);
    } else {
      result = std::make_pair(best_improvement_locally_observed, best_index_locally_observed);
    }

    // Release uncontraction locks, i.e. allow the nodes to be uncontracted again once the best prefix has been applied
    for (auto & local_move : localMoves) {
      uncontraction_locks->strongReleaseLock(local_move.node, contraction_group_id);
    }

    return result;
  }

//    template<typename FMStrategy>
//    void AsyncKWayFM<FMStrategy>::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {
//
//      auto delta_func = [&](const HyperedgeID he,
//                            const HyperedgeWeight edge_weight,
//                            const HypernodeID edge_size,
//                            const HypernodeID pin_count_in_from_part_after,
//                            const HypernodeID pin_count_in_to_part_after) {
//          _km1_delta += km1Delta(he, edge_weight, edge_size,pin_count_in_from_part_after, pin_count_in_to_part_after);
//      };
//
//      runStats.local_reverts += localMoves.size() - bestGainIndex;
//      while (localMoves.size() > bestGainIndex) {
//        Move& m = localMoves.back();
//        if constexpr (FMStrategy::uses_gain_cache) {
//          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from,context.partition.max_part_weights[m.from],[]{}, delta_func, true);
//        } else {
//          phg.changeNodePart(m.node, m.to, m.from,context.partition.max_part_weights[m.from],[]{}, delta_func, true);
//        }
//        m.invalidate();
//
//        sharedData.nodeTracker.releaseNode(m.node, contraction_group_id);
//        uncontraction_locks->strongReleaseLock(m.node, contraction_group_id);
//
//        localMoves.pop_back();
//      }
//
//     Release uncontraction locks, i.e. allow the nodes to be uncontracted again once the best prefix has been applied
//    for (size_t i = 0; i < localMoves.size(); ++i) {
//    uncontraction_locks->strongReleaseLock(localMoves[i].node, contraction_group_id);
//}
//
//    }

  template<typename FMStrategy>
  void AsyncKWayFM<FMStrategy>::memoryConsumption(utils::MemoryTreeNode *parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode *localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode *deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(neighborDeduplicator.capacity() * sizeof(HypernodeID));
    utils::MemoryTreeNode *edges_to_activate_node = localized_fm_node->addChild("edgesWithGainChanges");
    edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() * sizeof(HyperedgeID));

    utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
    local_moves_node->updateSize(localMoves.capacity() * sizeof(std::pair<Move, MoveID>));

    fm_strategy.memoryConsumption(localized_fm_node);
    deltaPhg.memoryConsumption(localized_fm_node);
  }

  template<typename FMStrategy>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void AsyncKWayFM<FMStrategy>::releaseMoveLocksForLocalMovedNodes() {
    for (const auto& local_move : localMoves) {
      sharedData.nodeTracker.releaseNode(local_move.node, contraction_group_id);
    }
  }

    template<typename FMStrategy>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void AsyncKWayFM<FMStrategy>::lockUncontractionLockWithWaiting(const HypernodeID& hn) {
      while (!uncontraction_locks->tryToAcquireLock(hn, contraction_group_id)){/* keep trying to lock */}
    }

}   // namespace mt_kahypar


// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"

namespace mt_kahypar {
  template class AsyncKWayFM<GainCacheStrategy<AsyncFMSharedData>>;
  template class AsyncKWayFM<GainDeltaStrategy<AsyncFMSharedData>>;
  template class AsyncKWayFM<RecomputeGainStrategy<AsyncFMSharedData>>;
  template class AsyncKWayFM<GainCacheOnDemandStrategy<AsyncFMSharedData>>;
}