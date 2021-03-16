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

#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"

namespace mt_kahypar {

  template<typename FMStrategy>
  bool LocalizedKWayFM<FMStrategy>::findMoves(PartitionedHypergraph& phg, size_t taskID, size_t numSeeds) {
    localMoves.clear();
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    HypernodeID seedNode;
    while (runStats.pushes < numSeeds && sharedData.refinementNodes.try_pop(seedNode, taskID)) {
      SearchID previousSearchOfSeedNode = sharedData.nodeTracker.searchOfNode[seedNode].load(std::memory_order_relaxed);
      if (sharedData.nodeTracker.tryAcquireNode(seedNode, thisSearch)) {
        fm_strategy.insertIntoPQ(phg, seedNode, previousSearchOfSeedNode);
      }
    }

    if (runStats.pushes > 0) {
      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
      }

      if (context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints) {
        internalFindMoves<false>(phg);
      } else {
        deltaPhg.clear();
        deltaPhg.setPartitionedHypergraph(&phg);
        internalFindMoves<true>(phg);
        if (deltaPhg.combinedMemoryConsumption() > sharedData.deltaMemoryLimitPerThread) {
          sharedData.deltaExceededMemoryConstraints = true;
        }
      }
      return true;
    } else {
      return false;
    }
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
  void LocalizedKWayFM<FMStrategy>::acquireOrUpdateNeighbors(PHG& phg, const Move& move) {
    // Note: In theory we should acquire/update all neighbors. It just turned out that this works fine
    // Actually: only vertices incident to edges with gain changes can become new boundary vertices.
    // Vertices that already were boundary vertices, can still be considered later since they are in the task queue
    // --> actually not that bad
    for (HyperedgeID e : edgesWithGainChanges) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (neighborDeduplicator[v] != deduplicationTime) {
            SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(std::memory_order_acq_rel);
            if (searchOfV == thisSearch) {
              fm_strategy.updateGain(phg, v, move);
            } else if (sharedData.nodeTracker.tryAcquireNode(v, thisSearch)) {
              fm_strategy.insertIntoPQ(phg, v, searchOfV);
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


  template<typename FMStrategy>
  template<bool use_delta>
  void LocalizedKWayFM<FMStrategy>::internalFindMoves(PartitionedHypergraph& phg) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

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

      if constexpr (use_delta) {
        fm_strategy.deltaGainUpdates(deltaPhg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
      } else {
        fm_strategy.deltaGainUpdates(phg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
      }

    };

    // we can almost make this function take a generic partitioned hypergraph
    // we would have to add the success func to the interface of DeltaPhg (and then ignore it there...)
    // and do the local rollback outside this function


    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit) {

      if constexpr (use_delta) {
        if (!fm_strategy.findNextMove(deltaPhg, move)) break;
      } else {
        if (!fm_strategy.findNextMove(phg, move)) break;
      }

      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      bool moved = false;
      if (move.to != kInvalidPartition) {
        if constexpr (use_delta) {
          heaviestPartWeight = heaviestPartAndWeight(deltaPhg).second;
          fromWeight = deltaPhg.partWeight(move.from);
          toWeight = deltaPhg.partWeight(move.to);
          moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
                                          context.partition.max_part_weights[move.to], delta_func);
        } else {
          heaviestPartWeight = heaviestPartAndWeight(phg).second;
          fromWeight = phg.partWeight(move.from);
          toWeight = phg.partWeight(move.to);
          moved = phg.changeNodePart(move.node, move.from, move.to,
                                     context.partition.max_part_weights[move.to],
                                     [&] { move_id = sharedData.moveTracker.insertMove(move); }, delta_func);
        }
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += move.gain;
        localMoves.emplace_back(move, move_id);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();

          if constexpr (use_delta) {
            applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, true /* apply all moves */);
            bestImprovementIndex = 0;
            localMoves.clear();
            deltaPhg.clear();   // clear hashtables, save memory :)
          }
        }

        if constexpr (use_delta) {
          acquireOrUpdateNeighbors(deltaPhg, move);
        } else {
          acquireOrUpdateNeighbors(phg, move);
        }
      }

    }

    if constexpr (use_delta) {
      std::tie(bestImprovement, bestImprovementIndex) =
              applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, false);
    } else {
      revertToBestLocalPrefix(phg, bestImprovementIndex);
    }

    runStats.estimated_improvement = bestImprovement;
    fm_strategy.clearPQs(bestImprovementIndex);
    runStats.merge(stats);
  }


  template<typename FMStrategy>
  std::pair<Gain, size_t> LocalizedKWayFM<FMStrategy>::applyBestLocalPrefixToSharedPartition(
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
      assert(i < localMoves.size());
      Move& local_move = localMoves[i].first;
      MoveID& move_id = localMoves[i].second;
      attributed_gain = 0;

      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(local_move.node, local_move.from, local_move.to,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [&] { move_id = sharedData.moveTracker.insertMove(local_move); },
                                              delta_gain_func);
      } else {
        phg.changeNodePart(local_move.node, local_move.from, local_move.to,
                           std::numeric_limits<HypernodeWeight>::max(),
                           [&] { move_id = sharedData.moveTracker.insertMove(local_move); },
                           delta_gain_func);
      }

      attributed_gain = -attributed_gain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      improvement_from_attributed_gains += attributed_gain;
      ASSERT(move_id != std::numeric_limits<MoveID>::max());
      if (improvement_from_attributed_gains >= best_improvement_from_attributed_gains) {
        best_improvement_from_attributed_gains = improvement_from_attributed_gains;
        best_index_from_attributed_gains = i;
      }
    }

    runStats.local_reverts += localMoves.size() - best_index_locally_observed;
    if (!apply_all_moves && best_index_from_attributed_gains != best_index_locally_observed) {
      runStats.best_prefix_mismatch++;
    }

    // kind of double rollback, if attributed gains say we overall made things worse
    if (!apply_all_moves && improvement_from_attributed_gains < 0) {
      // always using the if-branch gave similar results
      runStats.local_reverts += best_index_locally_observed - best_index_from_attributed_gains + 1;
      for (size_t i = best_index_from_attributed_gains + 1; i < best_index_locally_observed; ++i) {
        Move& m = sharedData.moveTracker.getMove(localMoves[i].second);

        if constexpr (FMStrategy::uses_gain_cache) {
          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
        } else {
          phg.changeNodePart(m.node, m.to, m.from);
        }

        m.invalidate();
      }
      return std::make_pair(best_improvement_from_attributed_gains, best_index_from_attributed_gains);
    } else {
      return std::make_pair(best_improvement_locally_observed, best_index_locally_observed);
    }
  }

  template<typename FMStrategy>
  void LocalizedKWayFM<FMStrategy>::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {
    runStats.local_reverts += localMoves.size() - bestGainIndex;
    while (localMoves.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localMoves.back().second);
      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      } else {
        phg.changeNodePart(m.node, m.to, m.from);
      }
      m.invalidate();
      localMoves.pop_back();
    }
  }

  template<typename FMStrategy>
  void LocalizedKWayFM<FMStrategy>::memoryConsumption(utils::MemoryTreeNode *parent) const {
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

}   // namespace mt_kahypar


// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include <mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h>

namespace mt_kahypar {
  template class LocalizedKWayFM<GainCacheStrategy>;
  template class LocalizedKWayFM<GainDeltaStrategy>;
  template class LocalizedKWayFM<RecomputeGainStrategy>;
  template class LocalizedKWayFM<GainCacheOnDemandStrategy>;
}