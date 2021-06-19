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
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    for (const auto& seedNode : refinement_nodes ) {
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
      return _km1_delta;
    } else {
      return _km1_delta;
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
  void AsyncKWayFM<FMStrategy>::internalFindMoves(PartitionedHypergraph& phg) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID edge_size,
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
        _km1_delta += km1Delta(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
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

    while (!stopRule.searchShouldStop()) {

      auto reinsert_nodes = [&](const std::vector<HypernodeID>& nodes) {
          for (const auto& node : nodes) {
            SearchID previousSearchOfNode = sharedData.nodeTracker.searchOfNode[node].load(std::memory_order_relaxed);
            fm_strategy.insertIntoPQ(phg, node, previousSearchOfNode);
          }
      };

      // Find the next best move and attempt to lock the moved node until it works for a move. Then reinsert all nodes
      // for which locking failed.
      bool move_found = true;
      bool acquired_lock = false;
      std::vector<HypernodeID> nodes_to_reinsert;
      while (!acquired_lock) {
        if constexpr (use_delta) {
          move_found = fm_strategy.findNextMove(deltaPhg, move);
        } else {
          move_found = fm_strategy.findNextMove(phg, move);
        }

        if (!move_found) {
          if (nodes_to_reinsert.empty()) {
            // This means actually no moves are left
            break;
          } else {
            // This means for all moves the attempt to lock the node failed. Reinsert all and try again.
            // This is equivalent to waiting for a node from any move to become lockable.
            reinsert_nodes(nodes_to_reinsert);
            nodes_to_reinsert.clear();
            continue;
          }
        }

        ASSERT(move_found);
        acquired_lock = lock_manager->tryToAcquireLock(move.node, contraction_group_id);
        if (!acquired_lock) {
          nodes_to_reinsert.push_back(move.node);
        }
      }

      reinsert_nodes(nodes_to_reinsert);

      if (!move_found) {
        // Stop the search when no moves are left at all
        ASSERT(nodes_to_reinsert.empty());
        break;
      }

      ASSERT(move.isValid());
      ASSERT(lock_manager->isHeldBy(move.node, contraction_group_id));

      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
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
                                     [&] {}, delta_func);
        }
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
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();

          if constexpr (use_delta) {
            applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, true /* apply all moves */);
            bestImprovementIndex = 0;
            releaseLocksForLocalMovedNodes();
            localMoves.clear();
            deltaPhg.clear();   // clear hashtables, save memory :)
          }
        }

        if constexpr (use_delta) {
          acquireOrUpdateNeighbors(deltaPhg, move);
        } else {
          acquireOrUpdateNeighbors(phg, move);
        }
      } else {
        // Release lock on move node if the move does not go through
        lock_manager->strongReleaseLock(move.node, contraction_group_id);
      }
    }

    HEAVY_REFINEMENT_ASSERT(
        std::all_of(localMoves.begin(), localMoves.end(),
                    [&](const auto& local_move) {
                        return lock_manager->isHeldBy(local_move.node, contraction_group_id);
                    }));

    if constexpr (use_delta) {
      std::tie(bestImprovement, bestImprovementIndex) =
              applyBestLocalPrefixToSharedPartition(phg, bestImprovementIndex, bestImprovement, false);
    } else {
      revertToBestLocalPrefix(phg, bestImprovementIndex);
    }

    releaseLocksForLocalMovedNodes();

    runStats.estimated_improvement = bestImprovement;
    fm_strategy.clearPQs(bestImprovementIndex);
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
      assert(i < localMoves.size());
      Move& local_move = localMoves[i];
      attributed_gain = 0;

      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(local_move.node, local_move.from, local_move.to,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [&] {},
                                              delta_gain_func);
      } else {
        phg.changeNodePart(local_move.node, local_move.from, local_move.to,
                           std::numeric_limits<HypernodeWeight>::max(),
                           [&] {},
                           delta_gain_func);
      }

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

    // kind of double rollback, if attributed gains say we overall made things worse
    if (!apply_all_moves && improvement_from_attributed_gains < 0) {
      // always using the if-branch gave similar results
      runStats.local_reverts += best_index_locally_observed - best_index_from_attributed_gains + 1;
      for (size_t i = best_index_from_attributed_gains + 1; i < best_index_locally_observed; ++i) {
//        Move& m = sharedData.moveTracker.getMove(localMoves[i].second);
        Move& m = localMoves[i];
        attributed_gain = 0;

        if constexpr (FMStrategy::uses_gain_cache) {
          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from,
                                                std::numeric_limits<HypernodeWeight>::max(),
                                                [&] {},
                                                delta_gain_func);
        } else {
          phg.changeNodePart(m.node, m.to, m.from,
                             std::numeric_limits<HypernodeWeight>::max(),
                             [&] {},
                             delta_gain_func);
        }

        _km1_delta += attributed_gain;

        m.invalidate();
      }
      return std::make_pair(best_improvement_from_attributed_gains, best_index_from_attributed_gains);
    } else {
      return std::make_pair(best_improvement_locally_observed, best_index_locally_observed);
    }
  }

  template<typename FMStrategy>
  void AsyncKWayFM<FMStrategy>::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
        _km1_delta += km1Delta(he, edge_weight, edge_size,
                                    pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    runStats.local_reverts += localMoves.size() - bestGainIndex;
    while (localMoves.size() > bestGainIndex) {
//      Move& m = sharedData.moveTracker.getMove(localMoves.back().second);
      Move& m = localMoves.back();
      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [&] {},
                                              delta_gain_func);
      } else {
        phg.changeNodePart(m.node, m.to, m.from,
                           std::numeric_limits<HypernodeWeight>::max(),
                           [&] {},
                           delta_gain_func);
      }
      m.invalidate();
      localMoves.pop_back();
    }
  }

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
  void AsyncKWayFM<FMStrategy>::releaseLocksForLocalMovedNodes() {
    for (const auto& local_move : localMoves) {
      lock_manager->strongReleaseLock(local_move.node, contraction_group_id);
    }
  }

}   // namespace mt_kahypar


// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include <mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h>

namespace mt_kahypar {
  template class AsyncKWayFM<GainCacheStrategy>;
  template class AsyncKWayFM<GainDeltaStrategy>;
  template class AsyncKWayFM<RecomputeGainStrategy>;
  template class AsyncKWayFM<GainCacheOnDemandStrategy>;
}