/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
            SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(std::memory_order_relaxed);
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

      // skip if no target block available
      if (move.to == kInvalidPartition) {
        continue;
      }

      bool expect_improvement = estimatedImprovement + move.gain > bestImprovement;
      bool high_deg = phg.nodeDegree(move.node) >= PartitionedHypergraph::HIGH_DEGREE_THRESHOLD;

      // skip if high degree (unless it nets actual improvement; but don't apply on deltaPhg then)
      if (!expect_improvement && high_deg) {
        continue;
      }
      // less restrictive option: skip if negative gain (or < -5000 or smth).
      // downside: have to flush before improvement or run it through deltaPhg
      // probably quite similar since this only really matters in the first few moves where the stop rule
      // doesn't signal us to stop yet

      edgesWithGainChanges.clear(); // clear before move. delta_func feeds nets of moved vertex.
      MoveID move_id = std::numeric_limits<MoveID>::max();
      bool moved = false;
      if constexpr (use_delta) {
        heaviestPartWeight = heaviestPartAndWeight(deltaPhg).second;
        fromWeight = deltaPhg.partWeight(move.from);
        toWeight = deltaPhg.partWeight(move.to);
        if (expect_improvement) {
          // since we will flush the move sequence, don't bother running it through the deltaPhg
          // this is intended to allow moving high deg nodes (blow up hash tables) if they give an improvement.
          // The nets affected by a gain cache update are collected when we apply this improvement on the
          // global partition (used to expand the localized search and update the gain values).
          moved = toWeight + phg.nodeWeight(move.node) <= context.partition.max_part_weights[move.to];
        } else {
          moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
                                          context.partition.max_part_weights[move.to], delta_func);
        }
      } else {
        heaviestPartWeight = heaviestPartAndWeight(phg).second;
        fromWeight = phg.partWeight(move.from);
        toWeight = phg.partWeight(move.to);
        moved = phg.changeNodePart(move.node, move.from, move.to,
                                   context.partition.max_part_weights[move.to],
                                   [&] { move_id = sharedData.moveTracker.insertMove(move); }, delta_func);
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += move.gain;
        localMoves.emplace_back(move, move_id);
        stopRule.update(move.gain);
        bool improved_km1 = estimatedImprovement > bestImprovement;
        bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
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

        // no need to update our PQs if we stop anyways
        if (stopRule.searchShouldStop()
              || sharedData.finishedTasks.load(std::memory_order_relaxed) >= sharedData.finishedTasksLimit) {
          break;
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
          const bool apply_delta_improvement) {

    Gain improvement_from_attributed_gains = 0;
    Gain attributed_gain = 0;
    bool is_last_move = false;

    auto delta_gain_func = [&](const HyperedgeID he,
                                    const HyperedgeWeight edge_weight,
                                    const HypernodeID edge_size,
                                    const HypernodeID pin_count_in_from_part_after,
                                    const HypernodeID pin_count_in_to_part_after) {
      attributed_gain += km1Delta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);

      // Gains of the pins of a hyperedge can only change in the following situations.
      if ( is_last_move &&
           ( pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
             pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2 ) ) {
        // This vector is used by the acquireOrUpdateNeighbor function to expand to neighbors
        // or update the gain values of neighbors of the moved node and is cleared afterwards.
        // BEWARE. Adding the nets at this stage works, because the vector is cleared before the move,
        // and the expansion happens after applyBestLocalPrefixToSharedPartition.
        edgesWithGainChanges.push_back(he);
      }

      // TODO: We have different strategies to maintain the gain values during an FM search.
      // Some use the gain cache, others compute them each time from scratch or use delta gain updates.
      // In case the delta gain update strategy is used, we would have to call the deltaGainUpdate function
      // of the FM strategy here. However, the current strategy in our presets use the gain cache and calling the deltaGainUpdate
      // function would apply the updates on the thread-local partition, which we do not want here.
      // Keep in mind that the gain values of the FMGainDeltaStrategy might be incorrect afterwards.
      // However, this strategy is only experimental We should remove the
      // different strategies since we do not use them.
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain best_improvement_from_attributed_gains = 0;
    size_t best_index_from_attributed_gains = 0;
    for (size_t i = 0; i < best_index_locally_observed; ++i) {
      ASSERT(i < localMoves.size());
      Move& local_move = localMoves[i].first;
      MoveID& move_id = localMoves[i].second;
      attributed_gain = 0;
      // In a localized FM search, we apply all moves to a thread-local partition (delta_phg)
      // using hash tables. Once we find an improvement, we apply the corresponding move
      // sequence to the global partition. To save memory (in the hash tables), we do not apply
      // the last move that leads to the improvement to the thread-local partition as we reset them after
      // an improvement is found, anyways. However, when applying a move on the thread-local partition,
      // we collect all nets affected by a gain cache update and expand the search to pins
      // contained in these nets. Since, we do not apply last move on the thread-local partition we collect
      // these nets here.
      is_last_move = apply_delta_improvement && i == best_index_locally_observed - 1;

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
    if (!apply_delta_improvement && best_index_from_attributed_gains != best_index_locally_observed) {
      runStats.best_prefix_mismatch++;
    }

    // kind of double rollback, if attributed gains say we overall made things worse
    if (!apply_delta_improvement && improvement_from_attributed_gains < 0) {
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
