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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/unconstrained_strategy.h"

namespace mt_kahypar {

  template<typename TypeTraits, typename GainTypes>
  template<typename DispatchedFMStrategy>
  bool LocalizedKWayFM<TypeTraits, GainTypes>::findMoves(DispatchedFMStrategy& fm_strategy, PartitionedHypergraph& phg,
                                                         size_t taskID, size_t numSeeds) {
    localMoves.clear();
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    HypernodeID seedNode;
    while (runStats.pushes < numSeeds && sharedData.refinementNodes.try_pop(seedNode, taskID)) {
      if (sharedData.nodeTracker.tryAcquireNode(seedNode, thisSearch)) {
        fm_strategy.insertIntoPQ(phg, gain_cache, seedNode);
      }
    }

    if (runStats.pushes > 0) {
      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
        delta_gain_cache.dropMemory();
      }

      if (context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints) {
        if ( phg.hasFixedVertices() ) {
          internalFindMoves<false, true>(phg, fm_strategy);
        } else {
          internalFindMoves<false, false>(phg, fm_strategy);
        }
      } else {
        deltaPhg.clear();
        delta_gain_cache.clear();
        deltaPhg.setPartitionedHypergraph(&phg);
        if ( phg.hasFixedVertices() ) {
          internalFindMoves<true, true>(phg, fm_strategy);
        } else {
          internalFindMoves<true, false>(phg, fm_strategy);
        }
      }
      return true;
    } else {
      return false;
    }
  }

  template<typename Partition>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight>
  heaviestPartAndWeight(const Partition& partition, const PartitionID k) {
    PartitionID p = kInvalidPartition;
    HypernodeWeight w = std::numeric_limits<HypernodeWeight>::min();
    for (PartitionID i = 0; i < k; ++i) {
      if (partition.partWeight(i) > w) {
        w = partition.partWeight(i);
        p = i;
      }
    }
    return std::make_pair(p, w);
  }

  template<typename TypeTraits, typename GainTypes>
  template<bool has_fixed_vertices, typename PHG, typename CACHE, typename DispatchedFMStrategy>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void LocalizedKWayFM<TypeTraits, GainTypes>::acquireOrUpdateNeighbors(PHG& phg, CACHE& gain_cache, const Move& move,
                                                                        DispatchedFMStrategy& fm_strategy) {
    // Note: In theory we should acquire/update all neighbors. It just turned out that this works fine
    // Actually: only vertices incident to edges with gain changes can become new boundary vertices.
    // Vertices that already were boundary vertices, can still be considered later since they are in the task queue
    // --> actually not that bad
    for (HyperedgeID e : edgesWithGainChanges) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if constexpr ( has_fixed_vertices ) {
            if ( phg.isFixed(v) ) continue;
          }
          if (neighborDeduplicator[v] != deduplicationTime) {
            SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(std::memory_order_relaxed);
            if (searchOfV == thisSearch) {
              fm_strategy.updateGain(phg, gain_cache, v, move);
            } else if (sharedData.nodeTracker.tryAcquireNode(v, thisSearch)) {
              fm_strategy.insertIntoPQ(phg, gain_cache, v);
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


  template<typename TypeTraits, typename GainTypes>
  template<bool use_delta, bool has_fixed_vertices, typename DispatchedFMStrategy>
  void LocalizedKWayFM<TypeTraits, GainTypes>::internalFindMoves(PartitionedHypergraph& phg,
                                                                 DispatchedFMStrategy& fm_strategy) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    auto delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if ( GainCache::triggersDeltaGainUpdate(sync_update) ) {
        edgesWithGainChanges.push_back(sync_update.he);
      }

      if constexpr (use_delta) {
        fm_strategy.deltaGainUpdates(deltaPhg, delta_gain_cache, sync_update);
      } else {
        fm_strategy.deltaGainUpdates(phg, gain_cache, sync_update);
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
        if (!fm_strategy.findNextMove(deltaPhg, delta_gain_cache, move)) break;
      } else {
        if (!fm_strategy.findNextMove(phg, gain_cache, move)) break;
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
      const HypernodeWeight allowed_weight = DispatchedFMStrategy::is_unconstrained ? std::numeric_limits<HypernodeWeight>::max()
                                             : context.partition.max_part_weights[move.to];
      if constexpr (use_delta) {
        heaviestPartWeight = heaviestPartAndWeight(deltaPhg, context.partition.k).second;
        fromWeight = deltaPhg.partWeight(move.from);
        toWeight = deltaPhg.partWeight(move.to);
        if (expect_improvement) {
          // since we will flush the move sequence, don't bother running it through the deltaPhg
          // this is intended to allow moving high deg nodes (blow up hash tables) if they give an improvement.
          // The nets affected by a gain cache update are collected when we apply this improvement on the
          // global partition (used to expand the localized search and update the gain values).
          moved = toWeight + phg.nodeWeight(move.node) <= allowed_weight;
        } else {
          moved = deltaPhg.changeNodePart(move.node, move.from, move.to, allowed_weight, delta_func);
          fm_strategy.applyMove(deltaPhg, delta_gain_cache, move, false);
        }
      } else {
        heaviestPartWeight = heaviestPartAndWeight(phg, context.partition.k).second;
        fromWeight = phg.partWeight(move.from);
        toWeight = phg.partWeight(move.to);
        moved = phg.changeNodePart(move.node, move.from, move.to, allowed_weight,
                                   [&] { move_id = sharedData.moveTracker.insertMove(move); }, delta_func);
        fm_strategy.applyMove(phg, gain_cache, move, true);
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
            applyBestLocalPrefixToSharedPartition(phg, fm_strategy, bestImprovementIndex);
            bestImprovementIndex = 0;
            localMoves.clear();
            deltaPhg.clear();   // clear hashtables, save memory :)
            delta_gain_cache.clear();
          }
        }

        // no need to update our PQs if we stop anyways
        if (stopRule.searchShouldStop()
              || sharedData.finishedTasks.load(std::memory_order_relaxed) >= sharedData.finishedTasksLimit) {
          break;
        }

        if constexpr (use_delta) {
          acquireOrUpdateNeighbors<has_fixed_vertices>(deltaPhg, delta_gain_cache, move, fm_strategy);
        } else {
          acquireOrUpdateNeighbors<has_fixed_vertices>(phg, gain_cache, move, fm_strategy);
        }
      }

    }

    if constexpr (use_delta) {
      // in this case there is no improved local prefix to apply (was already applied in the loop)
      ASSERT(bestImprovementIndex == 0);
    } else {
      revertToBestLocalPrefix(phg, fm_strategy, bestImprovementIndex);
    }

    fm_strategy.reset();
    runStats.estimated_improvement = bestImprovement;
    runStats.merge(stats);
  }


  template<typename TypeTraits, typename GainTypes>
  template<typename DispatchedFMStrategy>
  void LocalizedKWayFM<TypeTraits, GainTypes>::applyBestLocalPrefixToSharedPartition(
          PartitionedHypergraph& phg,
          DispatchedFMStrategy& fm_strategy,
          const size_t best_index_locally_observed) {
    // Note: if this precondition does not hold, the call to fm_strategy.flushLocalChanges() would be incorrect
    ASSERT(best_index_locally_observed == localMoves.size());

    bool is_last_move = false;

    auto delta_gain_func = [&](const SynchronizedEdgeUpdate& sync_update) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if ( is_last_move && GainCache::triggersDeltaGainUpdate(sync_update) ) {
        // This vector is used by the acquireOrUpdateNeighbor function to expand to neighbors
        // or update the gain values of neighbors of the moved node and is cleared afterwards.
        // BEWARE. Adding the nets at this stage works, because the vector is cleared before the move,
        // and the expansion happens after applyBestLocalPrefixToSharedPartition.
        edgesWithGainChanges.push_back(sync_update.he);
      }
    };

    // Apply move sequence to original hypergraph
    for (size_t i = 0; i < best_index_locally_observed; ++i) {
      ASSERT(i < localMoves.size());
      Move& local_move = localMoves[i].first;
      MoveID& move_id = localMoves[i].second;
      // In a localized FM search, we apply all moves to a thread-local partition (delta_phg)
      // using hash tables. Once we find an improvement, we apply the corresponding move
      // sequence to the global partition. To save memory (in the hash tables), we do not apply
      // the last move that leads to the improvement to the thread-local partition as we reset them after
      // an improvement is found, anyways. However, when applying a move on the thread-local partition,
      // we collect all nets affected by a gain cache update and expand the search to pins
      // contained in these nets. Since, we do not apply last move on the thread-local partition we collect
      // these nets here.
      is_last_move = (i == best_index_locally_observed - 1);

      phg.changeNodePart(
        gain_cache, local_move.node, local_move.from, local_move.to,
        std::numeric_limits<HypernodeWeight>::max(),
        [&] { move_id = sharedData.moveTracker.insertMove(local_move); },
        delta_gain_func);
      ASSERT(move_id != std::numeric_limits<MoveID>::max());
    }
    fm_strategy.flushLocalChanges();
  }

  template<typename TypeTraits, typename GainTypes>
  template<typename DispatchedFMStrategy>
  void LocalizedKWayFM<TypeTraits, GainTypes>::revertToBestLocalPrefix(PartitionedHypergraph& phg,
                                                                       DispatchedFMStrategy& fm_strategy,
                                                                       size_t bestGainIndex) {
    runStats.local_reverts += localMoves.size() - bestGainIndex;
    while (localMoves.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localMoves.back().second);
      phg.changeNodePart(gain_cache, m.node, m.to, m.from);
      fm_strategy.revertMove(phg, gain_cache, m, true);
      m.invalidate();
      localMoves.pop_back();
    }
  }

  template<typename TypeTraits, typename GainTypes>
  void LocalizedKWayFM<TypeTraits, GainTypes>::changeNumberOfBlocks(const PartitionID new_k) {
    deltaPhg.changeNumberOfBlocks(new_k);
    blockPQ.resize(new_k);
    for ( VertexPriorityQueue& pq : vertexPQs ) {
      pq.setHandle(sharedData.vertexPQHandles.data(), sharedData.numberOfNodes);
    }
    while ( static_cast<size_t>(new_k) > vertexPQs.size() ) {
      vertexPQs.emplace_back(sharedData.vertexPQHandles.data(), sharedData.numberOfNodes);
    }
  }

  template<typename TypeTraits, typename GainTypes>
  void LocalizedKWayFM<TypeTraits, GainTypes>::memoryConsumption(utils::MemoryTreeNode *parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode *localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode *deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(neighborDeduplicator.capacity() * sizeof(HypernodeID));
    utils::MemoryTreeNode *edges_to_activate_node = localized_fm_node->addChild("edgesWithGainChanges");
    edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() * sizeof(HyperedgeID));

    size_t vertex_pq_sizes = std::accumulate(
            vertexPQs.begin(), vertexPQs.end(), 0,
            [](size_t init, const VertexPriorityQueue& pq) { return init + pq.size_in_bytes(); }
    );
    localized_fm_node->addChild("PQs", blockPQ.size_in_bytes() + vertex_pq_sizes);

    utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
    local_moves_node->updateSize(localMoves.capacity() * sizeof(std::pair<Move, MoveID>));

    deltaPhg.memoryConsumption(localized_fm_node);
    delta_gain_cache.memoryConsumption(localized_fm_node);
  }

  namespace {
  #define LOCALIZED_KWAY_FM(X, Y) LocalizedKWayFM<X, Y>;                                                      \
    template bool LocalizedKWayFM<X, Y>::findMoves(LocalUnconstrainedStrategy&,                               \
                    typename LocalizedKWayFM<X, Y>::PartitionedHypergraph&, size_t, size_t);                  \
    template bool LocalizedKWayFM<X, Y>::findMoves(LocalGainCacheStrategy&,                                   \
                    typename LocalizedKWayFM<X, Y>::PartitionedHypergraph&, size_t, size_t)
  }

  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(LOCALIZED_KWAY_FM)

}   // namespace mt_kahypar
