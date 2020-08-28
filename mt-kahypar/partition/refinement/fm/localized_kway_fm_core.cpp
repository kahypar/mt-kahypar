#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"

namespace mt_kahypar {

  bool LocalizedKWayFM::findMovesLocalized(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t taskID) {
    localData.clear();
    validHyperedges.clear();

    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    const size_t nSeeds = context.refinement.fm.num_seed_nodes;
    HypernodeID seedNode;
    while (localData.runStats.pushes < nSeeds && sharedData.refinementNodes.try_pop(seedNode, taskID)) {
      if (insertPQ(phg, seedNode, sharedData)) {
        localData.seedVertices.push_back(seedNode);
      }
    }
    for (PartitionID i = 0; i < k; ++i) {
      updateBlock(i);
    }

    if (localData.runStats.pushes > 0) {
      if (!context.refinement.fm.perform_moves_global
          && deltaPhg.combinedMemoryConsumption() > sharedData.deltaMemoryLimitPerThread) {
        sharedData.deltaExceededMemoryConstraints = true;
      }

      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
      }

      if (context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints) {
        internalFindMovesOnGlobalHypergraph(phg, sharedData);
      } else {
        deltaPhg.clear();
        deltaPhg.setPartitionedHypergraph(&phg);
        internalFindMovesOnDeltaHypergraph(phg, sharedData);
      }
      return true;
    } else {
      return false;
    }

  }

  bool LocalizedKWayFM::findMovesUsingFullBoundary(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    localData.clear();
    validHyperedges.clear();
    thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;

    for (HypernodeID u : sharedData.refinementNodes.safely_inserted_range()) {
      insertPQ(phg, u, sharedData);
    }
    for (PartitionID i = 0; i < k; ++i) {
      updateBlock(i);
    }

    // this is boundary FM, so it's sequential. no need for delta hypergraph
    internalFindMovesOnGlobalHypergraph(phg, sharedData);
    return true;
  }


  template<typename Partition>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight> heaviestPartAndWeight(const Partition& partition) {
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


  void LocalizedKWayFM:: internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg,
                                                            FMSharedData& sharedData) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    auto hes_to_update_func = [&](const HyperedgeID he,
                                  const HyperedgeWeight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situation.
      // In such cases, we mark the hyperedge as invalid and update the gain of all
      // pins afterwards.
      if ( pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
           pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2 ) {
        validHyperedges[he] = false;
      }
    };

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
           && findNextMove(deltaPhg, move)) {
      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);

      bool moved = false;
      if (move.to != kInvalidPartition) {
        heaviestPartWeight = heaviestPartAndWeight(deltaPhg).second;
        fromWeight = deltaPhg.partWeight(move.from);
        toWeight = deltaPhg.partWeight(move.to);
        moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
                                        context.partition.max_part_weights[move.to], hes_to_update_func);
      }

      if (moved) {
        localData.runStats.moves++;
        estimatedImprovement += move.gain;
        localData.localMoves.push_back(move);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localData.localMoves.size();
        }

        insertOrUpdateNeighbors(deltaPhg, sharedData, move);
      }

      for (PartitionID i = 0; i < k; ++i) {
        updateBlock(i);
      }
    }

    std::tie(bestImprovement, bestImprovementIndex) =
            applyMovesOnGlobalHypergraph(phg, sharedData, bestImprovementIndex, bestImprovement);
    localData.runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    localData.runStats.merge(stats);
  }

  void LocalizedKWayFM::internalFindMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                                            FMSharedData& sharedData) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;


    auto hes_to_update_func = [&](const HyperedgeID he,
                                  const HyperedgeWeight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      // In such cases, we mark the hyperedge as invalid and update the gain of all
      // pins afterwards.
      if ( pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
           pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2 ) {
        validHyperedges[he] = false;
      }
    };

    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit
           && findNextMove(phg, move)) {
      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      auto report_success = [&] { move_id = sharedData.moveTracker.insertMove(move); };

      bool moved = false;
      if (move.to != kInvalidPartition) {
        heaviestPartWeight = heaviestPartAndWeight(phg).second;
        fromWeight = phg.partWeight(move.from);
        toWeight = phg.partWeight(move.to);
        moved = phg.changeNodePartFullUpdate(move.node, move.from, move.to,
                                             context.partition.max_part_weights[move.to],
                                             report_success, hes_to_update_func);
      }

      if (moved) {
        localData.runStats.moves++;
        estimatedImprovement += move.gain;
        localData.localMoveIDs.push_back(move_id);
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localData.localMoveIDs.size();
        }

        insertOrUpdateNeighbors(phg, sharedData, move);
      }

      for (PartitionID i = 0; i < k; ++i) {
        updateBlock(i);
      }
    }

    revertToBestLocalPrefix(phg, sharedData, bestImprovementIndex);
    localData.runStats.estimated_improvement = bestImprovement;
    clearPQs(sharedData, bestImprovementIndex);
    localData.runStats.merge(stats);
  }


  std::pair<Gain, size_t> LocalizedKWayFM::applyMovesOnGlobalHypergraph(PartitionedHypergraph& phg,
                                                       FMSharedData& sharedData,
                                                       const size_t bestGainIndex,
                                                       const Gain bestEstimatedImprovement) {
    ASSERT(localData.localMoveIDs.empty());
    Gain estimatedImprovement = 0;
    Gain lastGain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      lastGain += km1Delta(he, edge_weight, edge_size,
                           pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain bestImprovement = 0;
    size_t bestIndex = 0;
    for ( size_t i = 0; i < bestGainIndex; ++i ) {
      Move& m = localData.localMoves[i];
      MoveID m_id = std::numeric_limits<MoveID>::max();
      lastGain = 0;
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(),
                                   [&] { m_id = sharedData.moveTracker.insertMove(m); }, delta_gain_func);
      lastGain = -lastGain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      estimatedImprovement += lastGain;
      ASSERT(m_id != std::numeric_limits<MoveID>::max());
      Move& move = sharedData.moveTracker.getMove(m_id);
      move.gain = lastGain; // Update gain value based on hypergraph delta
      localData.localMoveIDs.push_back(m_id);
      if ( estimatedImprovement >= bestImprovement ) {  // TODO also incorporate balance into this
        bestImprovement = estimatedImprovement;
        bestIndex = i;
      }
    }

    // Kind of double rollback, if gain values are not correct
    ASSERT(localData.localMoveIDs.size() == bestGainIndex);
    if ( estimatedImprovement < 0 ) {
      localData.runStats.local_reverts += localData.localMoves.size() - bestIndex;
      for ( size_t i = bestIndex + 1; i < bestGainIndex; ++i ) {
        Move& m = sharedData.moveTracker.getMove(localData.localMoveIDs[i]);
        phg.changeNodePartFullUpdate(m.node, m.to, m.from);
        sharedData.moveTracker.invalidateMove(m);
      }
      localData.runStats.local_reverts += localData.localMoves.size() - bestGainIndex;
      return std::make_pair(bestImprovement, bestIndex);
    } else {
      return std::make_pair(bestEstimatedImprovement, bestGainIndex);
    }
  }

  void LocalizedKWayFM::revertToBestLocalPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData, size_t bestGainIndex) {
    localData.runStats.local_reverts += localData.localMoveIDs.size() - bestGainIndex;
    while (localData.localMoveIDs.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localData.localMoveIDs.back());
      phg.changeNodePartFullUpdate(m.node, m.to, m.from);
      sharedData.moveTracker.invalidateMove(m);
      localData.localMoveIDs.pop_back();
    }
  }

  void LocalizedKWayFM::clearPQs(FMSharedData& sharedData, const size_t bestImprovementIndex) {
    // release all nodes that were not moved
    // reinsert into task queue only if we're doing multitry and at least one node was moved
    // unless a node was moved, only seed nodes are in the pqs
    const bool release = context.refinement.fm.release_nodes
                         && context.refinement.fm.algorithm == FMAlgorithm::fm_multitry
                         && localData.runStats.moves > 0;
    const bool reinsert_seeds = bestImprovementIndex > 0;

    if (release) {
      if (!reinsert_seeds) {
        for (HypernodeID u : localData.seedVertices) {
          sharedData.fruitlessSeed.set(u, true);
        }
      }

      auto release_node = [&](const HypernodeID& v) {
        sharedData.nodeTracker.releaseNode(v);
        if (!sharedData.fruitlessSeed[v] && sharedData.refinementNodes.was_pushed_and_removed(v)) {
          sharedData.refinementNodes.concurrent_push(v);
          localData.runStats.task_queue_reinsertions++;
        }
      };

      // Release all nodes contained in PQ
      for (PartitionID i = 0; i < k; ++i) {
        for (PosT j = 0; j < vertexPQs[i].size(); ++j) {
          const HypernodeID v = vertexPQs[i].at(j);
          release_node(v);
        }
      }
    }

    for (PartitionID i = 0; i < k; ++i) {
      vertexPQs[i].clear();
    }
    blockPQ.clear();
  }

  void LocalizedKWayFM::memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode* deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(updateDeduplicator.size_in_bytes());
    utils::MemoryTreeNode* valid_hyperedges_node = localized_fm_node->addChild("Valid Hyperedges");
    valid_hyperedges_node->updateSize(validHyperedges.size_in_bytes());
    utils::MemoryTreeNode* block_pq_node = localized_fm_node->addChild("Block PQ");
    block_pq_node->updateSize(blockPQ.size_in_bytes());
    utils::MemoryTreeNode* vertex_pq_node = localized_fm_node->addChild("Vertex PQ");
    for ( const VertexPriorityQueue& pq : vertexPQs ) {
      vertex_pq_node->updateSize(pq.size_in_bytes());
    }
    localData.memoryConsumption(localized_fm_node);
    deltaPhg.memoryConsumption(localized_fm_node);
  }

}