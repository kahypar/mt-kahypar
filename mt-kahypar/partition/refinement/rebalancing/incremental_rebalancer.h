//
// Created by tobias on 01.02.25.
//
#pragma once


#include <vector>
#include <mt-kahypar/datastructures/priority_queue.h>
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include <mt-kahypar/partition/refinement/gains/gain_cache_ptr.h>

namespace mt_kahypar {

namespace rebalancer {

} // namespace rebalancer

class IncrementalRebalancer final {
public:

    void init(ds::PartitionedHypergraph<ds::MutableHypergraph> &partitioned_hypergraph, Context &context, gain_cache_t &gain_cache, vec<Gain> &benefit_aggregator) {
      partitioned_hypergraph_m = &partitioned_hypergraph;
      _context = &context;
      _gain_cache = &gain_cache;
      _benefit_aggregator = &benefit_aggregator;
      populateBlockQueues();
    }

    void reset() {
      ASSERT(partitioned_hypergraph_m != nullptr);
      ASSERT(_context != nullptr);
      ASSERT(_gain_cache != nullptr);
      ASSERT(_benefit_aggregator != nullptr);
      if (partitioned_hypergraph_m == nullptr || _context == nullptr || _gain_cache == nullptr || _benefit_aggregator == nullptr) {
        std::cout << "IncrementalRebalancer not initialized properly!" << std::endl;
        exit(1);
      }
      // std::cout << "Resetting IncrementalRebalancer 1" << std::endl;
      if (_blocks.empty()) {
        std::cout << "IncrementalRebalancer blocks not initialized properly!" << std::endl;
        exit(1);
      }
      // std::cout << "Resetting IncrementalRebalancer 2" << std::endl;
      // _blocks.clear();
      // // std::cout << "Resetting IncrementalRebalancer" << std::endl;
      // populateBlockQueues();

      // build new blocks first, then swap
      std::vector<BlockQueues> new_blocks;
      new_blocks.reserve(_context->partition.k);
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        new_blocks.emplace_back(partitioned_hypergraph_m->initialNumNodes() * 2);
      }
      // populate new_blocks similarly to populateBlockQueues()...


      for (HypernodeID u = 0; u < partitioned_hypergraph_m->initialNumNodes(); ++u) {
        if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
          continue;
        }
        const PartitionID part = partitioned_hypergraph_m->partID(u);
        HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == part) {
            continue;
          }
          const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
          if (gain > highest_gain) {
            highest_gain = gain;
          }
          HyperedgeWeight weighted_pull_gain = (gain > 0) ? gain / partitioned_hypergraph_m->nodeWeight(u) : gain * partitioned_hypergraph_m->nodeWeight(u);
          new_blocks[b].pull.insert(u, weighted_pull_gain);
        }
        HyperedgeWeight weighted_push_gain = (highest_gain > 0) ? highest_gain * partitioned_hypergraph_m->nodeWeight(u) : highest_gain / partitioned_hypergraph_m->nodeWeight(u);
        new_blocks[part].push.insert(u, weighted_push_gain);
      }

      // finally swap
      _blocks.swap(new_blocks);
      // std::cout << "Resetting IncrementalRebalancer 3" << std::endl;
      ASSERT(checkBlockQueues());
      // std::cout << "Resetting IncrementalRebalancer 4" << std::endl;
    }

    // void updateAllForMove(Move move) {
    //   if (!partitioned_hypergraph_m->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(*_gain_cache)))
    //   {
    //     std::cout << "Gain cache is not valid" << std::endl;
    //     exit(1);
    //   }
    //   _context->dynamic.incremental_km1 -= GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(move.node, move.from, move.to);
    //   updateGainCacheForMove(move);
    //   updateHeapsForMove(move);
    // }

    std::tuple<HyperedgeWeight, std::vector<HypernodeID>> rebalanceAndUpdateGainCache() {
      HyperedgeWeight total_gain = 0;
      std::vector<HypernodeID> moved_nodes;
      while (!mt_kahypar::metrics::isBalanced(*partitioned_hypergraph_m, *_context)) {
        PartitionID imbalanced_block = kInvalidPartition;
        for (PartitionID i = 0; i < _context->partition.k; ++i) {
          if (partitioned_hypergraph_m->partWeight(i) > _context->partition.max_part_weights[i]) {
            imbalanced_block = i;
            break;
          }
        }
        ASSERT(imbalanced_block != kInvalidPartition);
        ASSERT(!_blocks[imbalanced_block].push.empty());
        std::vector<std::tuple<PartitionID, HyperedgeWeight>> violations;
        while (!_blocks[imbalanced_block].push.empty() && partitioned_hypergraph_m->partWeight(imbalanced_block) > _context->partition.max_part_weights[imbalanced_block]) {
          const HypernodeID u = _blocks[imbalanced_block].push.top();
          const HyperedgeWeight stored_gain = _blocks[imbalanced_block].push.topKey();
          _blocks[imbalanced_block].push.deleteTop();

          if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
            continue;
          }

          //TODO: Why is this happening with mixed queries?
          if (partitioned_hypergraph_m->partID(u) != imbalanced_block) {
            std::cout << "Node " << u << " is not in block " << imbalanced_block << " but in " << partitioned_hypergraph_m->partID(u) << std::endl;
            insertOrUpdateNode(u, partitioned_hypergraph_m->partID(u));
            continue;
          }

          ASSERT(partitioned_hypergraph_m->partID(u) == imbalanced_block);

          auto [best_push_block, gain] = findFittingBestPushBlock(u, imbalanced_block);
          ASSERT(best_push_block != kInvalidPartition && best_push_block != imbalanced_block);

          //reassure that the gain has not changed
          const HyperedgeWeight weighted_gain = (gain > 0) ? gain * partitioned_hypergraph_m->nodeWeight(u) : gain / partitioned_hypergraph_m->nodeWeight(u);
          if (weighted_gain != stored_gain) {
            _blocks[imbalanced_block].push.insert(u, weighted_gain);
            continue;
          }

          Move move = {imbalanced_block, best_push_block, u, gain};
          partitioned_hypergraph_m->changeNodePart(move.node, move.from, move.to);
          updateGainCacheForMove(move);
          updateHeapsForMove(move);

          moved_nodes.push_back(u);

          total_gain += gain;
        }
      }
      return std::make_tuple(total_gain, moved_nodes);
    }

    std::tuple<HyperedgeWeight, std::vector<HypernodeID>> pullAndUpdateGainCache(PartitionID block_pull_target) {
      HyperedgeWeight total_gain = 0;
      std::vector<HypernodeID> moved_nodes;
      std::vector<std::tuple<PartitionID, HyperedgeWeight>> violations;
      while (!_blocks[block_pull_target].pull.empty()) {
        const HypernodeID u = _blocks[block_pull_target].pull.top();
        const HyperedgeWeight stored_gain = _blocks[block_pull_target].pull.topKey();

        if (stored_gain <= 0) {
          break;
        }

        _blocks[block_pull_target].pull.deleteTop();

        if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
          continue;
        }
        PartitionID block_pull_source = partitioned_hypergraph_m->partID(u);

        //TODO: Why is this happening with mixed queries?
        //TODO: Because v_cycle changes are not reflected in the rebalancer
        if (block_pull_source == block_pull_target) {
          std::cout << "Node " << u << " is in block " << block_pull_source << " and should be moved to " << block_pull_target << std::endl;
          insertOrUpdateNode(u, block_pull_target);
          continue;
        }
        ASSERT(block_pull_source != block_pull_target);

        HyperedgeWeight gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, block_pull_source, block_pull_target);

        //reassure that the gain has not changed
        const HyperedgeWeight weighted_gain = (gain > 0) ? gain / partitioned_hypergraph_m->nodeWeight(u) : gain * partitioned_hypergraph_m->nodeWeight(u);
        if (weighted_gain != stored_gain) {
          _blocks[block_pull_target].pull.insert(u, weighted_gain);
          continue;
        }

        //prevent moves that violate max_part_weights
        if (partitioned_hypergraph_m->partWeight(block_pull_target) + partitioned_hypergraph_m->nodeWeight(u) >
            _context->partition.max_part_weights[block_pull_target]) {
          violations.push_back(std::make_tuple(u, weighted_gain));
          continue;
        }

        Move move = {block_pull_source, block_pull_target, u, gain};
        partitioned_hypergraph_m->changeNodePart(move.node, move.from, move.to);
        updateGainCacheForMove(move);
        updateHeapsForMove(move);

        moved_nodes.push_back(u);

        total_gain += gain;
      }

      //reinsert skipped nodes
      for (const auto& violation : violations) {
        _blocks[block_pull_target].pull.insert(std::get<0>(violation), std::get<1>(violation));
      }

      return std::make_tuple(total_gain, moved_nodes);
    }

    void updateHeapsForMove(Move move) {
      HypernodeID u = move.node;
      ASSERT(partitioned_hypergraph_m->nodeIsEnabled(u));

      //remove node from old push queue
      if (_blocks[move.from].push.contains(u)) {
        _blocks[move.from].push.remove(u);
      }

      //remove node from old pull queue
      if (_blocks[move.to].pull.contains(u)) {
        _blocks[move.to].pull.remove(u);
      }

      insertOrUpdateNode(u, move.to);
    }

    //insert node into new pull queues and push queue
    //if the node was deleted and reinserted, the queues are updated to prevent duplicate entries
    void insertOrUpdateNode(HypernodeID u, PartitionID part) {
      ASSERT(partitioned_hypergraph_m->nodeIsEnabled(u));
      if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
        //TODO Why is this happening with mixed queries?
        return;
      }
      ASSERT(partitioned_hypergraph_m->partID(u) == part);
      HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();

      if (_blocks[part].push.get_positions_size() <= u)
      {
        reset();
      }

      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        if (b == part) {
          continue;
        }
        const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
        if (gain > highest_gain) {
          highest_gain = gain;
        }
        HyperedgeWeight weighted_pull_gain = (gain > 0) ? gain / partitioned_hypergraph_m->nodeWeight(u) : gain * partitioned_hypergraph_m->nodeWeight(u);
        ASSERT(u < _blocks[b].pull.get_positions_size());
        _blocks[b].pull.insertOrAdjustKey(u, weighted_pull_gain);
      }
      HyperedgeWeight weighted_push_gain = (highest_gain > 0) ? highest_gain * partitioned_hypergraph_m->nodeWeight(u) : highest_gain / partitioned_hypergraph_m->nodeWeight(u);
      ASSERT(u < _blocks[part].push.get_positions_size());
      _blocks[part].push.insertOrAdjustKey(u, weighted_push_gain);
    }

    bool checkBlockQueues() {
      for (HypernodeID u = 0; u < partitioned_hypergraph_m->initialNumNodes(); ++u) {
        if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
          continue;
        }
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == partitioned_hypergraph_m->partID(u)) {
            if (!_blocks[b].push.contains(u)) {
              std::cout << "Push queue does not contain node " << u << " in block " << b << std::endl;
              return false;
            }
          } else if (!_blocks[b].pull.contains(u)) {
            std::cout << "Pull queue does not contain node " << u << " in block " << b << std::endl;
            return false;
          }
        }
      }
      return true;
    }

    bool checkPullQueueGains() {
      for (HypernodeID u = 0; u < partitioned_hypergraph_m->initialNumNodes(); ++u) {
        if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
          continue;
        }
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == partitioned_hypergraph_m->partID(u)) {
            continue;
          }
          if (!_blocks[b].pull.contains(u)) {
            std::cout << "Pull queue does not contain node " << u << " in block " << b << std::endl;
            return false;
          }
          const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, partitioned_hypergraph_m->partID(u), b);
          const HyperedgeWeight weighted_gain = (gain > 0) ? gain / partitioned_hypergraph_m->nodeWeight(u) : gain * partitioned_hypergraph_m->nodeWeight(u);
          if (weighted_gain != _blocks[b].pull.keyOf(u)) {
            std::cout << "Pull queue gain for node " << u << " in block " << b << " does not match gain cache" << std::endl;
            std::cout << "Pull queue gain: " << _blocks[b].pull.keyOf(u) << " Gain cache: " << weighted_gain << std::endl;
            std::cout << "Gain: " << gain << std::endl;
            std::cout << "Node weight: " << partitioned_hypergraph_m->nodeWeight(u) << std::endl;
            return false;
          }
        }
      }
      return true;
    }

    //update gain cache for all changed nodes
    //account for the fact that later moves in the vector affect thresholds of earlier moves
    void updateGainForMoves(const std::vector<Move>& moves) {
      for (size_t i = 0; i < moves.size(); ++i) {
        Move move = moves[i];
        HypernodeID hn = move.node;
        ds::MutableHypergraph& hypergraph = partitioned_hypergraph_m->hypergraph();

        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          int nodes_in_removed_partition_post_removal = partitioned_hypergraph_m->pinCountInPart(he, move.from);
          int nodes_in_added_partition =partitioned_hypergraph_m->pinCountInPart(he,move.to);
          for (size_t later_move_ID = 0; later_move_ID < moves.size(); ++later_move_ID) {
            Move later_move = moves[later_move_ID];
            if (later_move.from == move.to) {
              nodes_in_removed_partition_post_removal -= 1;
              nodes_in_added_partition -= 1;
            }
            if (later_move.from == move.from) {
              nodes_in_removed_partition_post_removal -= 1;
              nodes_in_added_partition -= 1;
            }
            if (later_move.to == move.from) {
              nodes_in_removed_partition_post_removal -= 1;
              nodes_in_added_partition -= 1;
            }
            if (later_move.to == move.to) {
              nodes_in_removed_partition_post_removal -= 1;
              nodes_in_added_partition -= 1;
            }
          }
          if (nodes_in_removed_partition_post_removal <=1 || nodes_in_added_partition <= 2) {
            for (const HypernodeID& hn2 : hypergraph.pins(he)) {
              insertOrUpdateNode(hn2, partitioned_hypergraph_m->partID(hn2));
            }
          }
        }
      }
    }

private:

    ds::PartitionedHypergraph<ds::MutableHypergraph>* partitioned_hypergraph_m;
    gain_cache_t *_gain_cache;
    Context *_context;
    vec<Gain> *_benefit_aggregator;
    std::vector<uint32_t> shared_push_queue_locations;

//     struct BlockQueues {
// //        ds::Heap<Gain, HypernodeID> push;
//         ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> push;
//         ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> pull;
//
//         BlockQueues(uint32_t *, size_t positions_size):
//         push(positions_size), pull(positions_size) {}
//     };

  struct BlockQueues {
    ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> push;
    ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> pull;

    // prefer a size-only ctor (pointer parameter unused in your snippet)
    explicit BlockQueues(size_t positions_size)
      : push(positions_size), pull(positions_size) {}

    // disable copying (avoid shallow copies/double frees)
    BlockQueues(const BlockQueues&) = delete;
    BlockQueues& operator=(const BlockQueues&) = delete;

    // default move ctor/assign are fine if ExclusiveHandleHeap implements moves correctly
    BlockQueues(BlockQueues&&) noexcept = default;
    BlockQueues& operator=(BlockQueues&&) noexcept = default;
  };


    std::vector<BlockQueues> _blocks;

    std::tuple<PartitionID, HyperedgeWeight> findFittingBestPushBlock(HypernodeID u, PartitionID part) {
      HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
      PartitionID best_push_block = kInvalidPartition;
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
        if (gain > highest_gain && b != part && partitioned_hypergraph_m->partWeight(b) + partitioned_hypergraph_m->nodeWeight(u) <= _context->partition.max_part_weights[b]) {
          highest_gain = gain;
          best_push_block = b;
        }
      }
      return std::make_tuple(best_push_block, highest_gain);
    }

    void populateBlockQueues() {
      // shared_push_queue_locations = std::vector<uint32_t>(partitioned_hypergraph_m->initialNumNodes(), partitioned_hypergraph_m->initialNumNodes());
      // _blocks.reserve(_context->partition.k);
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        // _blocks.emplace_back(shared_push_queue_locations.data(), partitioned_hypergraph_m->initialNumNodes());
        _blocks.emplace_back(partitioned_hypergraph_m->initialNumNodes());
      }
      for (HypernodeID u = 0; u < partitioned_hypergraph_m->initialNumNodes(); ++u) {
        if (!partitioned_hypergraph_m->nodeIsEnabled(u)) {
          continue;
        }
        const PartitionID part = partitioned_hypergraph_m->partID(u);
        HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == part) {
            continue;
          }
          const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
          if (gain > highest_gain) {
            highest_gain = gain;
          }
          HyperedgeWeight weighted_pull_gain = (gain > 0) ? gain / partitioned_hypergraph_m->nodeWeight(u) : gain * partitioned_hypergraph_m->nodeWeight(u);
          _blocks[b].pull.insert(u, weighted_pull_gain);
        }
        HyperedgeWeight weighted_push_gain = (highest_gain > 0) ? highest_gain * partitioned_hypergraph_m->nodeWeight(u) : highest_gain / partitioned_hypergraph_m->nodeWeight(u);
        _blocks[part].push.insert(u, weighted_push_gain);
      }
      ASSERT(checkBlockQueues());
    }

    void updateGainCacheForMove(Move move) {

      HypernodeID hn = move.node;
      ds::MutableHypergraph& hypergraph = partitioned_hypergraph_m->hypergraph();

      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        size_t nodes_in_removed_partition_post_removal = partitioned_hypergraph_m->pinCountInPart(he, move.from);
        size_t nodes_in_added_partition =partitioned_hypergraph_m->pinCountInPart(he,move.to);
        if (nodes_in_removed_partition_post_removal == 0 || nodes_in_added_partition == 1 || nodes_in_added_partition == 2) {
          for (const HypernodeID& hn2 : hypergraph.pins(he)) {
            GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(
                    *partitioned_hypergraph_m, hn2, *_benefit_aggregator);
            insertOrUpdateNode(hn2, partitioned_hypergraph_m->partID(hn2));
          }
        } else if (nodes_in_removed_partition_post_removal == 1) {
          for (const HypernodeID &hn2: hypergraph.pins(he)) {
            if (hn2 != hn && partitioned_hypergraph_m->partID(hn2) == move.from) {
              GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(
                      *partitioned_hypergraph_m, hn2, *_benefit_aggregator);
              insertOrUpdateNode(hn2, partitioned_hypergraph_m->partID(hn2));
              break;
            }
          }
        }
      }

      GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(*partitioned_hypergraph_m, hn, *_benefit_aggregator);

      // ASSERT(partitioned_hypergraph_m->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(*_gain_cache)));
    }

};

} // namespace mt_kahypar
