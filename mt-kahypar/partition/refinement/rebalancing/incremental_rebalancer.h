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

    void init(ds::PartitionedHypergraph<ds::StaticHypergraph> &partitioned_hypergraph, Context &context, gain_cache_t &gain_cache, vec<Gain> &benefit_aggregator) {
      _partitioned_hypergraph_s = &partitioned_hypergraph;
      _context = &context;
      _gain_cache = &gain_cache;
      _benefit_aggregator = &benefit_aggregator;
      populateBlockQueues();
    }

    std::tuple<HyperedgeWeight, std::vector<HypernodeID>> rebalanceAndUpdateGainCache() {
      HyperedgeWeight total_gain = 0;
      std::vector<HypernodeID> moved_nodes;
      while (!mt_kahypar::metrics::isBalanced(*_partitioned_hypergraph_s, *_context)) {
        PartitionID imbalanced_block = kInvalidPartition;
        for (PartitionID i = 0; i < _context->partition.k; ++i) {
          if (_partitioned_hypergraph_s->partWeight(i) > _context->partition.max_part_weights[i]) {
            imbalanced_block = i;
            break;
          }
        }
        ASSERT(imbalanced_block != kInvalidPartition);
        ASSERT(!_blocks[imbalanced_block].push.empty());
        std::vector<std::tuple<PartitionID, HyperedgeWeight>> violations;
        while (!_blocks[imbalanced_block].push.empty()) {
          const HypernodeID u = _blocks[imbalanced_block].push.top();
          const HyperedgeWeight stored_gain = _blocks[imbalanced_block].push.topKey();
          ASSERT(_partitioned_hypergraph_s->partID(u) == imbalanced_block);
          _blocks[imbalanced_block].push.deleteTop();

          if (!_partitioned_hypergraph_s->nodeIsEnabled(u)) {
            continue;
          }

          auto [best_push_block, gain] = findFittingBestPushBlock(u, imbalanced_block);
          ASSERT(best_push_block != kInvalidPartition && best_push_block != imbalanced_block);

          //reassure that the gain has not changed
          const HyperedgeWeight weighted_gain = (gain > 0) ? gain * _partitioned_hypergraph_s->nodeWeight(u) : gain / _partitioned_hypergraph_s->nodeWeight(u);
          if (weighted_gain != stored_gain) {
            _blocks[imbalanced_block].push.insert(u, weighted_gain);
            continue;
          }

          Move move = {imbalanced_block, best_push_block, u, gain};
          executeMoveAndUpdateGainCache(move);
          applyMove(move);

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

        if (!_partitioned_hypergraph_s->nodeIsEnabled(u)) {
          continue;
        }
        PartitionID block_pull_source = _partitioned_hypergraph_s->partID(u);
        ASSERT(block_pull_source != block_pull_target);

        HyperedgeWeight gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, block_pull_source, block_pull_target);

        //reassure that the gain has not changed
        const HyperedgeWeight weighted_gain = (gain > 0) ? gain / _partitioned_hypergraph_s->nodeWeight(u) : gain * _partitioned_hypergraph_s->nodeWeight(u);
        if (weighted_gain != stored_gain) {
          _blocks[block_pull_target].pull.insert(u, weighted_gain);
          continue;
        }

        //prevent moves that violate max_part_weights
        if (_partitioned_hypergraph_s->partWeight(block_pull_target) + _partitioned_hypergraph_s->nodeWeight(u) >
            _context->partition.max_part_weights[block_pull_target]) {
          violations.push_back(std::make_tuple(u, weighted_gain));
          continue;
        }

        //reinsert skipped nodes
        for (const auto& violation : violations) {
          _blocks[block_pull_target].pull.insert(std::get<0>(violation), std::get<1>(violation));
        }

        Move move = {block_pull_source, block_pull_target, u, gain};
        executeMoveAndUpdateGainCache(move);
        applyMove(move);

        moved_nodes.push_back(u);

        total_gain += gain;
      }
      return std::make_tuple(total_gain, moved_nodes);
    }

    Move pop_pull_checked(PartitionID block_pull_target) {
      while (!_blocks[block_pull_target].pull.empty()) {
        const HypernodeID u = _blocks[block_pull_target].pull.top();
        if (!_partitioned_hypergraph_s->nodeIsEnabled(u)) {
          _blocks[block_pull_target].pull.deleteTop();
          continue;
        }
        const HyperedgeWeight stored_gain = _blocks[block_pull_target].pull.topKey();
        ASSERT(_partitioned_hypergraph_s->partID(u) != block_pull_target);
        _blocks[block_pull_target].pull.deleteTop();

        PartitionID block_pull_source = _partitioned_hypergraph_s->partID(u);

        HyperedgeWeight gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).benefitTerm(u, block_pull_target);

        if (gain != stored_gain) {
          _blocks[block_pull_target].pull.insert(u, gain);
          continue;
        }
        Move move = {block_pull_source, block_pull_target, u, gain};
        return move;
      }
    }

    void applyMove(Move move) {
      HypernodeID u = move.node;
      ASSERT(_partitioned_hypergraph_s->nodeIsEnabled(u));

      //remove node from old push queue
      if (_blocks[move.from].push.contains(u)) {
        _blocks[move.from].push.remove(u);
      }

      //remove node from old pull queue
      if (_blocks[move.to].pull.contains(u)) {
        _blocks[move.to].pull.remove(u);
      }

      insertNewNode(u, move.to);
    }

    //insert node into new pull queues and push queue
    //if the node was deleted and reinserted, the queues are updated to prevent duplicate entries
    void insertNewNode(HypernodeID u, PartitionID part) {
      ASSERT(_partitioned_hypergraph_s->nodeIsEnabled(u));
      HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        if (b == part) {
          continue;
        }
        const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
        if (gain > highest_gain) {
          highest_gain = gain;
        }
        HyperedgeWeight weighted_pull_gain = (gain > 0) ? gain / _partitioned_hypergraph_s->nodeWeight(u) : gain * _partitioned_hypergraph_s->nodeWeight(u);
        _blocks[b].pull.insertOrAdjustKey(u, weighted_pull_gain);
      }
      HyperedgeWeight weighted_push_gain = (highest_gain > 0) ? highest_gain * _partitioned_hypergraph_s->nodeWeight(u) : highest_gain / _partitioned_hypergraph_s->nodeWeight(u);
      _blocks[part].push.insertOrAdjustKey(u, weighted_push_gain);
    }

    bool checkBlockQueues() {
      for (HypernodeID u = 0; u < _partitioned_hypergraph_s->initialNumNodes(); ++u) {
        if (!_partitioned_hypergraph_s->nodeIsEnabled(u)) {
          continue;
        }
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == _partitioned_hypergraph_s->partID(u)) {
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

private:

    ds::PartitionedHypergraph<ds::StaticHypergraph> *_partitioned_hypergraph_s;
    gain_cache_t *_gain_cache;
    Context *_context;
    vec<Gain> *_benefit_aggregator;
    std::vector<uint32_t> shared_push_queue_locations;

    struct BlockQueues {
//        ds::Heap<Gain, HypernodeID> push;
        ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> push;
        ds::ExclusiveHandleHeap<ds::Heap<Gain, HypernodeID>> pull;

        BlockQueues(uint32_t * push_positions, size_t positions_size):
        push(positions_size), pull(positions_size) {}
    };

    std::vector<BlockQueues> _blocks;

    std::tuple<PartitionID, HyperedgeWeight> findFittingBestPushBlock(HypernodeID u, PartitionID part) {
      HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
      PartitionID best_push_block = kInvalidPartition;
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
        if (gain > highest_gain && b != part && _partitioned_hypergraph_s->partWeight(b) + _partitioned_hypergraph_s->nodeWeight(u) <= _context->partition.max_part_weights[b]) {
          highest_gain = gain;
          best_push_block = b;
        }
      }
      return std::make_tuple(best_push_block, highest_gain);
    }

    void populateBlockQueues() {
      shared_push_queue_locations = std::vector<uint32_t>(_partitioned_hypergraph_s->initialNumNodes(), _partitioned_hypergraph_s->initialNumNodes());
      for (PartitionID b = 0; b < _context->partition.k; ++b) {
        _blocks.emplace_back(shared_push_queue_locations.data(), _partitioned_hypergraph_s->initialNumNodes());
      }
      for (HypernodeID u = 0; u < _partitioned_hypergraph_s->initialNumNodes(); ++u) {
        if (!_partitioned_hypergraph_s->nodeIsEnabled(u)) {
          continue;
        }
        const PartitionID part = _partitioned_hypergraph_s->partID(u);
        HyperedgeWeight highest_gain = std::numeric_limits<HyperedgeWeight>::min();
        for (PartitionID b = 0; b < _context->partition.k; ++b) {
          if (b == part) {
            continue;
          }
          const Gain gain = GainCachePtr::cast<Km1GainCache>(*_gain_cache).gain(u, part, b);
          if (gain > highest_gain) {
            highest_gain = gain;
          }
          HyperedgeWeight weighted_pull_gain = (gain > 0) ? gain / _partitioned_hypergraph_s->nodeWeight(u) : gain * _partitioned_hypergraph_s->nodeWeight(u);
          _blocks[b].pull.insert(u, weighted_pull_gain);
        }
        HyperedgeWeight weighted_push_gain = (highest_gain > 0) ? highest_gain * _partitioned_hypergraph_s->nodeWeight(u) : highest_gain / _partitioned_hypergraph_s->nodeWeight(u);
        _blocks[part].push.insert(u, weighted_push_gain);
      }
      ASSERT(checkBlockQueues());
    }

    void executeMoveAndUpdateGainCache(Move move) {

      _partitioned_hypergraph_s->changeNodePart(move.node, move.from, move.to);
      HypernodeID hn = move.node;
      ds::StaticHypergraph &hypergraph = _partitioned_hypergraph_s->hypergraph();

      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        size_t nodes_in_removed_partition_post_removal = _partitioned_hypergraph_s->pinCountInPart(he, move.from);
        size_t nodes_in_added_partition =_partitioned_hypergraph_s->pinCountInPart(he,move.to);
        if (nodes_in_removed_partition_post_removal == 0 || nodes_in_added_partition == 1 || nodes_in_added_partition == 2) {
          for (const HypernodeID& hn2 : hypergraph.pins(he)) {
            GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(
                    *_partitioned_hypergraph_s, hn2, *_benefit_aggregator);
          }
        } else if (nodes_in_removed_partition_post_removal == 1) {
          for (const HypernodeID &hn2: hypergraph.pins(he)) {
            if (hn2 != hn && _partitioned_hypergraph_s->partID(hn2) == move.from) {
              GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(
                      *_partitioned_hypergraph_s, hn2, *_benefit_aggregator);
              break;
            }
          }
        }
      }

      GainCachePtr::cast<Km1GainCache>(*_gain_cache).initializeGainCacheEntryForNode(*_partitioned_hypergraph_s, hn, *_benefit_aggregator);

      ASSERT(_partitioned_hypergraph_s->checkTrackedPartitionInformation(GainCachePtr::cast<Km1GainCache>(*_gain_cache)));
    }

};

} // namespace mt_kahypar
