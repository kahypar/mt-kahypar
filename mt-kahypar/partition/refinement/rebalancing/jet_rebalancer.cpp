/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "mt-kahypar/partition/refinement/rebalancing/jet_rebalancer.h"

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

  template <typename TypeTraits, typename GainTypes>
  bool JetRebalancer<TypeTraits, GainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                        const vec<HypernodeID>&,
                                                        Metrics& best_metrics,
                                                        double) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    resizeDataStructuresForCurrentK();
    // If partition is imbalanced, rebalancer is activated
    bool improvement = false;
    if ( !metrics::isBalanced(phg, _context) ) {
      _gain.reset();
      initializeDataStructures(phg);

      rebalancingRound<true>(phg);

      // Update metrics statistics
      Gain delta = _gain.delta();
      // HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(_gain_cache));
      HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
        V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
      best_metrics.quality += delta;
      best_metrics.imbalance = metrics::imbalance(phg, _context);
      DBG << "[REBALANCE] " << V(delta) << "  imbalance=" << best_metrics.imbalance;
      improvement = delta < 0;
    }
    return improvement;
  }

  template <typename TypeTraits, typename GainTypes>
  template<bool ensure_balanced_moves>
  void JetRebalancer<TypeTraits, GainTypes>::rebalancingRound(PartitionedHypergraph& phg) {
    DBG << "[REBALANCE] Weights before rebalancing:";
    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      DBG << V(k) << "  weight=" << phg.partWeight(k) << "  max=" << _context.partition.max_part_weights[k];
    }
    insertNodesIntoBuckets(phg, [&](const HypernodeID hn) {
      return computeGainAndTargetPart(phg, hn, false).first;
    });

    processBuckets(phg, [&](const HypernodeID hn, const PartitionID from, bool is_retry) {
      auto [_, to] = computeGainAndTargetPart(phg, hn, true, is_retry);
      return changeNodePart(phg, hn, from, to, ensure_balanced_moves);
    }, ensure_balanced_moves, false);

    DBG << "[REBALANCE] Weights after rebalancing:";
    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      DBG << V(k) << "  weight=" << phg.partWeight(k) << "  max=" << _context.partition.max_part_weights[k];
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template<typename F>
  void JetRebalancer<TypeTraits, GainTypes>::insertNodesIntoBuckets(PartitionedHypergraph& phg, F compute_gain_fn) {
    parallel::scalable_vector<AtomicWeight> bucket_counts(NUM_BUCKETS); // TODO: remove

    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      const PartitionID from = phg.partID(hn);
      if (imbalance(from) > 0) {
        auto& local_weights = _local_bucket_weights.local();
        const size_t bucket = getBucketID(compute_gain_fn(hn));
        bucket_counts[bucket].fetch_add(1); // TODO: remove
        _buckets[bucket].insert(hn, HypernodeID(hn));
        local_weights[from * NUM_BUCKETS + bucket] += phg.nodeWeight(hn);
      }
    });

    // sum the bucket weight, so we know which buckets should be moved completely
    auto add_range_fn = [&](size_t start, size_t end) {
      for (const auto& local_weights: _local_bucket_weights) {
        ASSERT(_bucket_weights.size() == local_weights.size());
        for (size_t i = start; i < end; ++i) {
          ASSERT(i < local_weights.size());
          _bucket_weights[i] += local_weights[i];
        }
      }
    };
    if (_context.partition.k < 64) {
      add_range_fn(0, _bucket_weights.size());
    } else {
      tbb::parallel_for(static_cast<PartitionID>(0), _context.partition.k, [&](const PartitionID k) {
        add_range_fn(k * NUM_BUCKETS, (k + 1) * NUM_BUCKETS);
      });
    }
    for (size_t i = 0; i < bucket_counts.size(); ++i) {
      DBG << V(i) << "  num nodes in bucket: " << bucket_counts[i].load();
    }
  }

  template <typename TypeTraits, typename GainTypes>
  template<typename F>
  void JetRebalancer<TypeTraits, GainTypes>::processBuckets(PartitionedHypergraph& phg,
                                                            F move_node_fn,
                                                            bool retry_invalid_moves,
                                                            bool update_local_part_weights) {

    // process buckets in ascending order of loss value, moving the nodes to a block with free capacity
    tbb::enumerable_thread_specific<vec<HypernodeID>> retry_nodes;
    for (size_t bucket = 0; bucket < _buckets.size(); ++bucket) {
      auto process_current_node = [&](HypernodeID hn, PartitionID from, bool is_retry) {
        bool success = true;
        if (_bucket_weights[from * NUM_BUCKETS + bucket] < imbalance(from)) {
          // all nodes of the current bucket should be moved
          success = move_node_fn(hn, from, is_retry);
          if (update_local_part_weights && success) {
            _part_weights[from].fetch_sub(phg.nodeWeight(hn), std::memory_order_relaxed);
          }
        } else if (imbalance(from) > 0) {
          // the nodes of the current bucket will restore the balance of the block, thus we need to
          // directly apply the weight updates so that we don't remove more nodes than necessary
          const HypernodeWeight hn_weight = phg.nodeWeight(hn);
          const HypernodeWeight old_block_weight = _part_weights[from].fetch_sub(hn_weight, std::memory_order_relaxed);
          bool moved = false;
          if (old_block_weight > _context.partition.max_part_weights[from]) {
            // block is still imbalanced, try to move the node
            moved = move_node_fn(hn, from, is_retry);
            success = moved;
          }
          if (!moved) {
            _part_weights[from].fetch_add(hn_weight, std::memory_order_relaxed);
          }
        }
        return success;
      };

      BucketMap& map = _buckets[bucket];
      map.doParallelForAllBuckets([&](const size_t j) {
        vec<HypernodeID>& local_retry_nodes = retry_nodes.local();
        const auto& batch = map.getBucket(j);
        // handle all nodes in the bucket
        for (const HypernodeID hn: batch) {
          const PartitionID from = phg.partID(hn);
          bool success = process_current_node(hn, from, false);
          if (retry_invalid_moves && !success) {
            local_retry_nodes.push_back(hn);
          }
        }

        if (retry_invalid_moves) {
          for (const HypernodeID hn: local_retry_nodes) {
            const PartitionID from = phg.partID(hn);
            bool success = false;
            // retry ten times (in practice, first try should succeed in almost all cases)
            for (int i = 0; !success && i < 10; ++i) {
              success = process_current_node(hn, from, true);
            }
          }
          local_retry_nodes.clear();
        }
      });

      updateImbalance(phg, !update_local_part_weights);
      if (_num_imbalanced_blocks == 0) {
        break;
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  std::pair<Gain, PartitionID> JetRebalancer<TypeTraits, GainTypes>::computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
                                                                                              const HypernodeID hn,
                                                                                              bool non_adjacent_blocks,
                                                                                              bool use_precise_part_weights) {
    const HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
    RatingMap& tmp_scores = _gain.localScores();
    Gain isolated_block_gain = 0;
    _gain.precomputeGains(hypergraph, hn, tmp_scores, isolated_block_gain);

    Gain best_gain = isolated_block_gain;
    PartitionID best_target = kInvalidPartition;
    for (const auto& entry : tmp_scores) {
      const PartitionID to = entry.key;
      const Gain gain = _gain.gain(entry.value, isolated_block_gain);
      if (isValidTarget(hypergraph, to, hn_weight, use_precise_part_weights) && gain <= best_gain) {
        best_gain = gain;
        best_target = to;
      }
    }

    // if no adjacent block with free capacity exists, we need to consider non-adjacent blocks
    if (non_adjacent_blocks && best_target == kInvalidPartition) {
      utils::Randomize& rand = utils::Randomize::instance();
      // we start with a block that is chosen by random, to ensure a reasonable distribution of nodes
      // to target blocks (note: this does not always result in a uniform distribution since some blocks
      // are not an acceptable target, but it should be good enough)
      const PartitionID start = rand.getRandomInt(0, static_cast<int>(_context.partition.k - 1), SCHED_GETCPU);
      PartitionID to = start;
      do {
        if (!tmp_scores.contains(to) && isValidTarget(hypergraph, to, hn_weight, use_precise_part_weights)) {
          best_target = to;
          break;
        }

        ++to;
        if (to == _context.partition.k) {
          to = 0;
        }
      } while (to != start);
      // assertion does not always hold with tight balance constraint or large node weights
      // ASSERT(best_target != kInvalidPartition);
    }
    tmp_scores.clear();
    return std::make_pair(best_gain, best_target);
  }

  template <typename TypeTraits, typename GainTypes>
  Gain JetRebalancer<TypeTraits, GainTypes>::computeAverageGain(const PartitionedHypergraph& hypergraph,
                                                                const HypernodeID hn) {
    ASSERT(_num_valid_targets > 0);
    const HypernodeWeight hn_weight = hypergraph.nodeWeight(hn);
    RatingMap& tmp_scores = _gain.localScores();
    Gain isolated_block_gain = 0;
    _gain.precomputeGains(hypergraph, hn, tmp_scores, isolated_block_gain);

    Gain sum = 0;
    PartitionID connected_valid_targets = 0;
    for (const auto& entry : tmp_scores) {
      if (isValidTarget(hypergraph, entry.key, hn_weight, false)) {
        sum += _gain.gain(entry.value, isolated_block_gain);
        ++connected_valid_targets;
      }
    }
    tmp_scores.clear();

    PartitionID remaining_valid_targets = _num_valid_targets - connected_valid_targets;
    sum += remaining_valid_targets * isolated_block_gain;
    return sum / _num_valid_targets;
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRebalancer<TypeTraits, GainTypes>::updateImbalance(const PartitionedHypergraph& hypergraph,
                                                             bool read_weights_from_graph) {
    _num_imbalanced_blocks = 0;
    _num_valid_targets = 0;
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      if (read_weights_from_graph) {
        _part_weights[block] = hypergraph.partWeight(block);
      }
      if (imbalance(block) > 0) {
        ++_num_imbalanced_blocks;
      } else if (isValidTarget(hypergraph, block, 0, false)) {
        ++_num_valid_targets;
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  void JetRebalancer<TypeTraits, GainTypes>::initializeDataStructures(const PartitionedHypergraph& hypergraph) {
    ASSERT(_buckets.size() == NUM_BUCKETS && _part_weights.size() == static_cast<size_t>(_context.partition.k));
    updateImbalance(hypergraph);
    _bucket_weights.assign(_current_k * NUM_BUCKETS, 0);
    for (auto& local_weights: _local_bucket_weights) {
      local_weights.assign(_current_k * NUM_BUCKETS, 0);
    }

    size_t estimated_insertions = _num_imbalanced_blocks * hypergraph.initialNumNodes() / _context.partition.k;
    // buckets zero and one are probably mostly empty
    for (size_t i = 2; i < _buckets.size(); ++i) {
      estimated_insertions /= 2;
      _buckets[i].reserve_for_estimated_number_of_insertions(estimated_insertions);
    }
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  namespace {
  #define JET_REBALANCER(X, Y) JetRebalancer<X, Y>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(JET_REBALANCER)
}