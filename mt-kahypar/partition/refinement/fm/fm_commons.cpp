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

#include <algorithm>
#include <atomic>
#include <limits>

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"


namespace mt_kahypar {
  template<typename L, typename R>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  uint64_t pairToKey(L left, R right) {
    ASSERT(left >= 0 && static_cast<uint64_t>(left) <= std::numeric_limits<uint32_t>::max());
    ASSERT(right >= 0 && static_cast<uint64_t>(right) <= std::numeric_limits<uint32_t>::max());
    return (static_cast<uint64_t>(left) << 32) + static_cast<uint64_t>(right);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<uint32_t, uint32_t> keyToPair(uint64_t key) {
    return {key >> 32, key & std::numeric_limits<uint32_t>::max()};
  }

  Gain UnconstrainedFMData::estimatePenaltyForImbalancedMove(PartitionID to,
                                                             HNWeightConstRef initial_imbalance,
                                                             HNWeightConstRef moved_weight) const {
    ASSERT(initialized && to != kInvalidPartition);
    Gain result = 0;
    const Dimension dimension = moved_weight.dimension();
    for (Dimension d = 0; d < dimension; ++d) {
      const HNWeightScalar wu = moved_weight.at(d);
      BucketID bucketId = 0;
      while (bucketId < NUM_BUCKETS
            && initial_imbalance.at(d) + wu > bucket_weights[indexForBucket(to, d, bucketId)]) {
        ++bucketId;
      }

      if (bucketId >= NUM_BUCKETS) {
        // fallback case (it should be very unlikely that fallback_bucket_weights contains elements)
        const auto& fallback = fallback_bucket_weights[fallbackIndex(to, d)];
        while (bucketId < NUM_BUCKETS + fallback.size()
              && initial_imbalance.at(d) + wu > fallback[bucketId - NUM_BUCKETS]) {
          ++bucketId;
        }

        if (bucketId == NUM_BUCKETS + fallback.size()) {
          return std::numeric_limits<Gain>::max();  // invalid
        }
      }
      result += std::ceil(wu * gainPerWeightForBucket(bucketId));
    }
    return result;
  }


  template<typename GraphAndGainTypes>
  void UnconstrainedFMData::InitializationHelper<GraphAndGainTypes>::initialize(
            UnconstrainedFMData& data, const Context& context,
            const typename GraphAndGainTypes::PartitionedHypergraph& phg,
            const typename GraphAndGainTypes::GainCache& gain_cache) {
    auto get_node_stats = [&](const HypernodeID hypernode) {
      // TODO(maas): we might want to save the total incident weight in the hypergraph data structure
      // at some point in the future
      HyperedgeWeight total_incident_weight = 0;
      for (const HyperedgeID& he : phg.incidentEdges(hypernode)) {
        total_incident_weight += phg.edgeWeight(he);
      }
      HyperedgeWeight internal_weight = gain_cache.penaltyTerm(hypernode, phg.partID(hypernode));
      ASSERT(internal_weight == gain_cache.recomputePenaltyTerm(phg, hypernode));
      return std::make_pair(internal_weight, total_incident_weight);
    };

    const Dimension dimension = phg.dimension();
    const double bn_treshold = context.refinement.fm.treshold_border_node_inclusion;
    const auto weight_init_fn = [dimension]{ return AllocatedHNWeight(dimension, 0); };
    tbb::enumerable_thread_specific<AllocatedHNWeight> local_considered_weight(weight_init_fn);
    tbb::enumerable_thread_specific<AllocatedHNWeight> local_inserted_weight(weight_init_fn);

    // collect nodes and fill buckets
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      const HNWeightConstRef hn_weight = phg.nodeWeight(hn);
      if (weight::isZero(hn_weight)) return;

      auto [internal_weight, total_incident_weight] = get_node_stats(hn);
      if (static_cast<double>(internal_weight) >= bn_treshold * total_incident_weight) {
        local_considered_weight.local() += hn_weight;

        for (Dimension d = 0; d < dimension; ++d) {
          const BucketID bucketId = bucketForGainPerWeight(static_cast<double>(internal_weight) / hn_weight.at(d));
          if (bucketId < NUM_BUCKETS) {
            local_inserted_weight.local().at(d) += hn_weight.at(d);
            auto& local_weights = data.local_bucket_weights.local();
            local_weights[data.indexForBucket(phg.partID(hn), d, bucketId)] += hn_weight.at(d);
            data.rebalancing_nodes.set(hn, true);
          }
        }
      }
    });

    auto& bucket_weights = data.bucket_weights;
    // for each block compute prefix sum of bucket weights, which is later used for estimating penalties
    auto compute_prefix_sum_for_range = [&](size_t start, size_t end) {
      for (const auto& local_weights: data.local_bucket_weights) {
        ASSERT(bucket_weights.size() == local_weights.size());
        for (size_t i = start; i < end; ++i) {
          ASSERT(i < local_weights.size());
          bucket_weights[i] += local_weights[i];
        }
      }
      for (size_t i = start; i + 1 < end; ++i) {
        bucket_weights[i + 1] += bucket_weights[i];
      }
    };
    tbb::parallel_for(static_cast<PartitionID>(0), context.partition.k, [&](const PartitionID block) {
      for (Dimension d = 0; d < dimension; ++d) {
        const size_t first_bucket = data.indexForBucket(block, d, 0);
        compute_prefix_sum_for_range(first_bucket, first_bucket + NUM_BUCKETS);
      }
    }, tbb::static_partitioner());

    AllocatedHNWeight considered_weight(dimension, 0);
    for (const auto& val: local_considered_weight) {
      considered_weight += val;
    }
    AllocatedHNWeight inserted_weight(dimension, 0);
    for (const auto& val: local_inserted_weight) {
      inserted_weight += val;
    }

    for (Dimension d = 0; d < dimension; ++d) {
      if (static_cast<double>(inserted_weight.at(d)) / considered_weight.at(d) < FALLBACK_TRESHOLD) {
        // Use fallback if fixed number of buckets per block is not sufficient:
        // For unweighted instances or instances with reasonable weight distribution this should almost never
        // be necessary. We use more expensive precomputations (hash maps instead of arrays) here in order to
        // keep memory overhead low and still get fast queries for estimating imbalance penalties.
        using SparseMap = ds::DynamicSparseMap<uint64_t, HNWeightScalar>;

        // collect nodes into local hashmaps
        tbb::enumerable_thread_specific<SparseMap> local_accumulator;
        phg.doParallelForAllNodes([&](const HypernodeID hn) {
          const HNWeightScalar hn_weight = phg.nodeWeight(hn).at(d);
          if (hn_weight == 0) return;

          auto [internal_weight, total_incident_weight] = get_node_stats(hn);
          if (static_cast<double>(internal_weight) >= bn_treshold * total_incident_weight) {
            const BucketID bucketId = bucketForGainPerWeight(static_cast<double>(internal_weight) / hn_weight);
            if (bucketId >= NUM_BUCKETS) {
              auto& map = local_accumulator.local();
              // hash by block id and bucket id
              map[pairToKey(phg.partID(hn), bucketId - NUM_BUCKETS)] += hn_weight;
            }
          }
        });

        vec<AtomicBucketID> max_rank_per_block(context.partition.k, AtomicBucketID(0));
        vec<AtomicWeight> weight_per_block(context.partition.k, AtomicWeight(0));
        // sort resulting values and determine ranks (so that larger values are ignored)
        local_accumulator.combine_each([&](SparseMap& map) {
          // we sort the values in the dense part (note that this invalidates the map)
          std::sort(map.begin(), map.end(), [](const auto& l, const auto& r) {
            return l.key < r.key;
          });
          // compute summed weight and rank/bucketID for each block
          auto it = map.begin();
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            const uint32_t block = static_cast<uint32_t>(p);

            ASSERT(it == map.end() || keyToPair(it->key).first >= block);
            HNWeightScalar total_weight = 0;
            while (it < map.end() && keyToPair(it->key).first == block) {
              total_weight += it->value;
              ++it;
            }
            // scan backwards to find (approximately) an element with rank according to the fallback treshold
            HNWeightScalar remaining_upper_weight = std::floor((1.0 - FALLBACK_TRESHOLD) * total_weight);
            auto backwards_it = it;
            while (total_weight > 0 && remaining_upper_weight >= (--backwards_it)->value) {
              ASSERT(keyToPair(backwards_it->key).first == block);
              remaining_upper_weight -= backwards_it->value;
            }
            // write result to global arrays
            weight_per_block[block].fetch_add(total_weight, std::memory_order_relaxed);
            const auto [curr_block, new_rank] = keyToPair(backwards_it->key);
            if (curr_block == block) {
              AtomicBucketID& global_rank = max_rank_per_block[block];
              BucketID current = global_rank.load(std::memory_order_relaxed);
              while (current < new_rank
                    && !global_rank.compare_exchange_strong(current, new_rank, std::memory_order_relaxed)) { /* try again */ }
            } else {
              ASSERT(total_weight == 0);
            }
          }
        });

        // resize vectors accordingly, set rank to zero if no fallback is required for this block
        tbb::parallel_for(static_cast<PartitionID>(0), context.partition.k, [&](const PartitionID block) {
          const HNWeightScalar handled_weight = bucket_weights[data.indexForBucket(block, d, NUM_BUCKETS - 1)];
          const HNWeightScalar fallback_weight = weight_per_block[block];
          if (static_cast<double>(handled_weight) / (handled_weight + fallback_weight) >= FALLBACK_TRESHOLD) {
            max_rank_per_block[block].store(0);
          } else {
            data.fallbackWeights(block, d).resize(max_rank_per_block[block] + 1, 0);
          }
        }, tbb::static_partitioner());

        // accumulate results in fallback_bucket_weights
        local_accumulator.combine_each([&](SparseMap& map) {
          auto it = map.begin();
          for (PartitionID p = 0; p < context.partition.k; ++p) {
            const uint32_t block = static_cast<uint32_t>(p);
            const size_t upper_limit = data.fallbackWeights(block, d).size();
            ASSERT(upper_limit == 0 || upper_limit == max_rank_per_block[block] + 1);

            ASSERT(it == map.end() || keyToPair(it->key).first >= block);
            while (it < map.end() && keyToPair(it->key).first == block) {
              BucketID current_rank = keyToPair(it->key).second;
              if (current_rank < upper_limit) {
                std::atomic_ref(data.fallbackWeights(block, d)[current_rank])
                    .fetch_add(it->value, std::memory_order::relaxed);
              }
              ++it;
            }
          }
        });

        // compute prefix sums
        tbb::parallel_for(static_cast<PartitionID>(0), context.partition.k, [&](const PartitionID block) {
          auto& weights = data.fallbackWeights(block, d);
          if (!weights.empty()) {
            weights[0] += bucket_weights[data.indexForBucket(block, d, NUM_BUCKETS - 1)];
            for (size_t  i = 0; i + 1 < weights.size(); ++i) {
              weights[i + 1] += weights[i];
            }
          }
        }, tbb::static_partitioner());
      }
    }

    data.initialized = true;
  }

  void UnconstrainedFMData::reset() {
    // TODO: parallel
    rebalancing_nodes.reset();
    bucket_weights.assign(current_k * dimension * NUM_BUCKETS, 0);
    if (virtual_weight_delta.empty()) {
      virtual_weight_delta.replaceWith(current_k, dimension, 0, false);
    } else {
      virtual_weight_delta.assign(current_k, 0, false);
    }
    for (auto& local_weights: local_bucket_weights) {
      local_weights.assign(current_k * dimension * NUM_BUCKETS, 0);
    }
    fallback_bucket_weights.assign(current_k * dimension, {});
    initialized = false;
  }

  namespace {
  #define UNCONSTRAINED_FM_INITIALIZATION(X) UnconstrainedFMData::InitializationHelper<X>
  }

  INSTANTIATE_CLASS_WITH_VALID_TRAITS(UNCONSTRAINED_FM_INITIALIZATION)
}
