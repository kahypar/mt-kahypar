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

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"



namespace mt_kahypar {

  Gain UnconstrainedFMData::estimatedPenaltyForImbalancedMove(PartitionID to, HypernodeWeight weight) const {
    ASSERT(initialized);
    size_t bucketId = 0;
    HyperedgeWeight free_capacity = 0;
    while (free_capacity <= 0 && bucketId < NUM_BUCKETS) {
      ++bucketId;
      free_capacity = (bucketId == NUM_BUCKETS) ? 0 : bucket_weights[indexForBucket(to, bucketId)] -
                      consumed_bucket_weights[indexForBucket(to, bucketId)].load(std::memory_order_relaxed);
    }
    if (free_capacity < weight) {
      ++bucketId;
    }
    return (bucketId >= NUM_BUCKETS) ? std::numeric_limits<Gain>::max()
              : std::ceil(weight * gainPerWeightForBucket(bucketId));
  }

  Gain UnconstrainedFMData::applyEstimatedPenaltyForImbalancedMove(PartitionID to, HypernodeWeight weight) {
    ASSERT(initialized);
    HypernodeWeight remaining = weight;
    double penalty = 0;
    while (remaining > 0) {
      size_t bucketId = 0;
      while (consumed_bucket_weights[indexForBucket(to, bucketId)].load(std::memory_order_relaxed)
              >= bucket_weights[indexForBucket(to, bucketId)]) {
        ++bucketId;
        if (bucketId == NUM_BUCKETS) {
          return std::numeric_limits<Gain>::max(); // rebalancing not possible with considered nodes
        }
      }

      const HypernodeWeight old = consumed_bucket_weights[indexForBucket(to, bucketId)].fetch_add(
                                      remaining, std::memory_order_relaxed);
      const HypernodeWeight max_weight = bucket_weights[indexForBucket(to, bucketId)];
      // might have been updated concurrently by another thread
      if (old + remaining <= max_weight) {
        penalty += remaining * gainPerWeightForBucket(bucketId);
        remaining = 0;
      } else if (old < max_weight) {
        // set consumed weight equal to max weight (but use fetch_sub to accomodate concurrent updates)
        consumed_bucket_weights[indexForBucket(to, bucketId)].fetch_sub(
            old + remaining - max_weight, std::memory_order_relaxed);
        penalty += (max_weight - old) * gainPerWeightForBucket(bucketId);
        remaining -= (max_weight - old);
      }
      // else: try again
    }
    return std::ceil(penalty);
  }

  void UnconstrainedFMData::revertImbalancedMove(PartitionID to, HypernodeWeight weight) {
    ASSERT(initialized);
    HypernodeWeight remaining = weight;
    while (remaining > 0) {
      size_t bucketId = 0;
      while (bucketId + 1 < NUM_BUCKETS && consumed_bucket_weights[
                indexForBucket(to, bucketId + 1)].load(std::memory_order_relaxed) > 0) {
        ++bucketId;
      }

      const HypernodeWeight old = consumed_bucket_weights[indexForBucket(to, bucketId)].fetch_sub(
                                      remaining, std::memory_order_relaxed);
      // might have been updated concurrently by another thread
      if (old - remaining >= 0) {
        remaining = 0;
      } else if (old > 0) {
        // set consumed weight equal to zero (but use fetch_add to accomodate concurrent updates)
        consumed_bucket_weights[indexForBucket(to, bucketId)].fetch_add(
            remaining - old, std::memory_order_relaxed);
        remaining -= old;
      }
      // else: try again
    }
  }

  template<typename PartitionedHypergraphT>
  void UnconstrainedFMData::initialize(const Context& context, const PartitionedHypergraphT& phg) {
    if (!initialized) {
      local_bucket_weights = tbb::enumerable_thread_specific<vec<HypernodeWeight>>(context.partition.k * NUM_BUCKETS);
    }
    if (buckets.empty()) {
      for (size_t i = 0; i < NUM_BUCKETS; ++i) {
        buckets.emplace_back(BUCKET_FACTOR);
      }
    } else {
      ASSERT(buckets.size() == NUM_BUCKETS);
      for (size_t i = 0; i < NUM_BUCKETS; ++i) {
        buckets[i].clearParallel();
      }
    }
    bucket_weights.assign(context.partition.k * NUM_BUCKETS, 0);
    consumed_bucket_weights.assign(context.partition.k * NUM_BUCKETS, AtomicWeight(0));
    for (auto& local_weights: local_bucket_weights) {
      local_weights.assign(context.partition.k * NUM_BUCKETS, 0);
    }

    // collect nodes and fill buckets
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      const HypernodeWeight hn_weight = phg.nodeWeight(hn);
      if (hn_weight > 0 && !phg.isBorderNode(hn)) {
        // TODO(maas): only non border nodes does not seem like a good strategy for hypergraphs
        const PartitionID from = phg.partID(hn);
        auto& local_weights = local_bucket_weights.local();
        HyperedgeWeight incident_weight = 0;
        for (HyperedgeID he: phg.incidentEdges(hn)) {
          // TODO: try using gain cache here instead
          incident_weight += phg.edgeWeight(he);
        }
        const size_t bucketId = bucketForGainPerWeight(static_cast<double>(incident_weight) / hn_weight);
        if (bucketId < NUM_BUCKETS) {
          buckets[bucketId].insert(hn, HypernodeID(hn));
          local_weights[indexForBucket(from, bucketId)] += hn_weight;
          rebalancing_nodes.set(hn, true);
        }
      }
    });

    // sum the bucket weight, so we know which buckets should be moved completely
    auto add_range_fn = [&](size_t start, size_t end) {
      for (const auto& local_weights: local_bucket_weights) {
        ASSERT(bucket_weights.size() == local_weights.size());
        for (size_t i = start; i < end; ++i) {
          ASSERT(i < local_weights.size());
          bucket_weights[i] += local_weights[i];
        }
      }
    };
    if (context.partition.k < 64) {
      add_range_fn(0, bucket_weights.size());
    } else {
      tbb::parallel_for(static_cast<PartitionID>(0), context.partition.k, [&](const PartitionID k) {
        add_range_fn(k * NUM_BUCKETS, (k + 1) * NUM_BUCKETS);
      });
    }
    initialized = true;
  }

  namespace {
  #define UNCONSTRAINED_FM_INITIALIZE(X) void UnconstrainedFMData::initialize(const Context& context, const X& phg)
  }

  INSTANTIATE_FUNC_WITH_PARTITIONED_HG(UNCONSTRAINED_FM_INITIALIZE)

}
