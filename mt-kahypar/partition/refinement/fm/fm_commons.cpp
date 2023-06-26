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
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  double UnconstrainedFMData::estimatedPenaltyFromIndex(PartitionID to, size_t bucketId,
                                                        HypernodeWeight remaining) const {
    double penalty = 0;
    while (bucketId < NUM_BUCKETS && remaining > 0 && bucket_weights[indexForBucket(to, bucketId)] <= remaining) {
      const HypernodeWeight bucket_weight = bucket_weights[indexForBucket(to, bucketId)];
      penalty += gainPerWeightForBucket(bucketId) * bucket_weight;
      remaining -= bucket_weight;
      ++bucketId;
    }
    if (bucketId == NUM_BUCKETS && remaining > 0) {
      return static_cast<double>(std::numeric_limits<Gain>::max()); // rebalancing not possible with considered nodes
    }
    if (remaining > 0) {
      penalty += gainPerWeightForBucket(bucketId) * remaining;
    }
    return penalty;
  }

  Gain UnconstrainedFMData::estimatedPenaltyForImbalance(PartitionID to, HypernodeWeight total_imbalance) const {
    return std::ceil(estimatedPenaltyFromIndex(to, 0, total_imbalance));
  }

  Gain UnconstrainedFMData::estimatedPenaltyForDelta(PartitionID to, HypernodeWeight old_weight,
                                                     HypernodeWeight new_weight) const {
    ASSERT(old_weight <= new_weight);
    size_t bucketId = 0;
    HypernodeWeight common = 0;
    while (bucketId < NUM_BUCKETS && common + bucket_weights[indexForBucket(to, bucketId)] <= old_weight) {
      common += bucket_weights[indexForBucket(to, bucketId)];
      ++bucketId;
    }
    const double penalty_new = estimatedPenaltyFromIndex(to, bucketId, new_weight - common);
    const double penalty_old = estimatedPenaltyFromIndex(to, bucketId, old_weight - common);
    return std::ceil(penalty_new - penalty_old);
  }

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
      if (old >= remaining) {
        remaining = 0;
      } else if (old > 0) {
        // set consumed weight equal to zero (but use fetch_add to accomodate concurrent updates)
        consumed_bucket_weights[indexForBucket(to, bucketId)].fetch_add(
            remaining - old, std::memory_order_relaxed);
        remaining -= old;
      } else if (bucketId == 0) {
        // nothing to revert here
        return;
      }
      // else: try again
    }
  }

  template<typename PartitionedHypergraphT>
  void UnconstrainedFMData::precomputeForLevel(const PartitionedHypergraphT& phg) {
    disabled = false;
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      HyperedgeWeight incident_weight = 0;
      for (HyperedgeID he: phg.incidentEdges(hn)) {
        incident_weight += phg.edgeWeight(he);
      }
      incident_weight_of_node[hn] = incident_weight;
    });
  }

  template<typename PartitionedHypergraphT>
  void UnconstrainedFMData::initialize(const Context& context, const PartitionedHypergraphT& phg) {
    ASSERT(!initialized);
    changeNumberOfBlocks(context.partition.k);

    double bn_treshold = context.refinement.fm.treshold_border_node_inclusion;
    // collect nodes and fill buckets
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      const HypernodeWeight hn_weight = phg.nodeWeight(hn);
      if (hn_weight == 0) return;

      auto include_node = [&](const HyperedgeWeight used_weight) {
        const size_t bucketId = bucketForGainPerWeight(static_cast<double>(used_weight) / hn_weight);
        if (bucketId < NUM_BUCKETS) {
          auto& local_weights = local_bucket_weights.local();
          local_weights[indexForBucket(phg.partID(hn), bucketId)] += hn_weight;
          rebalancing_nodes.set(hn, true);
        }
      };

      if (bn_treshold == 1.0 && !phg.isBorderNode(hn)) {
        // TODO(maas): only non border nodes does not seem like a good strategy for hypergraphs
        include_node(incident_weight_of_node[hn]);
      } else if (bn_treshold < 1.0) {
        HyperedgeWeight internal_weight = 0;
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
          if (phg.connectivity(he) == 1) {
            internal_weight += phg.edgeWeight(he);
          }
        }
        if (static_cast<double>(internal_weight) >= bn_treshold * incident_weight_of_node[hn]) {
          include_node(internal_weight);
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
    ASSERT(upper_weight_limits.size() == static_cast<size_t>(context.partition.k));
    for (PartitionID block = 0; block < context.partition.k; ++block) {
      HypernodeWeight sum = 0;
      for (size_t bucketId = 0; bucketId < NUM_BUCKETS; ++ bucketId) {
        sum += bucket_weights[indexForBucket(block, bucketId)];
      }
      upper_weight_limits[block] = sum;
    }
    initialized = true;
  }

  namespace {
  #define UNCONSTRAINED_FM_PRECOMPUTE(X) void UnconstrainedFMData::precomputeForLevel(const X& phg)
  #define UNCONSTRAINED_FM_INITIALIZE(X) void UnconstrainedFMData::initialize(const Context& context, const X& phg)
  }

  INSTANTIATE_FUNC_WITH_PARTITIONED_HG(UNCONSTRAINED_FM_PRECOMPUTE)
  INSTANTIATE_FUNC_WITH_PARTITIONED_HG(UNCONSTRAINED_FM_INITIALIZE)

}
