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
  Gain UnconstrainedFMData::estimatePenaltyForImbalancedMove(PartitionID to,
                                                             HypernodeWeight initial_imbalance,
                                                             HypernodeWeight moved_weight) const {
    ASSERT(initialized);
    size_t bucketId = 0;
    while (initial_imbalance + moved_weight > bucket_weights[indexForBucket(to, bucketId)]
           && bucketId < NUM_BUCKETS) {
      ++bucketId;
      initial_imbalance -= bucket_weights[indexForBucket(to, bucketId)];
    }
    return (bucketId >= NUM_BUCKETS) ? std::numeric_limits<Gain>::max()
              : std::ceil(moved_weight * gainPerWeightForBucket(bucketId));
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

      HyperedgeWeight total_incident_weight = 0;
      HyperedgeWeight internal_weight = 0;
      for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        total_incident_weight += phg.edgeWeight(he);
        if (phg.connectivity(he) == 1) {
          internal_weight += phg.edgeWeight(he);
        }
      }
      if (static_cast<double>(internal_weight) >= bn_treshold * total_incident_weight) {
        const size_t bucketId = bucketForGainPerWeight(static_cast<double>(internal_weight) / hn_weight);
        if (bucketId < NUM_BUCKETS) {
          auto& local_weights = local_bucket_weights.local();
          local_weights[indexForBucket(phg.partID(hn), bucketId)] += hn_weight;
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

    // TODO: test this
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
