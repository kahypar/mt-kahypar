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

#pragma once

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"


namespace mt_kahypar {

class TwoHopClustering {
  using IncidenceMap = ds::SparseMap<HypernodeID, float>;

  struct MatchingEntry {
    HypernodeID key;
    HypernodeID hn;
  };

 public:
  TwoHopClustering(const HypernodeID num_nodes, const Context& context):
    _context(context),
    _degree_one_map(),
    _twins_map(),
    _local_incidence_map([=] {
      return IncidenceMap(num_nodes);
    }) { }

  TwoHopClustering(const TwoHopClustering&) = delete;
  TwoHopClustering(TwoHopClustering&&) = delete;

  TwoHopClustering & operator= (const TwoHopClustering &) = delete;
  TwoHopClustering & operator= (TwoHopClustering &&) = delete;

  template<bool has_fixed_vertices, typename Hypergraph, typename DegreeSimilarityPolicy>
  void performClustering(const Hypergraph& hg,
                         const parallel::scalable_vector<HypernodeID>& node_mapping,
                         const HypernodeID hierarchy_contraction_limit,
                         vec<HypernodeID>& cluster_ids,
                         MultilevelVertexPairRater& rater,
                         ConcurrentClusteringData& clustering_data,
                         NumNodesTracker& num_nodes_tracker,
                         ds::FixedVertexSupport<Hypergraph>& fixed_vertices,
                         const DegreeSimilarityPolicy& similarity_policy) {
    ASSERT(_context.coarsening.twin_required_similarity >= 0.5);
    _degree_one_map.reserve_for_estimated_number_of_insertions(num_nodes_tracker.currentNumNodes() / 3);
    _twins_map.reserve_for_estimated_number_of_insertions(num_nodes_tracker.currentNumNodes() / 3);
    _total_incident_weight.resize(hg.initialNumNodes());

    hg.doParallelForAllNodes([&](const HypernodeID hn) {  // can this be done on the fly?
      if (clustering_data.vertexIsUnmatched(hn)) {
        HyperedgeWeight incident_weight_sum = 0;
        for (const HyperedgeID& he: hg.incidentEdges(hn)) {
          incident_weight_sum += hg.edgeWeight(he);
        }
        _total_incident_weight[hn] = incident_weight_sum;
      }
    });

    // insert degree one nodes and candidates for twins into buckets
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < node_mapping.size());
      const HypernodeID hn = node_mapping[id];
      if (hg.nodeIsEnabled(hn) && clustering_data.vertexIsUnmatched(hn)) {
        // TODO: can we do this more efficiently for graphs?
        IncidenceMap& incidence_map = _local_incidence_map.local();
        for (const HyperedgeID& he : hg.incidentEdges(hn)) {
          for (const HypernodeID& pin: hg.pins(he)) {
            if (pin != hn) {
              PartitionID target_cluster = cluster_ids[hg.edgeTarget(he)];
              ASSERT(target_cluster != cluster_ids[hn]);  // holds since we only consider unmatched nodes
              incidence_map[target_cluster] += static_cast<double>(hg.edgeWeight(he)) / (hg.edgeSize(he) - 1);
            }
          }
        }

        HypernodeID key = 0;
        HyperedgeWeight considered_connectivity = 0;
        bool is_degree_one_node = false;
        const HypernodeWeight incident_weight_sum = _total_incident_weight[hn];
        for (const auto& entry: incidence_map) {
          const PartitionID target_cluster = entry.key;
          const double connectivity = entry.value;
          ASSERT(target_cluster >= 0 && target_cluster != kInvalidPartition);
          if (connectivity >= _context.coarsening.twin_required_similarity * incident_weight_sum) {
            // we consider this to be a degree one node
            _degree_one_map.insert(target_cluster, MatchingEntry{target_cluster, hn});
            is_degree_one_node = true;
            break;
          } else if (connectivity >= _context.coarsening.twin_min_relative_connectivity * incident_weight_sum) {
            // use sum of squares as hash value (must be commutative)
            key += target_cluster * target_cluster;
            considered_connectivity += connectivity;
          }
        }
        if (!is_degree_one_node && considered_connectivity >= _context.coarsening.twin_required_similarity * incident_weight_sum) {
          _twins_map.insert(key, MatchingEntry{key, hn});
        }
        incidence_map.clear();
      }
    });

    tbb::enumerable_thread_specific<HypernodeID> matched_d1_nodes(0);
    tbb::enumerable_thread_specific<HypernodeID> matched_twins(0);
    tbb::parallel_for(UL(0), _degree_one_map.numBuckets(), [&](const size_t bucket_id) {
      auto& bucket = _degree_one_map.getBucket(bucket_id);
      std::sort(bucket.begin(), bucket.end(), [&](const MatchingEntry& lhs, const MatchingEntry& rhs) {
        return lhs.key < rhs.key ||
          (lhs.key == rhs.key && similarity_policy.weightRatioForNode(hg, lhs.hn)
           < similarity_policy.weightRatioForNode(hg, rhs.hn));
      });
      for (size_t i = 0; i + 1 < bucket.size(); ++i) {
        // TODO: match more than 2 nodes?
        if (bucket[i].key == bucket[i + 1].key) {
          ASSERT(clustering_data.vertexIsUnmatched(bucket[i].hn)
                 && clustering_data.vertexIsUnmatched(bucket[i + 1].hn));
          bool success = clustering_data.template matchVertices<has_fixed_vertices>(
            hg, bucket[i].hn, bucket[i + 1].hn, cluster_ids, rater, fixed_vertices);
          if (success) {
            // update the number of nodes in a way that minimizes synchronization overhead
            num_nodes_tracker.subtractNode(_context.shared_memory.original_num_threads, hierarchy_contraction_limit);
            matched_d1_nodes.local()++;
            ++i;
          }
        }
      }
    });

    // tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
    //   ASSERT(id < node_mapping.size());
    //   const HypernodeID hn = node_mapping[id];
    //   if (hg.nodeIsEnabled(hn)
    //       && num_nodes_tracker.currentNumNodes() > hierarchy_contraction_limit
    //       && clustering_data.vertexIsUnmatched(hn)) {
    //     const Rating rating = rater.template rate<ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy, has_fixed_vertices>(
    //                                  hg, hn, cluster_ids, clustering_data.clusterWeight(), fixed_vertices,
    //                                  similarity_policy, _context.coarsening.max_allowed_node_weight);
    //     if (rating.target != kInvalidHypernode) {
    //       bool success = clustering_data.template matchVertices<has_fixed_vertices>(
    //         hg, hn, rating.target, cluster_ids, rater, fixed_vertices);
    //       if (success) {
    //         // update the number of nodes in a way that minimizes synchronization overhead
    //         num_nodes_tracker.subtractNode(_context.shared_memory.original_num_threads, hierarchy_contraction_limit);
    //       }
    //     }
    //   }
    // });
  }

 private:
  const Context& _context;
  parallel::scalable_vector<HyperedgeWeight> _total_incident_weight;
  ds::ConcurrentBucketMap<MatchingEntry> _degree_one_map;
  ds::ConcurrentBucketMap<MatchingEntry> _twins_map;
  tbb::enumerable_thread_specific<IncidenceMap> _local_incidence_map;
};

}  // namespace mt_kahypar
