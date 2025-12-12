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

#include <array>

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/coarsening/multilevel/clustering_context.h"
#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"


namespace mt_kahypar {

class TwoHopClustering {
  using CacheEfficienIncidenceMap = ds::FixedSizeSparseMap<HypernodeID, RatingType>;

  struct MatchingEntry {
    HypernodeID key;
    HypernodeID hn;
  };

 public:
  TwoHopClustering(const HypernodeID num_nodes, const Context& context):
    _context(context),
    _favorite_clusters(),
    _local_incidence_map([=] {
      return CacheEfficienIncidenceMap(3UL * std::min(UL(num_nodes), _context.coarsening.two_hop_degree_threshold), 0.0);
    }) {
      if (_context.coarsening.two_hop_cluster_size < 2) {
        throw InvalidParameterException("Value for c-two-hop-cluster-size too small, must be at least 2");
      }
    }

  TwoHopClustering(const TwoHopClustering&) = delete;
  TwoHopClustering(TwoHopClustering&&) = delete;

  TwoHopClustering & operator= (const TwoHopClustering &) = delete;
  TwoHopClustering & operator= (TwoHopClustering &&) = delete;

  template<typename Hypergraph>
  void performClustering(const Hypergraph& hg,
                         const parallel::scalable_vector<HypernodeID>& node_mapping,
                         ClusteringContext<Hypergraph>& cc,
                         bool has_fixed_vertices) {
    _favorite_clusters.clearParallel();
    _favorite_clusters.reserve_for_estimated_number_of_insertions(cc.currentNumNodes() / 3);

    auto fill_incidence_map_for_node = [&](auto& incidence_map, const HypernodeID hn, bool& too_many_accesses) {
      // TODO: specialized version for graphs
      size_t num_accesses = 0;
      HyperedgeWeight incident_weight_sum = 0;
      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        const HypernodeID considered_pins = hg.edgeSize(he) - 1;
        if (num_accesses + considered_pins > _context.coarsening.two_hop_degree_threshold) {
          too_many_accesses = true;
          break;
        }

        for (const HypernodeID& pin: hg.pins(he)) {
          if (pin != hn) {
            HypernodeID target_cluster = cc.clusterID(pin);
            ASSERT(target_cluster != cc.clusterID(hn));  // holds since we only consider unmatched nodes
            incidence_map[target_cluster] += static_cast<double>(hg.edgeWeight(he)) / considered_pins;
          }
        }
        num_accesses += considered_pins;
        incident_weight_sum += hg.edgeWeight(he);
      }
      return incident_weight_sum;
    };

    // insert unmatched (low degree) nodes into buckets based on their favorite cluster
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < node_mapping.size());
      const HypernodeID hn = node_mapping[id];
      if (hg.nodeIsEnabled(hn) && cc.vertexIsUnmatched(hn)
          && hg.nodeWeight(hn) <= _context.coarsening.max_allowed_node_weight / 2
          && hg.nodeDegree(hn) <= _context.coarsening.two_hop_degree_threshold) {
        CacheEfficienIncidenceMap& incidence_map = _local_incidence_map.local();
        incidence_map.clear();

        bool too_many_accesses = false;
        const HyperedgeWeight incident_weight_sum = fill_incidence_map_for_node(incidence_map, hn, too_many_accesses);

        if (!too_many_accesses) {
          const float required_connectivity = _context.coarsening.two_hop_required_connectivity * incident_weight_sum;
          float max_connectivity = 0;
          HypernodeID best_target = kInvalidHypernode;
          for (const auto& [target_cluster, connectivity]: incidence_map) {
            if (connectivity >= required_connectivity && connectivity > max_connectivity) {
              max_connectivity = connectivity;
              best_target = target_cluster;
            }
          }
          if (best_target != kInvalidHypernode) {
            _favorite_clusters.insert(best_target, MatchingEntry{best_target, hn});
          }
        }
      }
    });

    auto bucket_comparator = [&](const MatchingEntry& lhs, const MatchingEntry& rhs) {
      if (lhs.key == rhs.key) {
        return hg.communityID(lhs.hn) < hg.communityID(rhs.hn);
      } else {
        return lhs.key < rhs.key;
      }
    };
    auto accept_community = [&](const HypernodeID u, const HypernodeID v) {
      return cc.may_ignore_communities || (hg.communityID(u) == hg.communityID(v));
    };

    // match nodes that have the same favorite cluster
    tbb::parallel_for(UL(0), _favorite_clusters.numBuckets(), [&](const size_t bucket_id) {
      auto& bucket = _favorite_clusters.getBucket(bucket_id);
      std::sort(bucket.begin(), bucket.end(), bucket_comparator);

      const HypernodeID max_size = _context.coarsening.two_hop_cluster_size;
      for (size_t i = 0; i + 1 < bucket.size() && cc.shouldContinue(); ++i) {
        HypernodeID num_matches = 0;
        for (size_t j = i + 1; j < i + max_size && j < bucket.size() && cc.shouldContinue()
             && bucket[i].key == bucket[j].key && accept_community(bucket[i].hn, bucket[j].hn); ++j) {
          ASSERT((j > i + 1 || cc.vertexIsUnmatched(bucket[i].hn)) && cc.vertexIsUnmatched(bucket[j].hn));
          // Note: cluster weight and fixed vertices are checked by `matchVertices`
          bool success = cc.matchVertices(hg, bucket[i].hn, bucket[j].hn, has_fixed_vertices);
          if (success) {
            num_matches++;
          }
        }
        i += num_matches;
      }
    });
  }

 private:
  const Context& _context;
  ds::ConcurrentBucketMap<MatchingEntry> _favorite_clusters;
  tbb::enumerable_thread_specific<CacheEfficienIncidenceMap> _local_incidence_map;
};

}  // namespace mt_kahypar
