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

namespace {
  class NeighborhoodData {
   public:
    static constexpr size_t MAX_NEIGHBORHOOD_SIZE = 32;

    NeighborhoodData() {
      reset();
    }

    void insert(HypernodeID hn, float weight) {
      for (auto& entry: _data) {
        if (entry.first == kInvalidHypernode) {
          entry.first = hn;
          entry.second = weight;
          return;
        }
      }
      ASSERT(false, "capacity overflow!");
    }

    float get(HypernodeID hn) const {
      for (auto& entry: _data) {
        if (entry.first == hn) {
          return entry.second;
        }
      }
      return 0;
    }

    void reset() {
      _data.fill(std::make_pair(kInvalidHypernode, 0.0));
    }

    std::array<std::pair<HypernodeID, float>, MAX_NEIGHBORHOOD_SIZE> _data;
  };
} // namespace

class TwoHopClustering {
  using IncidenceMap = ds::SparseMap<HypernodeID, float>;

  static constexpr size_t TWIN_MATCHING_MAX_ATTEMPTS = 8;
  // degree threshold where it is extremely unlikely that two-hop coarsening is applicable
  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = 500;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

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
    }) {
      if (std::floor(1.0 / _context.coarsening.twin_min_relative_connectivity)
          > NeighborhoodData::MAX_NEIGHBORHOOD_SIZE) {
        ERR("Value for twin-min-relative-connectivity too small, must be at least"
             << (1.0 / NeighborhoodData::MAX_NEIGHBORHOOD_SIZE));
      }
      if (_context.coarsening.degree_one_node_cluster_size < 2) {
        ERR("Value for c-degree-one-node-cluster-size too small, must be at least 2");
      }
    }

  TwoHopClustering(const TwoHopClustering&) = delete;
  TwoHopClustering(TwoHopClustering&&) = delete;

  TwoHopClustering & operator= (const TwoHopClustering &) = delete;
  TwoHopClustering & operator= (TwoHopClustering &&) = delete;

  template<bool has_fixed_vertices, typename Hypergraph, typename DegreeSimilarityPolicy, typename F>
  void performClustering(const Hypergraph& hg,
                         const parallel::scalable_vector<HypernodeID>& node_mapping,
                         const DegreeSimilarityPolicy& similarity_policy,
                         ClusteringContext<Hypergraph>& cc,
                         F weight_ratio_for_node_fn,
                         int pass_nr = 0) {
    ASSERT(_context.coarsening.twin_required_similarity >= 0.5);
    _degree_one_map.reserve_for_estimated_number_of_insertions(cc.currentNumNodes() / 3);
    _twins_map.reserve_for_estimated_number_of_insertions(cc.currentNumNodes() / 3);

    auto fill_incidence_map_for_node = [&](IncidenceMap& incidence_map, const HypernodeID hn) {
      // TODO: can we do this more efficiently for graphs?
      HyperedgeWeight incident_weight_sum = 0;
      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        incident_weight_sum += hg.edgeWeight(he);
        for (const HypernodeID& pin: hg.pins(he)) {
          if (pin != hn) {
            HypernodeID target_cluster = cc.clusterID(pin);
            ASSERT(target_cluster != cc.clusterID(hn));  // holds since we only consider unmatched nodes
            incidence_map[target_cluster] += static_cast<double>(hg.edgeWeight(he)) / (hg.edgeSize(he) - 1);
          }
        }
      }
      return incident_weight_sum;
    };

    // insert degree one nodes and candidates for twins into buckets
    const double required_similarity = cc.contract_aggressively ?
        _context.coarsening.twin_reduced_required_similarity : _context.coarsening.twin_required_similarity;
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < node_mapping.size());
      const HypernodeID hn = node_mapping[id];
      if (hg.nodeIsEnabled(hn) && cc.vertexIsUnmatched(hn)
          && hg.nodeDegree(hn) <= HIGH_DEGREE_THRESHOLD) {
        IncidenceMap& incidence_map = _local_incidence_map.local();
        const HyperedgeWeight incident_weight_sum = fill_incidence_map_for_node(incidence_map, hn);

        HypernodeID key = 0;
        HyperedgeWeight considered_connectivity = 0;
        bool is_degree_one_node = false;
        for (const auto& entry: incidence_map) {
          const HypernodeID target_cluster = entry.key;
          const double connectivity = entry.value;
          if (connectivity >= required_similarity * incident_weight_sum) {
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
        if (!is_degree_one_node && incident_weight_sum > 0
            && considered_connectivity >= required_similarity * incident_weight_sum) {
          _twins_map.insert(key, MatchingEntry{key, hn});
        }
        incidence_map.clear();
      }
    });

    auto bucket_comparator = [&](const MatchingEntry& lhs, const MatchingEntry& rhs) {
      if (lhs.key == rhs.key) {
        const PartitionID community_lhs = hg.communityID(lhs.hn);
        const PartitionID community_rhs = hg.communityID(rhs.hn);
        return community_lhs < community_rhs ||
          (community_lhs == community_rhs && weight_ratio_for_node_fn(lhs.hn) < weight_ratio_for_node_fn(rhs.hn));
        } else {
          return lhs.key < rhs.key;
        }
    };
    auto accept_contraction = [&](const HypernodeID u, const HypernodeID v) {  // TODO: lazy evalution
      bool accept_community = cc.may_ignore_communities || (hg.communityID(u) == hg.communityID(v));
      bool weight_allowed = cc.clusterWeight(u) +  hg.nodeWeight(v) <= cc.max_allowed_node_weight;
      bool accept_similarity = similarity_policy.acceptContraction(hg, _context, u, v);
      // TODO fixed vertices
      return accept_community && weight_allowed && accept_similarity;
    };

    tbb::enumerable_thread_specific<HypernodeID> matched_d1_nodes(0);
    tbb::enumerable_thread_specific<HypernodeID> matched_twins(0);
    // match degree one nodes
    tbb::parallel_for(UL(0), _degree_one_map.numBuckets(), [&](const size_t bucket_id) {
      auto& bucket = _degree_one_map.getBucket(bucket_id);
      std::sort(bucket.begin(), bucket.end(), bucket_comparator);

      const HypernodeID max_size = _context.coarsening.degree_one_node_cluster_size;
      for (size_t i = 0; i + 1 < bucket.size() && cc.shouldContinue(); ++i) {
        HypernodeID num_matches = 0;
        for (size_t j = i + 1; j < i + max_size && j < bucket.size() && cc.shouldContinue()
             && bucket[i].key == bucket[j].key && accept_contraction(bucket[i].hn, bucket[j].hn); ++j) {
          ASSERT((j > i + 1 || cc.vertexIsUnmatched(bucket[i].hn)) && cc.vertexIsUnmatched(bucket[j].hn));
          bool success = cc.template matchVertices<has_fixed_vertices>(hg, bucket[i].hn, bucket[j].hn);
          if (success) {
            num_matches++;
          } else {
            break;
          }
        }
        i += num_matches;
        matched_d1_nodes.local() += num_matches;
      }
    });

    // match twins
    if (cc.shouldContinue()) {
      tbb::parallel_for(UL(0), _twins_map.numBuckets(), [&](const size_t bucket_id) {
        auto& bucket = _twins_map.getBucket(bucket_id);
        IncidenceMap& incidence_map = _local_incidence_map.local();
        std::sort(bucket.begin(), bucket.end(), bucket_comparator);

        for (size_t i = 0; i + 1 < bucket.size() && cc.shouldContinue(); ++i) {
          if (!cc.vertexIsUnmatched(bucket[i].hn)) {
            continue;
          }

          // prepare neighborhood data for current node
          const HyperedgeWeight incident_weight_sum_u = fill_incidence_map_for_node(incidence_map, bucket[i].hn);
          NeighborhoodData neighbor_data;
          for (const auto& entry: incidence_map) {
            const HypernodeID target_cluster = entry.key;
            const double connectivity = entry.value;
            if (connectivity >= _context.coarsening.twin_min_relative_connectivity * incident_weight_sum_u) {
              neighbor_data.insert(target_cluster, connectivity);
            }
          }
          incidence_map.clear();

          for (size_t j = i + 1; j < i + TWIN_MATCHING_MAX_ATTEMPTS + 1 && j < bucket.size()
              && bucket[i].key == bucket[j].key && accept_contraction(bucket[i].hn, bucket[j].hn); ++j) {
            if (!cc.vertexIsUnmatched(bucket[j].hn)) {
              continue;
            }

            // we compute a "relative" jaccard similarity of the two nodes
            const HyperedgeWeight incident_weight_sum_v = fill_incidence_map_for_node(incidence_map, bucket[j].hn);
            double relative_weight_of_intersection = 0;
            for (const auto& entry: incidence_map) {
              const double connectivity_v = entry.value;
              if (connectivity_v >= _context.coarsening.twin_min_relative_connectivity * incident_weight_sum_v) {
                double connectivity_u = neighbor_data.get(entry.key);
                if (connectivity_u > 0) {
                  relative_weight_of_intersection += std::min(connectivity_u / incident_weight_sum_u,
                                                              connectivity_v / incident_weight_sum_v);
                } else {
                  break;
                }
              }
            }
            ASSERT(relative_weight_of_intersection <= 1.0001);
            incidence_map.clear();

            // decide whether the neighborhood similarity is sufficiently high
            if (relative_weight_of_intersection >= required_similarity) {
              bool success = cc.template matchVertices<has_fixed_vertices>(hg, bucket[i].hn, bucket[j].hn);
              if (success) {
                matched_twins.local()++;
              }
            }
          }
        }
      });
    }

    _degree_one_map.clearParallel();
    _twins_map.clearParallel();

    if (_context.type == ContextType::main) {
      utils::Stats& stats = utils::Utilities::instance().getStats(_context.utility_id);
      auto report = [&](auto name, HypernodeID val) {
        std::stringstream ss;
        ss << "level_" << pass_nr << "_" << name;
        stats.add_stat<int32_t>(ss.str(), val);
      };
      report("matched_d1", matched_d1_nodes.combine(std::plus<HypernodeID>()));
      report("matched_twins", matched_twins.combine(std::plus<HypernodeID>()));
    }
  }

 private:
  const Context& _context;
  ds::ConcurrentBucketMap<MatchingEntry> _degree_one_map;
  ds::ConcurrentBucketMap<MatchingEntry> _twins_map;
  tbb::enumerable_thread_specific<IncidenceMap> _local_incidence_map;
};

}  // namespace mt_kahypar
