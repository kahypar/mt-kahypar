/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <tbb/enumerable_thread_specific.h>

#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/spanning_tree.h"

namespace mt_kahypar {
namespace star_partitioning {

using ds::SeparatedNodes;
using ds::Array;

enum class SNodesCoarseningStage : uint8_t {
  D1_TWINS = 0,
  D1_D2_TWINS = 1,
  D1_D2_TWINS_SIMILARITY = 2,
  D1_D2_TWINS_SIMILARITY_RELAXED = 3,
  D1_D2_DX_TWINS_SIMILARITY = 4,
  ANYTHING = 5
};

class SNodesCoarseningPass {
  using HashFunc = kahypar::math::MurmurHash<HypernodeID>;
  using HashValue = typename HashFunc::HashValue;

  struct NodeInfo {
    NodeInfo(): node(kInvalidHypernode),
                degree(0),
                density(0),
                tree_path(kInvalidPath),
                assigned_graph_node(kInvalidHypernode)  { }

    HypernodeID node;
    HyperedgeID degree;
    double density;
    TreePath tree_path;
    HypernodeID assigned_graph_node;
  };

  struct LocalizedData {
    tbb::enumerable_thread_specific<HypernodeID> match_counter;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_one_nodes;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_two_nodes;
    tbb::enumerable_thread_specific<vec<HypernodeID>> high_degree_nodes;
  };

  struct Params {
    double accepted_density_diff;
    HypernodeWeight max_node_weight;
    HypernodeID degree_one_cluster_size;
    uint32_t max_tree_distance;
  };

  template<size_t NUM_HASHES>
  struct Footprint {
    bool operator==(const Footprint& other) const {
      for ( size_t i = 0; i < NUM_HASHES; ++i ) {
        if ( footprint[i] != other.footprint[i] ) {
          return false;
        }
      }
      return true;
    }

    bool operator<(const Footprint& other) const {
      for ( size_t i = 0; i < NUM_HASHES; ++i ) {
        if ( footprint[i] < other.footprint[i] ) {
          return true;
        } else if ( footprint[i] > other.footprint[i] ) {
          return false;
        }
      }
      return false;
    }

    HashValue footprint[NUM_HASHES];
  };

  // a similarity hashing function (see also hypergraph_sparsifier.h)
  template<size_t NUM_HASHES>
  struct SimilarityHash {
    static constexpr bool requires_check = true;
    using HashResult = Footprint<NUM_HASHES>;

    explicit SimilarityHash() {
      for (size_t i = 0; i < NUM_HASHES; ++i) {
        int seed = utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu());
        hash_functions[i] = HashFunc(seed);
      }
    }

    HashResult calculateHash(const vec<HypernodeID>& nodes) const {
      HashResult footprint;
      for ( size_t i = 0; i < NUM_HASHES; ++i ) {
        footprint.footprint[i] = minHash(hash_functions[i], nodes);
      }
      return footprint;
    }

    size_t combineHash(const HashResult& footprint) const {
      HashValue hash_value = kEdgeHashSeed;
      for ( const HashValue& value : footprint.footprint ) {
        hash_value ^= value;
      }
      return hash_value;
    }

  private:
    HashValue minHash(const HashFunc& hash_function,
                      const vec<HypernodeID>& nodes ) const {
      HashValue hash_value = std::numeric_limits<HashValue>::max();
      for ( const HypernodeID& node : nodes ) {
        hash_value = std::min(hash_value, hash_function(node));
      }
      return hash_value;
    }

    HashFunc hash_functions[NUM_HASHES];
  };

  struct EqualityHash;

  // tuning constants
  static const HypernodeID DEGREE_ZERO_CLUSTER_SIZE = 64;
  static const HypernodeID MAX_CLUSTER_SIZE = 4;
  static constexpr size_t NUM_SIMILARITY_ROUNDS = 3;
  static constexpr size_t SIMILARITY_SEARCH_RANGE = 8;
  static constexpr double SIMILARITY_RATIO = 0.4;
  static constexpr double SIMILARITY_RATIO_RELAXED = 0.25;
  static constexpr size_t D2_SEARCH_RANGE = 3;

 public:
  SNodesCoarseningPass(const Hypergraph& hg, const Context& context,
                       const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage);

  // note: communities are allowed to be an empty vec
  void run(vec<HypernodeID>& communities);

  SNodesCoarseningStage stage() const {
    return _stage;
  }

  // only for testing
  const Array<HypernodeID>& nodeInfoBegin() const {
    return _node_info_begin;
  }

  // only for testing
  const Array<NodeInfo>& nodeInfo() const {
    return _node_info;
  }

  // public for testing
  static std::pair<HyperedgeWeight, HyperedgeWeight> intersection_and_union(vec<std::pair<HypernodeID, HyperedgeWeight>>& lhs,
                                                                            vec<std::pair<HypernodeID, HyperedgeWeight>>& rhs);

 private:
  void setupNodeInfo();

  HypernodeID runCurrentStage(vec<HypernodeID>& communities, bool first = false);

  void applyDegreeZeroCoarsening(const Params& params, vec<HypernodeID>& communities, LocalizedData& data);

  void applyCoarseningForNode(const Params& params, vec<HypernodeID>& communities,
                              LocalizedData& data, const HypernodeID& node);

  void sortByDensity(vec<HypernodeID>& nodes);

  void matchPairs(const Params& params, vec<HypernodeID>& communities,
                  vec<HypernodeID>& nodes, HypernodeID& counter);
  
  template<typename Hasher>
  void applyHashingRound(const Params& params, vec<HypernodeID>& communities,
                         LocalizedData& data, const HypernodeID min_degree,
                         const double required_similarity = 1.0) {
    Hasher hasher;
    ds::ConcurrentBucketMap<std::pair<HypernodeID, typename Hasher::HashResult>> bucket_map;
    tbb::enumerable_thread_specific<vec<HypernodeID>> adjacent_nodes;
    bucket_map.reserve_for_estimated_number_of_insertions(_node_info.size()); // TODO: more precise estimation?

    tbb::parallel_for(0UL, _node_info.size(), [&](const size_t& index) {
      const HypernodeID node = info(index).node;
      if (communities[node] == kInvalidHypernode && info(index).degree >= min_degree) {
        vec<HypernodeID>& local_adjacent_nodes = adjacent_nodes.local();
        local_adjacent_nodes.clear();
        for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(node)) {
          local_adjacent_nodes.push_back(e.target);
        }
        typename Hasher::HashResult hash = hasher.calculateHash(local_adjacent_nodes);
        bucket_map.insert(hasher.combineHash(hash), {index, hash});
      }
    });

    tbb::parallel_for(0UL, bucket_map.numBuckets(), [&](const size_t i) {
      auto& bucket = bucket_map.getBucket(i);
      if (bucket.size() > 1) {
        // sort by hash and afterwards by density
        std::sort(bucket.begin(), bucket.end(), [&](const auto& left, const auto& right) {
          if (left.second == right.second) {
            return info(left.first).density < info(right.first).density;
          }
          return left.second < right.second;
        });

        // calculate pair-wise matching
        size_t pos = 0;
        HypernodeID& counter = data.match_counter.local();
        vec<std::pair<HypernodeID, HyperedgeWeight>> buffer_left;
        vec<std::pair<HypernodeID, HyperedgeWeight>> buffer_right;
        while (pos + 1 < bucket.size()) {
          const auto [curr_index, curr_hash] = bucket[pos];
          const double curr_density = info(curr_index).density;

          size_t search_offset = 1;
          auto [next_index, next_hash] = bucket[pos + 1];
          double next_density = info(next_index).density;
          while (curr_hash == next_hash && next_density / curr_density <= params.accepted_density_diff
              && search_offset <= (Hasher::requires_check ? SIMILARITY_SEARCH_RANGE : 1)
              && pos + search_offset < bucket.size()) {
            const HypernodeID curr_node = info(curr_index).node;
            const HypernodeID next_node = info(next_index).node;
            bool matches = _s_nodes.nodeWeight(curr_node) + _s_nodes.nodeWeight(next_node) <= params.max_node_weight;
            if (Hasher::requires_check) {
              auto intersec_union = intersection_and_union(curr_node, next_node, buffer_left, buffer_right);
              const double similarity_target = std::max(required_similarity,
                  1.5 / static_cast<double>(std::min(_s_nodes.inwardDegree(curr_node), _s_nodes.inwardDegree(next_node))));
              matches = static_cast<double>(intersec_union.first) / static_cast<double>(intersec_union.second) >= similarity_target;
            }
            if (matches) {
              ASSERT(communities[curr_node] == kInvalidHypernode && communities[next_node] == kInvalidHypernode);
              ASSERT(curr_node != next_node);
              communities[curr_node] = curr_node;
              communities[next_node] = curr_node;
              std::swap(bucket[pos + search_offset], bucket[pos + 1]);
              ++counter;
              ++pos;
              break;
            }

            ++search_offset;
            if (pos + search_offset < bucket.size()) {
              next_index = bucket[pos + search_offset].first;
              next_hash = bucket[pos + search_offset].second;
              next_density = info(next_index).density;
            }
          }
          ++pos;
        }
      }
    });
  }

  std::pair<HyperedgeWeight, HyperedgeWeight> intersection_and_union(const HypernodeID& s_node_left, const HypernodeID& s_node_right,
                                                                     vec<std::pair<HypernodeID, HyperedgeWeight>>& buffer_left,
                                                                     vec<std::pair<HypernodeID, HyperedgeWeight>>& buffer_right);

  const NodeInfo& info(HypernodeID index) const {
    return _node_info[index];
  }

  const Hypergraph& _hg;
  const Context& _context;
  const SeparatedNodes& _s_nodes;
  Array<HypernodeID> _node_info_begin;
  Array<NodeInfo> _node_info;
  HypernodeID _current_num_nodes;
  HypernodeID _target_num_nodes;
  SNodesCoarseningStage _stage;
};

} // namepace star_partitioning
} // namespace mt_kahypar
