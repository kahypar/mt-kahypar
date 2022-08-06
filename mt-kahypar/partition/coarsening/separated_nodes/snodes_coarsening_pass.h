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

namespace mt_kahypar {
namespace star_partitioning {

using ds::SeparatedNodes;
using ds::Array;

enum class SNodesCoarseningStage : uint8_t {
  PREFERABLE_DEGREE_ONE = 0,
  DEGREE_ONE_AND_TWINS = 1,
  DEGREE_ONE_AND_DEGREE_TWO_AND_TWINS = 2,
  // TODO: mixing with high degree??
  ANY_DEGREE = 3,
  ANY_DEGREE_RELAXED = 4,
  ANYTHING = 5
};

bool allowsDegreeTwo(const SNodesCoarseningStage& stage);

bool allowsHighDegree(const SNodesCoarseningStage& stage);

bool appliesTwins(const SNodesCoarseningStage& stage);

class SNodesCoarseningPass {
  struct NodeInfo {
    NodeInfo(): density(0.0), assigned_graph_node(kInvalidHypernode) { }

    double density;
    HypernodeID assigned_graph_node;
    // TODO: position in spanning tree
  };

  struct FullNodeInfo {
    FullNodeInfo():
                 node(kInvalidHypernode),
                 degree(0),
                 density(0),
                 assigned_graph_node(kInvalidHypernode)  { }

    FullNodeInfo(HypernodeID node, HyperedgeID degree, const NodeInfo& info):
                 node(node),
                 degree(degree),
                 density(info.density),
                 assigned_graph_node(info.assigned_graph_node)  { }

    HypernodeID node;
    HyperedgeID degree;
    double density;
    HypernodeID assigned_graph_node;
    // TODO: position in spanning tree
  };

  struct LocalizedData {
    tbb::enumerable_thread_specific<HypernodeID> match_counter;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_one_nodes;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_two_nodes;
  };

  struct Params {
    double accepted_density_diff;
    HypernodeWeight max_node_weight;
    HypernodeID degree_one_cluster_size;
  };

  struct Footprint;
  struct SimilarityHash;
  struct EqualityHash;

  // tuning constants
  static const HypernodeID DEGREE_ZERO_CLUSTER_SIZE = 64;
  static const HypernodeID MAX_CLUSTER_SIZE = 4;
  static constexpr size_t NUM_HASHES = 3;
  static constexpr double PREFERRED_DENSITY_DIFF = 1.6;
  static constexpr double TOLERABLE_DENSITY_DIFF = 2.1;
  static constexpr double RELAXED_DENSITY_DIFF = 4.1;

 public:
  SNodesCoarseningPass(const Hypergraph& hg, const Context& context,
                       const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage);

  // note: communities are allowed to be an empty vec
  void run(vec<HypernodeID>& communities);

  // only for testing
  const Array<HypernodeID>& nodeInfoBegin() const {
    return _node_info_begin;
  }

  // only for testing
  const Array<FullNodeInfo>& nodeInfo() const {
    return _node_info;
  }

 private:
  void setupNodeInfo();

  HypernodeID runCurrentStage(vec<HypernodeID>& communities);

  void applyDegreeZeroCoarsening(const Params& params, vec<HypernodeID>& communities, LocalizedData& data);

  void applyCoarseningForNode(const Params& params, vec<HypernodeID>& communities,
                              LocalizedData& data, const HypernodeID& node);

  void sortByDensity(vec<HypernodeID>& nodes);
  
  template<typename Hasher>
  void applyHashingRound(const Params& params, vec<HypernodeID>& communities,
                         LocalizedData& data, const HypernodeID min_degree) {
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
        while (pos + 1 < bucket.size()) {
          const auto& [curr_index, curr_hash] = bucket[pos];
          const double curr_density = info(curr_index).density;
          const auto& [next_index, next_hash] = bucket[pos + 1];
          const double next_density = info(next_index).density;
          if (curr_hash == next_hash && next_density / curr_density <= params.accepted_density_diff
              && _s_nodes.nodeWeight(info(curr_index).node) + _s_nodes.nodeWeight(info(next_index).node) <= params.max_node_weight) {
            // TODO: check actual similarity?!
            const HypernodeID curr_node = info(curr_index).node;
            communities[curr_node] = curr_node;
            communities[info(next_index).node] = curr_node;
            ++counter;
            ++pos;
          }
          ++pos;
        }
      }
    });
  }

  const FullNodeInfo& info(HypernodeID index) const {
    return _node_info[index];
  }

  const Hypergraph& _hg;
  const Context& _context;
  const SeparatedNodes& _s_nodes;
  Array<HypernodeID> _node_info_begin;
  Array<FullNodeInfo> _node_info;
  HypernodeID _current_num_nodes;
  HypernodeID _target_num_nodes;
  SNodesCoarseningStage _stage;
};

} // namepace star_partitioning
} // namespace mt_kahypar
