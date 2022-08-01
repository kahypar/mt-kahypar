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

#include <vector>
#include <algorithm>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/partition/star_partitioning/simple_greedy.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;
using ds::StreamingVector;

class Approximate {
 public:
  Approximate(const PartitionID& k): _k(k) { }

  template<typename F, typename G, typename H>
  void partition(const HypernodeID num_nodes, Array<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight>& max_part_weights,
                 F get_edge_weights_of_node_fn, G get_node_weight_fn, H set_part_id_fn, parallel_tag_t) {
    Array<HyperedgeWeight> gains(num_nodes * _k);
    Array<PartitionID> preferred_part(num_nodes);

    tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
      get_edge_weights_of_node_fn(&gains[node * _k], node);

      HyperedgeWeight max_gain = 0;
      HyperedgeWeight min_gain = std::numeric_limits<HyperedgeWeight>::max();
      PartitionID max_part = 0;
      for (PartitionID part = 0; part < _k; ++part) {
        const HyperedgeWeight gain = gains[node * _k + part];
        min_gain = std::min(gain, min_gain);
        if (gain >= max_gain) {
          max_gain = gain;
          max_part = part;
        }
      }
      for (PartitionID part = 0; part < _k; ++part) {
        gains[node * _k + part] -= min_gain;
      }
      preferred_part[node] = max_part;
    });
    
    tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> node_index_in_part;
    tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> nodes;
    StreamingVector<HypernodeID> unassigned_stream;
    tbb::parallel_for(0, _k, [&](const PartitionID part) {
      auto get_gain = [&](HypernodeID node) { return gains[node * _k + part]; };

      auto& local_node_index_in_part = node_index_in_part.local();
      local_node_index_in_part.assign(num_nodes, 0);
      tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
        if (preferred_part[node] == part) {
          local_node_index_in_part[node] = 1;
        }
      });
      parallel::TBBPrefixSum<HypernodeID> index_prefix_sum(local_node_index_in_part);
      tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), num_nodes), index_prefix_sum);

      parallel::scalable_vector<HypernodeID>& local_nodes = nodes.local();
      local_nodes.clear();
      local_nodes.assign(index_prefix_sum.total_sum(), 0);
      tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
        if (preferred_part[node] == part) {
          local_nodes[index_prefix_sum[node]] = node;
        }
      });

      tbb::parallel_sort(local_nodes.begin(), local_nodes.end(),
                         compare_gain_weight_ratio(get_gain, get_node_weight_fn));
      std::vector<HypernodeID> excluded = minKnapsack(local_nodes, max_part_weights[part] - part_weights[part],
                                                      get_gain, get_node_weight_fn);

      size_t i = 0;
      for (const HypernodeID& node: local_nodes) {
        if (i < excluded.size() && node == excluded[i]) {
          unassigned_stream.stream(node);
          ++i;
        } else {
          set_part_id_fn(node, part);
          part_weights[part] += get_node_weight_fn(node);
        }
      }
      ALWAYS_ASSERT(i + 1 >= excluded.size());
    });

    // assign all currently unassigned nodes via the simple greedy algorithm
    SimpleGreedy sg(_k);
    auto unassigned = unassigned_stream.copy_parallel();
    sg.partition(unassigned.size(), part_weights, max_part_weights,
      [&](HyperedgeWeight* weights, const HypernodeID node) { get_edge_weights_of_node_fn(weights, unassigned[node]); },
      [&](const HypernodeID node) { return get_node_weight_fn(unassigned[node]); },
      [&](const HypernodeID node, const PartitionID part) { set_part_id_fn(unassigned[node], part); }
    );
  }

  template<typename F, typename G, typename H>
  void partition(const HypernodeID num_nodes, Array<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight>& max_part_weights,
                 F get_edge_weights_of_node_fn, G get_node_weight_fn, H set_part_id_fn) {
    Array<HyperedgeWeight> gains(num_nodes * _k);
    Array<PartitionID> preferred_part(num_nodes);
    Array<vec<HypernodeID>> nodes_per_part(_k);

    for (HypernodeID node = 0; node < num_nodes; ++node) {
      get_edge_weights_of_node_fn(&gains[node * _k], node);

      HyperedgeWeight max_gain = 0;
      HyperedgeWeight min_gain = std::numeric_limits<HyperedgeWeight>::max();
      PartitionID max_part = 0;
      for (PartitionID part = 0; part < _k; ++part) {
        const HyperedgeWeight gain = gains[node * _k + part];
        min_gain = std::min(gain, min_gain);
        if (gain >= max_gain) {
          max_gain = gain;
          max_part = part;
        }
      }
      for (PartitionID part = 0; part < _k; ++part) {
        gains[node * _k + part] -= min_gain;
      }
      nodes_per_part[max_part].push_back(node);
    }

    vec<HypernodeID> unassigned;
    for (PartitionID part = 0; part < _k; ++part) {
      auto get_gain = [&](HypernodeID node) { return gains[node * _k + part]; };

      vec<HypernodeID>& nodes = nodes_per_part[part];

      std::sort(nodes.begin(), nodes.end(), compare_gain_weight_ratio(get_gain, get_node_weight_fn));
      std::vector<HypernodeID> excluded = minKnapsack(nodes, max_part_weights[part] - part_weights[part],
                                                      get_gain, get_node_weight_fn);

      size_t i = 0;
      for (const HypernodeID& node: nodes) {
        if (i < excluded.size() && node == excluded[i]) {
          unassigned.push_back(node);
          ++i;
        } else {
          set_part_id_fn(node, part);
          part_weights[part] += get_node_weight_fn(node);
        }
      }
      ALWAYS_ASSERT(i + 1 >= excluded.size());
    }

    // assign all currently unassigned nodes via the simple greedy algorithm
    SimpleGreedy sg(_k);
    sg.partition(unassigned.size(), part_weights, max_part_weights,
      [&](HyperedgeWeight* weights, const HypernodeID node) { get_edge_weights_of_node_fn(weights, unassigned[node]); },
      [&](const HypernodeID node) { return get_node_weight_fn(unassigned[node]); },
      [&](const HypernodeID node, const PartitionID part) { set_part_id_fn(unassigned[node], part); }
    );
  }

 private:
  PartitionID _k;
};

template<typename F, typename G>
auto compare_gain_weight_ratio(F get_gain, G get_node_weight) {
  return [get_gain, get_node_weight](const HypernodeID& left, const HypernodeID& right) {
    const HypernodeWeight weight_left = get_node_weight(left);
    const HypernodeWeight weight_right = get_node_weight(right);
    if (weight_left == 0) {
      return false;
    } else if (weight_right == 0) {
      return true;
    }
    return static_cast<double>(get_gain(left)) / weight_left
            < static_cast<double>(get_gain(right)) / weight_right;
  };
}


/*
* Returns the nodes that are _not_ included in ascending order.
*/
template<typename F, typename G>
std::vector<HypernodeID> minKnapsack(const parallel::scalable_vector<HypernodeID>& sorted_nodes,
                                     const HypernodeWeight& capacity,
                                     F get_gain_fn, G get_node_weight_fn) {
  HypernodeWeight total_weight = 0;
  for (const HypernodeID& node: sorted_nodes) {
    total_weight += get_node_weight_fn(node);
  }
  std::vector<HypernodeID> result;
  const HypernodeWeight excluded_weight = total_weight - capacity;
  if (excluded_weight <= 0) {
    return result;
  }

  HyperedgeWeight min_gain = std::numeric_limits<HyperedgeWeight>::max();
  HypernodeID min_element = kInvalidHypernode;
  size_t min_num_elements = 0;
  HypernodeWeight current_weight = 0;
  HyperedgeWeight current_gain = 0;
  for (const HypernodeID& node: sorted_nodes) {
    const HypernodeWeight weight = get_node_weight_fn(node);
    const double gain = static_cast<double>(get_gain_fn(node));
    if (current_weight + weight <= excluded_weight) {
      current_weight += weight;
      current_gain += gain;
      result.push_back(node);
      if (current_weight == excluded_weight && current_gain <= min_gain) {
        return result;
      }
    } else if (weight == 0 || current_gain + (gain / weight) * (excluded_weight - current_weight)
               > static_cast<double>(min_gain)) {
      // bounding
      break;
    } else if (current_gain + gain < min_gain) {
      min_gain = current_gain + gain;
      min_element = node;
      min_num_elements = result.size();
    }
  }

  // reconstruct result
  result.resize(min_num_elements);
  result.push_back(min_element);
  return result;
}

} // namepace star_partitioning
} // namepace mt_kahypar

