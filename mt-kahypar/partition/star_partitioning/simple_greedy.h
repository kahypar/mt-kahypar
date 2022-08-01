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
#include "mt-kahypar/definitions.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;

/*!
 * Simple star partitioning implementation that sorts all leaf nodes
 * and assigns them in a global greedy order.
*/
class SimpleGreedy {
 public:
  SimpleGreedy(const PartitionID& k): _k(k), _tmp_edge_weights(k) { }

  template<typename F, typename G, typename H>
  void partition(const HypernodeID num_nodes, Array<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight>& max_part_weights,
                 F get_edge_weights_of_node_fn, G get_node_weight_fn, H set_part_id_fn, bool parallel = true) {
    Array<HyperedgeWeight> max_gains(num_nodes);
    Array<HypernodeID> sorted_nodes(num_nodes);

    auto set_node_gain = [&](const HypernodeID& node) {
      Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
      local_edge_weights.assign(_k, 0, false);
      get_edge_weights_of_node_fn(local_edge_weights.data(), node);

      HyperedgeWeight max_gain = 0;
      for (PartitionID part = 0; part < _k; ++part) {
        if (local_edge_weights[part] >= max_gain) {
          max_gain = local_edge_weights[part];
        }
      }
      max_gains[node] = max_gain;
      sorted_nodes[node] = node;
    };

    if (parallel) {
      tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
        set_node_gain(node);
      });
    } else {
      for (HypernodeID node = 0; node < num_nodes; ++node) {
        set_node_gain(node);
      }
    }
    
    auto compare = compare_gain_weight_ratio(
      [&](const HypernodeID& node) { return max_gains[node]; }, get_node_weight_fn
    );
    if (parallel) {
      tbb::parallel_sort(sorted_nodes.begin(), sorted_nodes.end(), compare);
    } else {
      std::sort(sorted_nodes.begin(), sorted_nodes.end(), compare);
    }

    for (size_t i = 0; i < sorted_nodes.size(); ++i) {
      const HypernodeID node = sorted_nodes[i];
      Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
      local_edge_weights.assign(_k, 0);
      get_edge_weights_of_node_fn(local_edge_weights.data(), node);

      // greedily assign separated nodes
      PartitionID max_part = kInvalidPartition;
      HyperedgeWeight max_gain = 0;
      for (PartitionID part = 0; part < _k; ++part) {
        if (local_edge_weights[part] >= max_gain &&
            part_weights[part] + get_node_weight_fn(node) <= max_part_weights[part]) {
          max_part = part;
          max_gain = local_edge_weights[part];
        }
      }

      set_part_id_fn(node, max_part);
      part_weights[max_part] += get_node_weight_fn(node);
    }

  }

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

 private:
  PartitionID _k;
  tbb::enumerable_thread_specific<Array<HyperedgeWeight>> _tmp_edge_weights;
};

} // namepace star_partitioning
} // namepace mt_kahypar

