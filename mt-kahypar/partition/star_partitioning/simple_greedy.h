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

#include <vector>

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
  SimpleGreedy(Context context):
    _k(context.partition.k), _tmp_edge_weights(_k) { }

  template<typename F, typename G, typename H>
  void partition(const HypernodeID num_nodes, Array<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight> max_part_weights,
                 F get_edge_weights_of_node_fn, G get_node_weight_fn, H set_part_id_fn) {
    Array<HyperedgeWeight> max_gains(num_nodes);

    tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
      Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
      local_edge_weights.assign(_k, 0, false);
      get_edge_weights_of_node_fn(local_edge_weights, node);

      HyperedgeWeight max_gain = 0;
      for (PartitionID part = 0; part < _k; ++part) {
        if (local_edge_weights[part] >= max_gain) {
          max_gain = local_edge_weights[part];
        }
      }
      max_gains[node] = max_gain;
    });

    Array<HypernodeID> sorted_nodes(num_nodes);
    tbb::parallel_for(ID(0), num_nodes, [&](const HypernodeID node) {
      sorted_nodes[node] = node;
    });

    tbb::parallel_sort(sorted_nodes.begin(), sorted_nodes.end(),
      [&](const HypernodeID& left, const HypernodeID& right) {
        const HypernodeWeight weight_left = get_node_weight_fn(left);
        const HypernodeWeight weight_right = get_node_weight_fn(right);
        if (weight_left == 0) {
          return true;
        } else if (weight_right == 0) {
          return false;
        }
        return max_gains[left] / weight_left > max_gains[right] / weight_right;
      }
    );

    for (size_t i = 0; i < sorted_nodes.size(); ++i) {
      const HypernodeID node = sorted_nodes[i];
      Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
      local_edge_weights.assign(_k, 0);
      get_edge_weights_of_node_fn(local_edge_weights, node);

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

 private:
  PartitionID _k;
  tbb::enumerable_thread_specific<Array<HyperedgeWeight>> _tmp_edge_weights;
};

} // namepace star_partitioning
} // namepace mt_kahypar

