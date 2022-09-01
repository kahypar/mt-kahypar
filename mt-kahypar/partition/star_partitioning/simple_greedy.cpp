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

#include "simple_greedy.h"

#include <vector>
#include <algorithm>

#include "mt-kahypar/partition/star_partitioning/star_partitioning.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;

void SimpleGreedy::partition(PartitionedHypergraph& phg, SeparatedNodes& s_nodes, const Context& context,
                             Array<HypernodeWeight>& part_weights, bool parallel) {
  Array<HyperedgeWeight> max_gains;
  max_gains.assign(s_nodes.numNodes(), 0, parallel);
  Array<HypernodeID> sorted_nodes;
  sorted_nodes.assign(s_nodes.numNodes(), 0, parallel);

  auto set_node_gain = [&](const HypernodeID& node) {
    Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
    local_edge_weights.assign(phg.k(), 0, false);
    getEdgeWeightsOfNode(phg, s_nodes, local_edge_weights, node);

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
    tbb::parallel_for(ID(0), s_nodes.numNodes(), [&](const HypernodeID node) {
      set_node_gain(node);
    });
  } else {
    for (HypernodeID node = 0; node < s_nodes.numNodes(); ++node) {
      set_node_gain(node);
    }
  }
  
  auto compare = compare_gain_weight_ratio(
    [&](const HypernodeID& node) { return max_gains[node]; },
    [&](const HypernodeID& node) { return s_nodes.nodeWeight(node); }
  );
  if (parallel) {
    tbb::parallel_sort(sorted_nodes.begin(), sorted_nodes.end(), compare);
  } else {
    std::sort(sorted_nodes.begin(), sorted_nodes.end(), compare);
  }

  for (size_t i = 0; i < sorted_nodes.size(); ++i) {
    const HypernodeID node = sorted_nodes[i];
    if (phg.separatedPartID(node) == kInvalidPartition) {
      Array<HyperedgeWeight>& local_edge_weights = _tmp_edge_weights.local();
      local_edge_weights.assign(phg.k(), 0, false);
      getEdgeWeightsOfNode(phg, s_nodes, local_edge_weights, node);

      // greedily assign separated nodes
      PartitionID max_part = 0;
      HyperedgeWeight max_gain = 0;
      const std::vector<HypernodeWeight>& max_part_weights = context.partition.max_part_weights;
      for (PartitionID part = 0; part < _k; ++part) {
        if (local_edge_weights[part] >= max_gain &&
            part_weights[part] + s_nodes.nodeWeight(node) <= max_part_weights[part]) {
          max_part = part;
          max_gain = local_edge_weights[part];
        }
      }

      phg.separatedSetOnlyNodePart(node, max_part);
      part_weights[max_part] += s_nodes.nodeWeight(node);
    }
  }

}

} // namepace star_partitioning
} // namepace mt_kahypar

