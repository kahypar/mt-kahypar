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

#include "approximate.h"

#include <algorithm>

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/partition/star_partitioning/simple_greedy.h"
#include "mt-kahypar/partition/star_partitioning/star_partitioning.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;
using ds::StreamingVector;

void Approximate::partition(PartitionedHypergraph& phg, const SeparatedNodes& s_nodes, const Context& context,
                            Array<HypernodeWeight>& part_weights, parallel_tag_t) {
  // TODO(maas): there's a bug in the parallel implementation (and it's not too performant anyways)
  Array<HyperedgeWeight> gains(s_nodes.numNodes() * _k);
  Array<PartitionID> preferred_part(s_nodes.numNodes());

  tbb::parallel_for(ID(0), s_nodes.numNodes(), [&](const HypernodeID node) {
    getEdgeWeightsOfNode(phg, s_nodes, gains, node, node * _k);

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
  
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> node_index_in_part;
  parallel::scalable_vector<parallel::TBBPrefixSum<HypernodeID>> prefix_sums;
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> nodes;
  node_index_in_part.resize(_k);
  nodes.resize(_k);
  tbb::parallel_for(0, _k, [&](const PartitionID part) {
    node_index_in_part[part].assign(s_nodes.numNodes(), 0);
  });

  // mapping the indizes: we first setup the data... 
  tbb::parallel_for(ID(0), s_nodes.numNodes(), [&](const HypernodeID node) {
    node_index_in_part[preferred_part[node]][node] = 1;
  });
  for (PartitionID part = 0; part < _k; ++part) {
    prefix_sums.emplace_back(node_index_in_part[part]);
  }
  // ...and calculate the mapped index via a prefix sum 
  tbb::parallel_for(0, _k, [&](const PartitionID part) {
    tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), s_nodes.numNodes()), prefix_sums[part]);
    nodes[part].assign(prefix_sums[part].total_sum(), 0);
  });

  // then, we can copy the nodes to the vector of their respective part
  tbb::parallel_for(ID(0), s_nodes.numNodes(), [&](const HypernodeID node) {
    const PartitionID part = preferred_part[node];
    nodes[part][prefix_sums[part][node]] = node;
  });

  // now, we solve the knapsack problem for each part in parallel
  tbb::parallel_for(0, _k, [&](const PartitionID part) {
    auto get_gain = [&](HypernodeID node) { return gains[node * _k + part]; };
    auto get_node_weight_fn = [&](const HypernodeID& node) { return s_nodes.nodeWeight(node); };

    const std::vector<HypernodeWeight>& max_part_weights = context.partition.max_part_weights;
    parallel::scalable_vector<HypernodeID>& local_nodes = nodes[part];
    if (context.shared_memory.num_threads >= 2 * static_cast<size_t>(_k)) {
      tbb::parallel_sort(local_nodes.begin(), local_nodes.end(),
                          compare_gain_weight_ratio(get_gain, get_node_weight_fn));
    } else {
      std::sort(local_nodes.begin(), local_nodes.end(),
                compare_gain_weight_ratio(get_gain, get_node_weight_fn));
    }
    std::vector<HypernodeID> excluded = minKnapsack(local_nodes, max_part_weights[part] - part_weights[part],
                                                    get_gain, get_node_weight_fn);

    size_t i = 0;
    for (const HypernodeID& node: local_nodes) {
      if (i < excluded.size() && node == excluded[i]) {
        ++i;
      } else {
        phg.separatedSetOnlyNodePart(node, part);
        part_weights[part] += get_node_weight_fn(node);
      }
    }
    ALWAYS_ASSERT(i + 1 >= excluded.size());
  });

  // assign all currently unassigned nodes via the simple greedy algorithm
  SimpleGreedy sg(_k);
  sg.partition(phg, s_nodes, context, part_weights, true);
}

void Approximate::partition(PartitionedHypergraph& phg, const SeparatedNodes& s_nodes, const Context& context,
                            Array<HypernodeWeight>& part_weights) {
  vec<vec<HypernodeID>> nodes_per_part(_k);
  Array<HyperedgeWeight> gains;
  gains.assign(s_nodes.numNodes() * _k, 0, false);
  Array<PartitionID> preferred_part;
  preferred_part.assign(s_nodes.numNodes(), 0, false);

  for (HypernodeID node = 0; node < s_nodes.numNodes(); ++node) {
    getEdgeWeightsOfNode(phg, s_nodes, gains, node, node * _k);

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

  for (PartitionID part = 0; part < _k; ++part) {
    auto get_gain = [&](HypernodeID node) { return gains[node * _k + part]; };

    vec<HypernodeID>& nodes = nodes_per_part[part];

    auto get_node_weight_fn = [&](const HypernodeID& node) { return s_nodes.nodeWeight(node); };
    const std::vector<HypernodeWeight>& max_part_weights = context.partition.max_part_weights;
    std::sort(nodes.begin(), nodes.end(), compare_gain_weight_ratio(get_gain, get_node_weight_fn));
    std::vector<HypernodeID> excluded = minKnapsack(nodes, max_part_weights[part] - part_weights[part],
                                                    get_gain, get_node_weight_fn);

    size_t i = 0;
    for (const HypernodeID& node: nodes) {
      if (i < excluded.size() && node == excluded[i]) {
        ++i;
      } else {
        phg.separatedSetOnlyNodePart(node, part);
        part_weights[part] += get_node_weight_fn(node);
      }
    }
    ALWAYS_ASSERT(i + 1 >= excluded.size());
  }

  // assign all currently unassigned nodes via the simple greedy algorithm
  SimpleGreedy sg(_k);
  sg.partition(phg, s_nodes, context, part_weights, false);
}

} // namepace star_partitioning
} // namepace mt_kahypar

