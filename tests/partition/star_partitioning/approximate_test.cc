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

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/star_partitioning/approximate.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"

#include <tbb/parallel_sort.h>

using ::testing::Test;

namespace mt_kahypar {
using ds::Array;
using ds::SeparatedNodes;
using ds::SepNodesStack;
using star_partitioning::Approximate;

class AMinKnapsack : public Test {
 public:
  void initialize(std::vector<HypernodeWeight>&& weights, std::vector<HyperedgeWeight>&& gains) {
    ASSERT(weights.size() == gains.size());
    _weights = std::move(weights);
    _gains = std::move(gains);

    _sorted_nodes.clear();
    _sorted_nodes.assign(_weights.size(), 0);
    tbb::parallel_for(0UL, _sorted_nodes.size(), [&](const size_t i) {
        _sorted_nodes[i] = i;
    });
    tbb::parallel_sort(_sorted_nodes.begin(), _sorted_nodes.end(),
      [&](const HypernodeID& left, const HypernodeID& right) {
        const HypernodeWeight weight_left = _weights[left];
        const HypernodeWeight weight_right = _weights[right];
        if (weight_left == 0) {
          return false;
        } else if (weight_right == 0) {
          return true;
        }
        return static_cast<double>(_gains[left]) / weight_left
               < static_cast<double>(_gains[right]) / weight_right;
      }
    );
  }

  void applyMinKnapsack(const HypernodeWeight& capacity, const std::vector<HypernodeID>& expected) {
    auto get_gain = [&](HypernodeID node) { return _gains[node]; };
    auto get_weight = [&](HypernodeID node) { return _weights[node]; };
    auto result = star_partitioning::minKnapsack(_sorted_nodes, capacity, get_gain, get_weight);

    ASSERT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        ASSERT_EQ(result[i], expected[i]);
    }
  }

 private:
  parallel::scalable_vector<HypernodeID> _sorted_nodes;
  std::vector<HypernodeWeight> _weights;
  std::vector<HyperedgeWeight> _gains;
};

class AnApproximate : public Test {

 public:
  AnApproximate() : graph(), context(), stack() {
  }

  void initialize(HypernodeID num_nodes, std::vector<HypernodeWeight> max_part_weights,
                  vec<vec<std::pair<HypernodeID, HyperedgeWeight>>> sep_nodes, vec<HypernodeWeight> node_weights = {}) {
    graph = HypergraphFactory::construct(num_nodes, 0, {});
    stack = SepNodesStack(num_nodes);
    graph.setSeparatedNodes(&stack);

    context.partition.k = max_part_weights.size();
    context.partition.max_part_weights = std::move(max_part_weights);

    SeparatedNodes& s_nodes = stack.onliest();
    vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> nodes;
    vec<SeparatedNodes::Edge> s_edges;
    size_t i = 0;
    for (const auto& edges_of_node: sep_nodes) {
      const HypernodeWeight weight = node_weights.empty() ? 1 : node_weights[i++];
      nodes.emplace_back(kInvalidHypernode, s_edges.size(), weight);
      for (const auto& [target, weight]: edges_of_node) {
        s_edges.emplace_back(target, weight);
      }
    }
    s_nodes.addNodes(graph, nodes, s_edges);
    s_nodes.revealAll();
    s_nodes.initializeOutwardEdges();

    phg = PartitionedHypergraph(context.partition.k, graph);
    phg.initializeSeparatedParts();
    phg.initializePartition();
  }

  void setup(vec<PartitionID> parts) {
    phg.resetPartition();
    ASSERT(parts.size() == phg.initialNumNodes());
    for (size_t i = 0; i < parts.size(); ++i) {
      phg.setOnlyNodePart(i, parts[i]);
    }
  }

  void partition(bool parallel) {
    Approximate ap(context.partition.k);
    Array<HypernodeWeight> part_weights;
    part_weights.assign(context.partition.k, 0, parallel);
    if (parallel) {
      ap.partition(phg, stack.onliest(), context, part_weights, parallel_tag_t());
    } else {
      ap.partition(phg, stack.onliest(), context, part_weights);
    }
  }

  PartitionedHypergraph phg;
  Hypergraph graph;
  Context context;
  SepNodesStack stack;
};

TEST_F(AMinKnapsack, EdgeCase) {
    initialize({}, {});
    applyMinKnapsack(0, {});
    applyMinKnapsack(10, {});

    initialize({1, 2}, {2, 1});
    applyMinKnapsack(0, {1, 0});
    applyMinKnapsack(3, {});
}

TEST_F(AMinKnapsack, SimpleCases) {
    initialize({4, 3, 2}, {5, 4, 3});
    applyMinKnapsack(2, {0, 1});
    applyMinKnapsack(3, {0, 2});
    applyMinKnapsack(4, {0, 2});
    applyMinKnapsack(5, {0});
    applyMinKnapsack(6, {1});
    applyMinKnapsack(7, {2});
    applyMinKnapsack(8, {2});
    applyMinKnapsack(9, {});
    applyMinKnapsack(10, {});
}

TEST_F(AMinKnapsack, WeightZero) {
    initialize({1, 1, 0, 0}, {1, 2, 1, 0});
    applyMinKnapsack(0, {0, 1});
    applyMinKnapsack(1, {0});
    applyMinKnapsack(2, {});
}

TEST_F(AnApproximate, TestMultipleParts) {
    initialize(2, {5, 5}, { {{0, 10}, {1, 9}}, {{0, 2}, {1, 1}}, {{0, 4}, {1, 0}},
                            {{0, 3}, {1, 0}}, {{0, 2}, {1, 0}} }, {3, 1, 3, 2, 1});
    setup({0, 1});
    partition(false);

    ASSERT_EQ(1, phg.separatedPartID(0));
    ASSERT_EQ(1, phg.separatedPartID(1));
    ASSERT_EQ(0, phg.separatedPartID(2));
    ASSERT_EQ(0, phg.separatedPartID(3));
    ASSERT_EQ(1, phg.separatedPartID(4));

    setup({0, 1});
    partition(true);

    ASSERT_EQ(1, phg.separatedPartID(0));
    ASSERT_EQ(1, phg.separatedPartID(1));
    ASSERT_EQ(0, phg.separatedPartID(2));
    ASSERT_EQ(0, phg.separatedPartID(3));
    ASSERT_EQ(1, phg.separatedPartID(4));
}

}  // namespace mt_kahypar
