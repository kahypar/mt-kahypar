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

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/spanning_tree.h"

using ::testing::Test;

namespace mt_kahypar {
using star_partitioning::SpanningTree;

class ASpanningTree : public Test {
 public:
  ASpanningTree() { }

  void assertChildrenEqual(SpanningTree& tree, HypernodeID node,
                           const std::vector<HypernodeID>& expected) {
    auto it = tree.children(node).begin();
    auto end = tree.children(node).end();
    size_t i = 0;
    for (; it < end; ++i) {
      ASSERT_EQ(expected[i], *it);
      ++it;
    }
    ASSERT_EQ(expected.size(), i);
  }
};

TEST_F(ASpanningTree, addsChildrenAndTracesDepth1) {
  SpanningTree tree(10, 0, 2);
  ASSERT_EQ(0, tree.depth(0));
  ASSERT_EQ(1, tree.requiredDepthForChild(0));
  ASSERT_TRUE(tree.mayAddChild(0));
  assertChildrenEqual(tree, 0, {});

  tree.addChild(0, 1);
  tree.addChild(0, 2);
  ASSERT_EQ(1, tree.depth(1));
  ASSERT_EQ(2, tree.requiredDepthForChild(1));
  ASSERT_FALSE(tree.mayAddChild(1));
  ASSERT_EQ(1, tree.depth(2));
  ASSERT_EQ(2, tree.requiredDepthForChild(2));
  ASSERT_FALSE(tree.mayAddChild(2));
  assertChildrenEqual(tree, 0, {1, 2});
}

TEST_F(ASpanningTree, addsChildrenAndTracesDepth2) {
  const HypernodeID num_nodes = 4 * SpanningTree::max_children + 10;
  std::vector<HypernodeID> first_pseudo;
  std::vector<HypernodeID> second_pseudo{ 1 };

  SpanningTree tree(num_nodes, 0, 3);
  HypernodeID current_node = 1;
  for (HypernodeID i = 1; i < SpanningTree::max_children; ++i) {
    ASSERT_EQ(1, tree.requiredDepthForChild(0));
    tree.addChild(0, current_node++);
  }
  ASSERT_EQ(2, tree.requiredDepthForChild(0));
  tree.addChild(0, current_node++);
  assertChildrenEqual(tree, num_nodes, {current_node - 1});
  first_pseudo.push_back(current_node - 1);

  // make nodes > 1 to non-leafs
  ASSERT_EQ(2, tree.requiredDepthForChild(0));
  for (HypernodeID i = 2; i < SpanningTree::max_children; ++i) {
    ASSERT_EQ(2, tree.requiredDepthForChild(i));
    tree.addChild(i, current_node++);
  }
  ASSERT_TRUE(tree.isLeaf(1));
  ASSERT_FALSE(tree.isLeaf(2));

  // fill the two pseudo nodes
  for (HypernodeID i = 2; i < SpanningTree::max_children; ++i) {
    tree.addChild(0, current_node++);
    first_pseudo.push_back(current_node - 1);
  }
  ASSERT_EQ(2, tree.requiredDepthForChild(0));
  ASSERT_TRUE(tree.mayAddChild(0));
  for (HypernodeID i = 2; i < SpanningTree::max_children; ++i) {
    tree.addChild(0, current_node++);
    second_pseudo.push_back(current_node - 1);
  }
  std::swap(second_pseudo[0], second_pseudo[1]);
  ASSERT_EQ(3, tree.requiredDepthForChild(0));
  ASSERT_FALSE(tree.mayAddChild(0));

  assertChildrenEqual(tree, num_nodes, first_pseudo);
  assertChildrenEqual(tree, num_nodes + 1, second_pseudo);
}

TEST_F(ASpanningTree, ConstructsSpanningTree) {
  std::vector<HyperedgeWeight> edge_weights{4, 3, 2, 4, 1};
  Hypergraph hg = HypergraphFactory::construct(4, 5, { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3} }, edge_weights.data());
  SpanningTree tree = star_partitioning::constructMaxSpanningTree(hg, 3);
  assertChildrenEqual(tree, 0, { 1, 3 });
  assertChildrenEqual(tree, 1, { 2 });
  tree = star_partitioning::constructMaxSpanningTree(hg, 2);
  assertChildrenEqual(tree, 0, { 1, 2, 3 });
}

}  // namespace mt_kahypar
