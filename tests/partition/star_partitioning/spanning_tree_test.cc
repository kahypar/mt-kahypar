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

}  // namespace mt_kahypar
