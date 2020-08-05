/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/contraction_tree.h"

namespace mt_kahypar {
namespace ds {

void verifyChilds(const ContractionTree& tree,
                  const HypernodeID u,
                  const std::set<HypernodeID>& expected_childs) {
  size_t num_childs = 0;
  for ( const HypernodeID v : tree.childs(u) ) {
    ASSERT_TRUE(expected_childs.find(v) != expected_childs.end())
      << "Child " << v << " not found!";
    ++num_childs;
  }
  ASSERT_EQ(num_childs, expected_childs.size());
  ASSERT_EQ(num_childs, tree.degree(u));
}

void verifyRoots(const parallel::scalable_vector<HypernodeID>& actual_roots,
                 const std::set<HypernodeID> expected_roots) {
  ASSERT_EQ(expected_roots.size(), actual_roots.size());
  for ( const HypernodeID& u : actual_roots ) {
    ASSERT_TRUE(expected_roots.find(u) != expected_roots.end())
      << "Root " << u << " not found!";
  }
}

TEST(AContractionTree, IsConstructedCorrectly1) {
  ContractionTree tree;
  tree.initialize(5);
  tree.setParent(1, 0);
  tree.setParent(2, 0);
  tree.setParent(3, 1);
  tree.setParent(4, 1);
  tree.finalize();

  verifyChilds(tree, 0, { 1, 2 });
  verifyChilds(tree, 1, { 3, 4 });
  verifyChilds(tree, 2, { });
  verifyChilds(tree, 3, { });
  verifyChilds(tree, 4, { });

  ASSERT_EQ(4, tree.subtreeSize(0));
  ASSERT_EQ(2, tree.subtreeSize(1));
  ASSERT_EQ(0, tree.subtreeSize(2));
  ASSERT_EQ(0, tree.subtreeSize(3));
  ASSERT_EQ(0, tree.subtreeSize(4));

  verifyRoots(tree.roots(), { 0 });
}

TEST(AContractionTree, IsConstructedCorrectly2) {
  ContractionTree tree;
  tree.initialize(10);
  tree.setParent(1, 0);
  tree.setParent(2, 0);
  tree.setParent(3, 1);
  tree.setParent(4, 1);
  tree.setParent(6, 5);
  tree.setParent(7, 5);
  tree.setParent(8, 5);
  tree.setParent(9, 5);
  tree.finalize();

  verifyChilds(tree, 0, { 1, 2 });
  verifyChilds(tree, 1, { 3, 4 });
  verifyChilds(tree, 2, { });
  verifyChilds(tree, 3, { });
  verifyChilds(tree, 4, { });
  verifyChilds(tree, 5, { 6, 7, 8, 9 });
  verifyChilds(tree, 6, { });
  verifyChilds(tree, 7, { });
  verifyChilds(tree, 8, { });
  verifyChilds(tree, 9, { });

  ASSERT_EQ(4, tree.subtreeSize(0));
  ASSERT_EQ(2, tree.subtreeSize(1));
  ASSERT_EQ(0, tree.subtreeSize(2));
  ASSERT_EQ(0, tree.subtreeSize(3));
  ASSERT_EQ(0, tree.subtreeSize(4));
  ASSERT_EQ(4, tree.subtreeSize(5));
  ASSERT_EQ(0, tree.subtreeSize(6));
  ASSERT_EQ(0, tree.subtreeSize(7));
  ASSERT_EQ(0, tree.subtreeSize(8));
  ASSERT_EQ(0, tree.subtreeSize(9));

  verifyRoots(tree.roots(), { 0, 5 });
}

TEST(AContractionTree, IsConstructedCorrectly3) {
  ContractionTree tree;
  tree.initialize(10);
  tree.setParent(1, 0);
  tree.setParent(2, 0);
  tree.setParent(3, 1);
  tree.setParent(4, 1);
  tree.setParent(5, 2);
  tree.setParent(6, 2);
  tree.setParent(7, 3);
  tree.setParent(8, 3);
  tree.setParent(9, 4);
  tree.finalize();

  verifyChilds(tree, 0, { 1, 2 });
  verifyChilds(tree, 1, { 3, 4 });
  verifyChilds(tree, 2, { 5, 6 });
  verifyChilds(tree, 3, { 7, 8 });
  verifyChilds(tree, 4, { 9 });
  verifyChilds(tree, 5, { });
  verifyChilds(tree, 6, { });
  verifyChilds(tree, 7, { });
  verifyChilds(tree, 8, { });
  verifyChilds(tree, 9, { });

  ASSERT_EQ(9, tree.subtreeSize(0));
  ASSERT_EQ(5, tree.subtreeSize(1));
  ASSERT_EQ(2, tree.subtreeSize(2));
  ASSERT_EQ(2, tree.subtreeSize(3));
  ASSERT_EQ(1, tree.subtreeSize(4));
  ASSERT_EQ(0, tree.subtreeSize(5));
  ASSERT_EQ(0, tree.subtreeSize(6));
  ASSERT_EQ(0, tree.subtreeSize(7));
  ASSERT_EQ(0, tree.subtreeSize(8));
  ASSERT_EQ(0, tree.subtreeSize(9));

  verifyRoots(tree.roots(), { 0 });
}

TEST(AContractionTree, IsConstructedCorrectly4) {
  ContractionTree tree;
  tree.initialize(21);
  tree.setParent(1, 0);
  tree.setParent(2, 0);
  tree.setParent(4, 3);
  tree.setParent(5, 3);
  tree.setParent(7, 6);
  tree.setParent(8, 6);
  tree.setParent(10, 9);
  tree.setParent(11, 9);
  tree.setParent(13, 12);
  tree.setParent(14, 12);
  tree.setParent(16, 15);
  tree.setParent(17, 15);
  tree.setParent(19, 18);
  tree.setParent(20, 18);
  tree.finalize();

  verifyChilds(tree, 0, { 1, 2 });
  verifyChilds(tree, 1, { });
  verifyChilds(tree, 2, { });
  verifyChilds(tree, 3, { 4, 5 });
  verifyChilds(tree, 4, { });
  verifyChilds(tree, 5, { });
  verifyChilds(tree, 6, { 7, 8 });
  verifyChilds(tree, 7, { });
  verifyChilds(tree, 8, { });
  verifyChilds(tree, 9, { 10, 11 });
  verifyChilds(tree, 10, { });
  verifyChilds(tree, 11, { });
  verifyChilds(tree, 12, { 13, 14 });
  verifyChilds(tree, 13, { });
  verifyChilds(tree, 14, { });
  verifyChilds(tree, 15, { 16, 17 });
  verifyChilds(tree, 16, { });
  verifyChilds(tree, 17, { });
  verifyChilds(tree, 18, { 19, 20 });
  verifyChilds(tree, 19, { });
  verifyChilds(tree, 20, { });

  ASSERT_EQ(2, tree.subtreeSize(0));
  ASSERT_EQ(0, tree.subtreeSize(1));
  ASSERT_EQ(0, tree.subtreeSize(2));
  ASSERT_EQ(2, tree.subtreeSize(3));
  ASSERT_EQ(0, tree.subtreeSize(4));
  ASSERT_EQ(0, tree.subtreeSize(5));
  ASSERT_EQ(2, tree.subtreeSize(6));
  ASSERT_EQ(0, tree.subtreeSize(7));
  ASSERT_EQ(0, tree.subtreeSize(8));
  ASSERT_EQ(2, tree.subtreeSize(9));
  ASSERT_EQ(0, tree.subtreeSize(10));
  ASSERT_EQ(0, tree.subtreeSize(11));
  ASSERT_EQ(2, tree.subtreeSize(12));
  ASSERT_EQ(0, tree.subtreeSize(13));
  ASSERT_EQ(0, tree.subtreeSize(14));
  ASSERT_EQ(2, tree.subtreeSize(15));
  ASSERT_EQ(0, tree.subtreeSize(16));
  ASSERT_EQ(0, tree.subtreeSize(17));
  ASSERT_EQ(2, tree.subtreeSize(18));
  ASSERT_EQ(0, tree.subtreeSize(19));
  ASSERT_EQ(0, tree.subtreeSize(20));

  verifyRoots(tree.roots(), { 0, 3, 6, 9, 12, 15, 18 });
}

TEST(AContractionTree, ContainsCorrectRootsInPresenceOfSingletonRoots) {
  ContractionTree tree;
  tree.initialize(10);
  tree.setParent(1, 0);
  tree.setParent(2, 0);
  tree.finalize();

  verifyRoots(tree.roots(), { 0 });
}

} // namespace ds
} // namespace mt_kahypar