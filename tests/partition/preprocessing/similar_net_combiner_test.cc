/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/sparsification/policies/similiar_net_combine.h"

using ::testing::Test;

namespace mt_kahypar {

namespace {
  using Hyperedge = parallel::scalable_vector<HypernodeID>;
} // namespace

Hypergraph createTestHypergraph() {
  const std::vector<HyperedgeWeight> hyperedge_weight = { 2, 3, 1, 4, 5, 6 };
  return HypergraphFactory::construct(7, 6, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6}, {0, 3, 5}, {0, 1, 2} },
    hyperedge_weight.data());
}

bool isEqual(const Hyperedge& actual, const Hyperedge& expected) {
  if ( actual.size() != expected.size() ) {
    return false;
  }

  for ( size_t i = 0; i < actual.size(); ++i ) {
    if ( actual[i] != expected[i] ) {
      return false;
    }
  }
  return true;
}

bool containsNoDuplicates(const Hyperedge& e) {
  for ( size_t i = 1; i < e.size(); ++i ) {
    if ( e[i - 1] == e[i] ) {
      return false;
    }
  }
  return true;
}

TEST(ASimilarNetCombiner, UnionTwoNets1) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 2, 3 };
  const Hyperedge rhs = { 2, 3, 4, 5 };
  const Hyperedge combined_he = UnionCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(std::is_sorted(combined_he.begin(), combined_he.end()));
  ASSERT_TRUE(containsNoDuplicates(combined_he));
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3, 4, 5 }));
}

TEST(ASimilarNetCombiner, UnionTwoNets2) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 2, 3 };
  const Hyperedge rhs = { 0, 1, 2, 3 };
  const Hyperedge combined_he = UnionCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(std::is_sorted(combined_he.begin(), combined_he.end()));
  ASSERT_TRUE(containsNoDuplicates(combined_he));
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3 }));
}

TEST(ASimilarNetCombiner, UnionTwoNets3) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 2, 3 };
  const Hyperedge rhs = { 4, 5, 6, 7 };
  const Hyperedge combined_he = UnionCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(std::is_sorted(combined_he.begin(), combined_he.end()));
  ASSERT_TRUE(containsNoDuplicates(combined_he));
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3, 4, 5, 6, 7 }));
}

TEST(ASimilarNetCombiner, UnionTwoNets4) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 4, 5, 6, 7 };
  const Hyperedge rhs = { 0, 1, 2, 3 };
  const Hyperedge combined_he = UnionCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(std::is_sorted(combined_he.begin(), combined_he.end()));
  ASSERT_TRUE(containsNoDuplicates(combined_he));
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3, 4, 5, 6, 7 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithMaxSizeStrategy1) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 2, 3 };
  const Hyperedge rhs = { 2, 3, 4, 5 };
  const Hyperedge combined_he = MaxSizeCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithMaxSizeStrategy2) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1 };
  const Hyperedge rhs = { 2, 3, 4, 5 };
  const Hyperedge combined_he = MaxSizeCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 2, 3, 4, 5 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithMaxSizeStrategy3) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 2, 3 };
  const Hyperedge rhs = { 2, 3 };
  const Hyperedge combined_he = MaxSizeCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2, 3 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithNetImportanceStrategy1) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 1, 3, 4 };
  const Hyperedge rhs = { 2, 5, 6 };
  const Hyperedge combined_he = NetImportanceCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 3, 4 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithNetImportanceStrategy2) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 3, 4, 6 };
  const Hyperedge rhs = { 0, 2 };
  const Hyperedge combined_he = NetImportanceCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 0, 2 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithNetImportanceStrategy3) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 0, 3, 5 };
  const Hyperedge rhs = { 0, 1, 2 };
  const Hyperedge combined_he = NetImportanceCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 0, 1, 2 }));
}

TEST(ASimilarNetCombiner, UnifiesTwoNetsWithNetImportanceStrategy4) {
  auto hypergraph = createTestHypergraph();
  const Hyperedge lhs = { 3, 4, 6 };
  const Hyperedge rhs = { 2, 5, 6 };
  const Hyperedge combined_he = NetImportanceCombiner::combine(hypergraph, lhs, rhs);
  ASSERT_TRUE(isEqual(combined_he, { 2, 5, 6 }));
}

}  // namespace mt_kahypar
