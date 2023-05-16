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

#include "mt-kahypar/datastructures/static_graph_factory.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"

using ::testing::Test;

namespace mt_kahypar {

class AProcessGraph : public Test {
 public:
  AProcessGraph() :
    graph(nullptr) {

    /**
     * Process Graph:
     *        1           2           4
     * 0  -------- 1  -------- 2  -------- 3
     * |           |           |           |
     * | 3         | 2         | 1         | 1
     * |      3    |      2    |      1    |
     * 4  -------- 5  -------- 6  -------- 7
     * |           |           |           |
     * | 1         | 1         | 3         | 2
     * |      2    |      4    |      2    |
     * 8  -------- 9  -------- 10 -------- 11
     * |           |           |           |
     * | 1         | 2         | 2         | 2
     * |      1    |      1    |      2    |
     * 12 -------- 13 -------- 14 -------- 15
    */
    vec<HyperedgeWeight> edge_weights =
      { 1, 2, 4,
        3, 2, 1, 1,
        3, 2, 1,
        1, 1, 3, 2,
        2, 4, 2,
        1, 2, 2, 2,
        1, 1, 2 };
    graph = std::make_unique<ProcessGraph>(
      ds::StaticGraphFactory::construct(16, 24,
        { { 0, 1 }, { 1, 2 }, { 2, 3 },
          { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
          { 4, 5 }, { 5, 6 }, { 6, 7 },
          { 4, 8 }, { 5, 9 }, { 6, 10 }, { 7, 11 },
          { 8, 9 }, { 9, 10 }, { 10, 11 },
          { 8, 12 }, { 9, 13 }, { 10, 14 }, { 11, 15 },
          { 12, 13 }, { 13, 14 }, { 14, 15 } },
          edge_weights.data()));

  }

  std::unique_ptr<ProcessGraph> graph;
};

TEST_F(AProcessGraph, HasCorrectNumberOfBlocks) {
  ASSERT_EQ(16, graph->numBlocks());
}

TEST_F(AProcessGraph, ComputesAllShortestPaths) {
  graph->precomputeDistances(2);
  ASSERT_EQ(1, graph->distance(0, 1));
  ASSERT_EQ(3, graph->distance(0, 2));
  ASSERT_EQ(6, graph->distance(0, 3));
  ASSERT_EQ(6, graph->distance(3, 0));
  ASSERT_EQ(6, graph->distance(1, 10));
  ASSERT_EQ(7, graph->distance(8, 11));
  ASSERT_EQ(7, graph->distance(12, 7));
  ASSERT_EQ(3, graph->distance(9, 14));
  ASSERT_EQ(5, graph->distance(6, 15));
  ASSERT_EQ(2, graph->distance(3, 6));
  ASSERT_EQ(7, graph->distance(4, 3));
}


}  // namespace mt_kahypar
