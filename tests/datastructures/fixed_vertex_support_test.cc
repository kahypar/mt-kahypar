/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {


class AFixedVertexSupport : public Test {

 public:
  using Hypergraph = StaticHypergraph;
  using Factory = typename Hypergraph::Factory;

  AFixedVertexSupport() :
    hypergraph(Factory::construct(
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    fixed_vertices(7, 3) {
    fixed_vertices.setHypergraph(&hypergraph);
    fixed_vertices.fixToBlock(0, 0);
    fixed_vertices.fixToBlock(2, 0);
    fixed_vertices.fixToBlock(4, 1);
    fixed_vertices.fixToBlock(6, 2);
  }

  Hypergraph hypergraph;
  FixedVertexSupport<Hypergraph> fixed_vertices;
};

TEST_F(AFixedVertexSupport, CheckIfNodesAreFixed) {
  ASSERT_TRUE(fixed_vertices.hasFixedVertices());
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_FALSE(fixed_vertices.isFixed(1));
  ASSERT_TRUE(fixed_vertices.isFixed(2));
  ASSERT_FALSE(fixed_vertices.isFixed(3));
  ASSERT_TRUE(fixed_vertices.isFixed(4));
  ASSERT_FALSE(fixed_vertices.isFixed(5));
  ASSERT_TRUE(fixed_vertices.isFixed(6));
}

TEST_F(AFixedVertexSupport, CheckFixedVertexBlocks) {
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_EQ(kInvalidPartition, fixed_vertices.fixedVertexBlock(1));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(2));
  ASSERT_EQ(kInvalidPartition, fixed_vertices.fixedVertexBlock(3));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlock(4));
  ASSERT_EQ(kInvalidPartition, fixed_vertices.fixedVertexBlock(5));
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlock(6));
}

TEST_F(AFixedVertexSupport, CheckFixedVertexBlockWeights) {
  ASSERT_EQ(4, fixed_vertices.totalFixedVertexWeight());
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlockWeight(0));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlockWeight(1));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlockWeight(2));
}

TEST_F(AFixedVertexSupport, ContractFreeOntoFreeVertex) {
  ASSERT_TRUE(fixed_vertices.contract(3, 5));
  ASSERT_FALSE(fixed_vertices.isFixed(3));
  ASSERT_FALSE(fixed_vertices.isFixed(5));
}

TEST_F(AFixedVertexSupport, ContractFreeOntoFixedVertex) {
  ASSERT_TRUE(fixed_vertices.contract(0, 3));
  ASSERT_TRUE(fixed_vertices.isFixed(3));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(3));
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_EQ(3, fixed_vertices.fixedVertexBlockWeight(0));
  ASSERT_EQ(5, fixed_vertices.totalFixedVertexWeight());
}

TEST_F(AFixedVertexSupport, ContractFixedOntoFixedVertex1) {
  ASSERT_TRUE(fixed_vertices.contract(0, 2));
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_TRUE(fixed_vertices.isFixed(2));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(2));
}

TEST_F(AFixedVertexSupport, ContractFixedOntoFixedVertex2) {
  ASSERT_FALSE(fixed_vertices.contract(0, 6));
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_TRUE(fixed_vertices.isFixed(6));
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlock(6));
}

TEST_F(AFixedVertexSupport, ContractFixedOntoFreeVertex) {
  ASSERT_TRUE(fixed_vertices.contract(1, 4));
  ASSERT_TRUE(fixed_vertices.isFixed(1));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlock(1));
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlockWeight(1));
  ASSERT_EQ(5, fixed_vertices.totalFixedVertexWeight());
}

TEST_F(AFixedVertexSupport, UnontractFreeOntoFreeVertex) {
  ASSERT_TRUE(fixed_vertices.contract(3, 5));
  fixed_vertices.uncontract(3, 5);
  ASSERT_FALSE(fixed_vertices.isFixed(3));
  ASSERT_FALSE(fixed_vertices.isFixed(5));
}

TEST_F(AFixedVertexSupport, UnontractFreeOntoFixedVertex) {
  ASSERT_TRUE(fixed_vertices.contract(0, 3));
  fixed_vertices.uncontract(0, 3);
  ASSERT_FALSE(fixed_vertices.isFixed(3));
  ASSERT_EQ(kInvalidPartition, fixed_vertices.fixedVertexBlock(3));
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlockWeight(0));
  ASSERT_EQ(4, fixed_vertices.totalFixedVertexWeight());
}

TEST_F(AFixedVertexSupport, UncontractFixedOntoFixedVertex) {
  ASSERT_TRUE(fixed_vertices.contract(0, 2));
  fixed_vertices.uncontract(0, 2);
  ASSERT_TRUE(fixed_vertices.isFixed(0));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(0));
  ASSERT_TRUE(fixed_vertices.isFixed(2));
  ASSERT_EQ(0, fixed_vertices.fixedVertexBlock(2));
  ASSERT_EQ(2, fixed_vertices.fixedVertexBlockWeight(0));
  ASSERT_EQ(4, fixed_vertices.totalFixedVertexWeight());
}

TEST_F(AFixedVertexSupport, UncontractFixedOntoFreeVertex) {
  ASSERT_TRUE(fixed_vertices.contract(1, 4));
  fixed_vertices.uncontract(1, 4);
  ASSERT_FALSE(fixed_vertices.isFixed(1));
  ASSERT_EQ(kInvalidPartition, fixed_vertices.fixedVertexBlock(1));
  ASSERT_TRUE(fixed_vertices.isFixed(4));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlock(4));
  ASSERT_EQ(1, fixed_vertices.fixedVertexBlockWeight(1));
  ASSERT_EQ(4, fixed_vertices.totalFixedVertexWeight());
}


}  // namespace ds
}  // namespace mt_kahypar
