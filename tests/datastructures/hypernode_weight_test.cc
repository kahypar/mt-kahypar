/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "tests/definitions.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"
#include "mt-kahypar/weight/hypernode_weight_buffer.h"

using ::testing::Test;

namespace mt_kahypar {

using weight::HNWeightScalar;

struct AHypernodeWeight : public Test {
};

TEST_F(AHypernodeWeight, BasicArithmetic) {
  HNWeightScalar vec_1[2] = {1, 2};
  HNWeightScalar vec_2[2] = {3, 4};
  HNWeightScalar vec_3[2] = {-1, 7};
  HNWeightScalar vec_sum[2];

  HNWeightConstRef left(vec_1, 2);
  HNWeightConstRef right(vec_2, 2);
  auto sum = left + right;

  HNWeightConstRef copy(left);
  ASSERT_EQ(copy, left);
  HNWeightConstRef assign = left;
  ASSERT_EQ(assign, left);

  HNWeightRef tmp_mut(vec_3, 2);
  HNWeightRef third_mut(std::move(tmp_mut));
  HNWeightRef sum_mut(vec_sum, 2);
  sum_mut = sum;

  HNWeightScalar buffer[2];
  HNWeightRef diff(buffer, 2);
  diff = sum_mut - third_mut;
  diff = 3 * diff;

  HNWeightScalar vec_expected[2] = {15, -3};
  HNWeightConstRef expected(vec_expected, 2);
  ASSERT_EQ(diff.at(0), 15);
  ASSERT_EQ(diff.at(1), -3);
  ASSERT_EQ(diff, expected);
  ASSERT_EQ(3 * (left + right) - 2 * third_mut, diff + third_mut);

  HNWeightAtomicCRef atomic_ref1 = left.load(std::memory_order_relaxed);
  ASSERT_EQ(atomic_ref1.at(0), 1);
  ASSERT_EQ(atomic_ref1.at(1), 2);
  HNWeightAtomicCRef atomic_ref2 = std::move(diff).load(std::memory_order_relaxed);
  ASSERT_EQ(atomic_ref2.at(0), 15);
  ASSERT_EQ(atomic_ref2.at(1), -3);
  auto atomic_sum = atomic_ref1 + atomic_ref2;
  ASSERT_EQ(atomic_sum.at(0), 16);
  ASSERT_EQ(atomic_sum.at(1), -1);
}

TEST_F(AHypernodeWeight, SimpleWeightArray) {
  HNWeightScalar chunk_1[3] = {1, 2, 3};
  HNWeightConstRef vec_1(chunk_1, 3);
  HNWeightScalar chunk_2[3] = {4, 5, 6};
  HNWeightConstRef vec_2(chunk_2, 3);
  HNWeightScalar chunk_3[3] = {7, 8, 9};
  HNWeightConstRef vec_3(chunk_3, 3);

  HypernodeWeightArray array(3, 3, 0);
  ASSERT_TRUE(weight::isZero(array[0]));
  ASSERT_TRUE(weight::isZero(array[1]));
  ASSERT_TRUE(weight::isZero(array[2]));

  array[0] = vec_1;
  array[1] = vec_2;
  array[2] = vec_3;
  ASSERT_FALSE(weight::isZero(array[0]));
  ASSERT_FALSE(weight::isZero(array[1]));
  ASSERT_FALSE(weight::isZero(array[2]));
  ASSERT_EQ(array[0], vec_1);
  ASSERT_EQ(array[1], vec_2);
  ASSERT_EQ(array[2], vec_3);
}

TEST_F(AHypernodeWeight, SimpleWeightBuffer) {
  weight::HypernodeWeightBuffer buffer(2);

  HNWeightRef first = weight::allocateAt(buffer, 0);
  first.set(0, 1);
  first.set(1, 2);

  HNWeightRef second = weight::allocateAt(buffer, 3);
  second.set(0, 4);
  second.set(1, 6);

  HNWeightRef third = weight::allocateAt(buffer, 42);
  third.set(0, 1);
  third.set(1, 2);

  ASSERT_NE(first, second);
  ASSERT_EQ(first, third);
}

TEST_F(AHypernodeWeight, SimpleAllocatedWeight) {
  AllocatedHNWeight weight_1(3, 1);
  AllocatedHNWeight weight_2(3, 2);

  HNWeightScalar chunk_1[3] = {1, 2, 3};
  HNWeightConstRef vec_1(chunk_1, 3);

  weight_1 = weight_2;
  ASSERT_EQ(weight_1, weight_2);
  weight_1 = vec_1;
  ASSERT_EQ(weight_1, vec_1);
}

}  // namespace mt_kahypar
