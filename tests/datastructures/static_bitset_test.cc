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

#include <numeric>
#include <algorithm>

#include "gmock/gmock.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/static_bitset.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using Block = StaticBitset::Block;

void set_one_bits(Block& block,
                  const vec<size_t>& one_bits) {
  for ( const size_t pos : one_bits ) {
    ASSERT(pos < 64);
    block |= ( UL(1) << pos );
  }
}

void verify_iterator(const StaticBitset& bitset,
                     const vec<PartitionID>& expected) {
  ASSERT_EQ(static_cast<size_t>(bitset.popcount()), expected.size());
  size_t cnt = 0;
  for ( const PartitionID& block : bitset ) {
    ASSERT_EQ(block, expected[cnt++]);
  }
  ASSERT_EQ(cnt, expected.size());
}

TEST(AStaticBitset, CountNumberOfOneBits1) {
  Block block = 0;
  set_one_bits(block, { 0, 1, 5, 7 });
  StaticBitset bitset(1, &block);
  ASSERT_EQ(4, bitset.popcount());
}

TEST(AStaticBitset, CountNumberOfOneBits2) {
  Block block = 0;
  set_one_bits(block, { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  StaticBitset bitset(1, &block);
  ASSERT_EQ(10, bitset.popcount());
}

TEST(AStaticBitset, CountNumberOfOneBits3) {
  vec<Block> blocks(2, 0);
  set_one_bits(blocks[0], { 0, 1, 5, 7 });
  set_one_bits(blocks[1], { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  StaticBitset bitset(2, blocks.data());
  ASSERT_EQ(14, bitset.popcount());
}

TEST(AStaticBitset, CountNumberOfOneBits4) {
  vec<Block> blocks(3, 0);
  set_one_bits(blocks[0], { 0, 1, 5, 7 });
  set_one_bits(blocks[1], { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  set_one_bits(blocks[2], { 23, 24, 25, 26, 27, 55 });
  StaticBitset bitset(3, blocks.data());
  ASSERT_EQ(20, bitset.popcount());
}

TEST(AStaticBitset, VerifyIterator1) {
  Block block = 0;
  set_one_bits(block, { 0, 1, 5, 7 });
  StaticBitset bitset(1, &block);
  verify_iterator(bitset, {0, 1, 5, 7});
}

TEST(AStaticBitset, VerifyIterator2) {
  Block block = 0;
  set_one_bits(block, { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  StaticBitset bitset(1, &block);
  verify_iterator(bitset, { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
}

TEST(AStaticBitset, VerifyIterator3) {
  vec<Block> blocks(2, 0);
  set_one_bits(blocks[0], { 0, 1, 5, 7 });
  set_one_bits(blocks[1], { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  StaticBitset bitset(2, blocks.data());
  verify_iterator(bitset, {0, 1, 5, 7, 64, 65, 69, 71, 80, 88, 94, 95, 110, 127});
}

TEST(AStaticBitset, VerifyIterator4) {
  vec<Block> blocks(3, 0);
  set_one_bits(blocks[0], { 0, 1, 5, 7 });
  set_one_bits(blocks[1], { 0, 1, 5, 7, 16, 24, 30, 31, 46, 63 });
  set_one_bits(blocks[2], { 23, 24, 25, 26, 27, 55 });
  StaticBitset bitset(3, blocks.data());
  verify_iterator(bitset, {0, 1, 5, 7, 64, 65, 69, 71, 80, 88,
                           94, 95, 110, 127, 151, 152, 153, 154, 155,  183});
}

}  // namespace ds
}  // namespace mt_kahypar
