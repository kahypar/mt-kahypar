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


#include <atomic>
#include <cstdlib>
#include <mt-kahypar/macros.h>

#include "gmock/gmock.h"
#include "tbb/task_group.h"

#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/datastructures/sparse_pin_counts.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

template <class F, class K>
void executeConcurrent(F f1, K f2) {
  std::atomic<int> cnt(0);
  tbb::task_group group;

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f1();
      });

  group.run([&] {
        cnt++;
        while (cnt < 2) { }
        f2();
      });

  group.wait();
}

template<typename PinCounts>
class APinCountDataStructure : public Test {

 public:
  APinCountDataStructure() :
    pin_count() { }

  void initialize(const HyperedgeID num_hyperedges,
                  const PartitionID k,
                  const HypernodeID max_value) {
    pin_count.initialize(num_hyperedges, k, max_value);
  }

  PinCounts pin_count;
};

using PinCountTestTypes =
  ::testing::Types<PinCountInPart, SparsePinCounts>;

TYPED_TEST_CASE(APinCountDataStructure, PinCountTestTypes);

TYPED_TEST(APinCountDataStructure, IsZeroInitialized_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(0, this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart1_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 2, 2);
  ASSERT_EQ(2, this->pin_count.pinCountInPart(4, 2));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart2_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 31, 1);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(4, 31));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart3_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(32, 4, 2);
  this->pin_count.setPinCountInPart(32, 5, 1);
  ASSERT_EQ(2, this->pin_count.pinCountInPart(32, 4));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(32, 5));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart4_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<HypernodeID> expected_pin_count(k, 0);
  for ( PartitionID block = 0; block < k; ++block ) {
    expected_pin_count[block] = rand() % max_value;
    this->pin_count.setPinCountInPart(16, block, expected_pin_count[block]);
  }

  for ( PartitionID block = 0; block < k; ++block ) {
    ASSERT_EQ(expected_pin_count[block], this->pin_count.pinCountInPart(16, block));
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart5_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<std::vector<HypernodeID>> expected_pin_count(
    num_hyperedges, std::vector<HypernodeID>(k, 0));
  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      expected_pin_count[he][block] = rand() % max_value;
      this->pin_count.setPinCountInPart(he, block, expected_pin_count[he][block]);
    }
  }

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(expected_pin_count[he][block], this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart1_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 31);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 31));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart2_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 31);
  this->pin_count.incrementPinCountInPart(5, 30);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 31));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 30));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart3_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 31);
  this->pin_count.incrementPinCountInPart(5, 31);
  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 31));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart1_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 31, 2);
  this->pin_count.decrementPinCountInPart(5, 31);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 31));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart2_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 31, 2);
  this->pin_count.setPinCountInPart(5, 30, 1);
  this->pin_count.decrementPinCountInPart(5, 31);
  this->pin_count.decrementPinCountInPart(5, 30);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 31));
  ASSERT_EQ(0, this->pin_count.pinCountInPart(5, 30));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart3_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 31, 2);
  this->pin_count.decrementPinCountInPart(5, 31);
  this->pin_count.decrementPinCountInPart(5, 31);
  ASSERT_EQ(0, this->pin_count.pinCountInPart(5, 31));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently1_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
  });

  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently2_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
    this->pin_count.decrementPinCountInPart(5, 4);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
    this->pin_count.incrementPinCountInPart(6, 1);
  });

  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(2, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently3_k32_Max2) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 32;
  const HypernodeID max_value = 2;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
    this->pin_count.decrementPinCountInPart(5, 4);
    this->pin_count.setPinCountInPart(7, 19, 2);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
    this->pin_count.incrementPinCountInPart(6, 1);
    this->pin_count.setPinCountInPart(4, 18, 1);
  });

  ASSERT_EQ(1, this->pin_count.pinCountInPart(4, 18));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(2, this->pin_count.pinCountInPart(6, 1));
  ASSERT_EQ(2, this->pin_count.pinCountInPart(7, 19));
}

TYPED_TEST(APinCountDataStructure, IsZeroInitialized_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(0, this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart1_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 2, 8);
  ASSERT_EQ(8, this->pin_count.pinCountInPart(4, 2));
}


TYPED_TEST(APinCountDataStructure, SetsPinCountPart2_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 19, 7);
  ASSERT_EQ(7, this->pin_count.pinCountInPart(4, 19));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart3_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(32, 4, 7);
  this->pin_count.setPinCountInPart(32, 5, 6);
  ASSERT_EQ(7, this->pin_count.pinCountInPart(32, 4));
  ASSERT_EQ(6, this->pin_count.pinCountInPart(32, 5));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart4_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<HypernodeID> expected_pin_count(k, 0);
  for ( PartitionID block = 0; block < k; ++block ) {
    expected_pin_count[block] = rand() % max_value;
    this->pin_count.setPinCountInPart(16, block, expected_pin_count[block]);
  }

  for ( PartitionID block = 0; block < k; ++block ) {
    ASSERT_EQ(expected_pin_count[block], this->pin_count.pinCountInPart(16, block));
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart5_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<std::vector<HypernodeID>> expected_pin_count(
    num_hyperedges, std::vector<HypernodeID>(k, 0));
  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      expected_pin_count[he][block] = rand() % max_value;
      this->pin_count.setPinCountInPart(he, block, expected_pin_count[he][block]);
    }
  }

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(expected_pin_count[he][block], this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart1_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart2_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 18);
  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 19));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 18));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart3_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  ASSERT_EQ(5, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart1_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 5);
  this->pin_count.decrementPinCountInPart(5, 19);
  ASSERT_EQ(4, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart2_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 5);
  this->pin_count.setPinCountInPart(5, 18, 4);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 18);
  ASSERT_EQ(4, this->pin_count.pinCountInPart(5, 19));
  ASSERT_EQ(3, this->pin_count.pinCountInPart(5, 18));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart3_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 5);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently1_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
  });

  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently2_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
    this->pin_count.decrementPinCountInPart(5, 4);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
    this->pin_count.incrementPinCountInPart(6, 1);
  });

  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(2, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently3_k20_Max8) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 20;
  const HypernodeID max_value = 8;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 8);
    this->pin_count.decrementPinCountInPart(5, 4);
    this->pin_count.setPinCountInPart(7, 19, 7);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 6);
    this->pin_count.incrementPinCountInPart(6, 1);
    this->pin_count.setPinCountInPart(4, 18, 4);
  });

  ASSERT_EQ(4, this->pin_count.pinCountInPart(4, 18));
  ASSERT_EQ(7, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(7, this->pin_count.pinCountInPart(6, 1));
  ASSERT_EQ(7, this->pin_count.pinCountInPart(7, 19));
}

TYPED_TEST(APinCountDataStructure, IsZeroInitialized_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(0, this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart1_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 2, 30);
  ASSERT_EQ(30, this->pin_count.pinCountInPart(4, 2));
}


TYPED_TEST(APinCountDataStructure, SetsPinCountPart2_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(4, 19, 29);
  ASSERT_EQ(29, this->pin_count.pinCountInPart(4, 19));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart3_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(32, 4, 23);
  this->pin_count.setPinCountInPart(32, 5, 22);
  ASSERT_EQ(23, this->pin_count.pinCountInPart(32, 4));
  ASSERT_EQ(22, this->pin_count.pinCountInPart(32, 5));
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart4_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<HypernodeID> expected_pin_count(k, 0);
  for ( PartitionID block = 0; block < k; ++block ) {
    expected_pin_count[block] = rand() % max_value;
    this->pin_count.setPinCountInPart(16, block, expected_pin_count[block]);
  }

  for ( PartitionID block = 0; block < k; ++block ) {
    ASSERT_EQ(expected_pin_count[block], this->pin_count.pinCountInPart(16, block));
  }
}

TYPED_TEST(APinCountDataStructure, SetsPinCountPart5_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  std::vector<std::vector<HypernodeID>> expected_pin_count(
    num_hyperedges, std::vector<HypernodeID>(k, 0));
  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      expected_pin_count[he][block] = rand() % max_value;
      this->pin_count.setPinCountInPart(he, block, expected_pin_count[he][block]);
    }
  }

  for ( HyperedgeID he = 0; he < num_hyperedges; ++he ) {
    for ( PartitionID block = 0; block < k; ++block ) {
      ASSERT_EQ(expected_pin_count[he][block], this->pin_count.pinCountInPart(he, block));
    }
  }
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart1_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart2_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 18);
  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 19));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 18));
}

TYPED_TEST(APinCountDataStructure, IncrementsPinCountInPart3_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  this->pin_count.incrementPinCountInPart(5, 19);
  ASSERT_EQ(9, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart1_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 20);
  this->pin_count.decrementPinCountInPart(5, 19);
  ASSERT_EQ(19, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart2_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 20);
  this->pin_count.setPinCountInPart(5, 18, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 18);
  ASSERT_EQ(19, this->pin_count.pinCountInPart(5, 19));
  ASSERT_EQ(18, this->pin_count.pinCountInPart(5, 18));
}

TYPED_TEST(APinCountDataStructure, DecrementsPinCountInPart3_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  this->pin_count.setPinCountInPart(5, 19, 20);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  this->pin_count.decrementPinCountInPart(5, 19);
  ASSERT_EQ(10, this->pin_count.pinCountInPart(5, 19));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently1_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
  });

  ASSERT_EQ(2, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(1, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently2_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 2);
    this->pin_count.decrementPinCountInPart(5, 4);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 1);
    this->pin_count.incrementPinCountInPart(6, 1);
  });

  ASSERT_EQ(1, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(2, this->pin_count.pinCountInPart(6, 1));
}

TYPED_TEST(APinCountDataStructure, ModifyTwoHyperedgesConcurrently3_k30_Max30) {
  const HyperedgeID num_hyperedges = 100;
  const PartitionID k = 30;
  const HypernodeID max_value = 30;
  this->initialize(num_hyperedges, k, max_value);

  executeConcurrent([&] {
    this->pin_count.setPinCountInPart(5, 4, 20);
    this->pin_count.decrementPinCountInPart(5, 4);
    this->pin_count.setPinCountInPart(7, 19, 30);
  }, [&] {
    this->pin_count.setPinCountInPart(6, 1, 26);
    this->pin_count.incrementPinCountInPart(6, 1);
    this->pin_count.setPinCountInPart(4, 18, 25);
  });

  ASSERT_EQ(25, this->pin_count.pinCountInPart(4, 18));
  ASSERT_EQ(19, this->pin_count.pinCountInPart(5, 4));
  ASSERT_EQ(27, this->pin_count.pinCountInPart(6, 1));
  ASSERT_EQ(30, this->pin_count.pinCountInPart(7, 19));
}

}  // namespace ds
}  // namespace mt_kahypar
