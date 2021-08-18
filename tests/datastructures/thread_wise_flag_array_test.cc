//
// Created by mlaupichler on 25.07.21.
//

#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/async/thread_wise_flag_array.h"

namespace mt_kahypar::ds {

    using IndexType = uint32_t;

  TEST(AThreadWiseFlagArray, SequentialOperations) {

      std::vector<size_t> thread_numbers = {8, 60, 128};
      const IndexType num_elements = 10;

      for (auto num_threads : thread_numbers) {

        auto array = ThreadWiseFlagArray<IndexType>(num_elements, num_threads);

        // Expect all false and setting to false to return false after initialization
        for (IndexType i = 0; i < num_elements; ++i) {
          ASSERT_FALSE(array.any_set(i));
          for (size_t t = 0; t < num_threads; ++t) {
            ASSERT_FALSE(array.is_set(i, t));
            ASSERT_FALSE(array.any_set_except_thread(i, t));
            bool set_block_to_zero = false;
            bool change = array.set_false(i, t, set_block_to_zero);
            ASSERT_FALSE(change);
            ASSERT_FALSE(set_block_to_zero);
          }
        }

        // Set to true and test is_set
        bool block_was_zero = false;
        bool change = array.set_true(0, 0, block_was_zero);
        ASSERT_TRUE(change);
        ASSERT_TRUE(block_was_zero);
        ASSERT(array.is_set(0, 0));

        // Make sure any_set is true but any_set_except_thread for the one just set is false still
        ASSERT_TRUE(array.any_set(0));
        ASSERT_FALSE(array.any_set_except_thread(0, 0));

        // Set same entry to true and assert return value is false and block_was_zero is false
        block_was_zero = false;
        change = array.set_true(0, 0, block_was_zero);
        ASSERT_FALSE(change);
        ASSERT_FALSE(block_was_zero);
        ASSERT(array.is_set(0, 0));

        // Change entry for other thread but same element and make sure block_was_zero is false
        block_was_zero = false;
        change = array.set_true(0, 1, block_was_zero);
        ASSERT_TRUE(change);
        ASSERT_FALSE(block_was_zero);

        // Make sure any_set is true and any_set_except_thread for each one just set is true too
        ASSERT_TRUE(array.any_set(0));
        ASSERT_TRUE(array.any_set_except_thread(0, 0));
        ASSERT_TRUE(array.any_set_except_thread(0, 1));

        // Set entry for one of the threads back to false and make sure set_block_to_zero is false
        bool set_block_to_zero = false;
        change = array.set_false(0, 1, set_block_to_zero);
        ASSERT_TRUE(change);
        ASSERT_FALSE(set_block_to_zero);
        ASSERT_FALSE(array.is_set(0, 1));
        ASSERT_TRUE(array.is_set(0, 0));

        // Set last entry for element to false and make sure set_block_to_zero is true
        set_block_to_zero = false;
        change = array.set_false(0, 0, set_block_to_zero);
        ASSERT_TRUE(change);
        ASSERT_TRUE(set_block_to_zero);
        ASSERT_FALSE(array.is_set(0, 0));

        ASSERT_FALSE(array.any_set(0));
      }
    }


} // end namespace mt_kahypar::ds
