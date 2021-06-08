//
// Created by mlaupichler on 08.06.21.
//

#include "gtest/gtest.h"
#include "mt-kahypar/datastructures/async/iterable_bit_set.h"


using ::testing::Test;

namespace mt_kahypar::ds {

    using IBS = IterableBitSet<PartitionID>;
    static constexpr PartitionID N = 10;

    void verifyBitSet(const IBS& bitset, const PartitionID n, const std::set<PartitionID>& expected_set) {
        for (PartitionID i = 0; i < n; ++i) {
            if (expected_set.find(i) != expected_set.end()) {
                ASSERT_TRUE(bitset.isSet(i)) << V(i);
            } else {
                ASSERT_FALSE(bitset.isSet(i)) << V(i);
            }
        }

        // Verify iterator
        size_t num = 0;
        for (const PartitionID& id : bitset) {
            ASSERT_TRUE(expected_set.find(id) != expected_set.end()) << V(id);
            ++num;
        }
        ASSERT_EQ(expected_set.size(), num);
    }

    TEST(AIterableBitSet, SetAndResetFunctions) {
        auto bitset = IBS(N);
        verifyBitSet(bitset, N, {});

        bitset.set_true(0); bitset.set_true(1); bitset.set_true(N-1);
        verifyBitSet(bitset, N, {0, 1, N-1});

        bitset.set_false(0); bitset.set_false(N-1);
        verifyBitSet(bitset, N, {1});
    }

    TEST(AIterableBitSet, RAOperatorBitReference) {
        auto bitset = IBS(N);
        bool bit_val = bitset.isSet(0);
        bool conv_from_bit_ref = bitset[0];
        ASSERT_EQ(bit_val, conv_from_bit_ref);

        bitset[0] = true;
        ASSERT_TRUE(bitset.isSet(0));
        bitset[0] = false;
        ASSERT_FALSE(bitset.isSet(0));

        bitset.set_true(1);
        auto bit_ref = bitset[1];
        ASSERT_TRUE(bit_ref);
        bitset.set_false(1);
        ASSERT_FALSE(bit_ref);

        verifyBitSet(bitset, N, {});
        bitset[0] = true; bitset[1] = true; bitset[N-1] = true;
        verifyBitSet(bitset, N, {0, 1, N-1});
    }

    TEST(AIterableBitSet, EmptyIteration) {
        auto bitset = IBS(N);

        auto begin = bitset.begin();
        auto end = bitset.end();
        ASSERT_EQ(begin, end);
        ASSERT_EQ(*begin, N);
        ASSERT_EQ(*end, N);

        HyperedgeID num = 0;
        for (const auto& i : bitset) {
            ++num;
        }
        ASSERT_EQ(num, 0);
    }

    TEST(AIterableBitSet, IsSetIndexOutOfBoundsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        auto bitset = IBS(1);
        ASSERT_DEATH(bitset.isSet(1), "");
    }

    TEST(AIterableBitSet, RAOperatorIndexOutOfBoundsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        auto bitset = IBS(1);
        ASSERT_DEATH(bitset[1], "");
    }

    TEST(AIterableBitSet, SetTrueIndexOutOfBoundsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        auto bitset = IBS(1);
        ASSERT_DEATH(bitset.set_true(1), "");
    }

    TEST(AIterableBitSet, SetFalseIndexOutOfBoundsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        auto bitset = IBS(1);
        ASSERT_DEATH(bitset.set_false(1), "");
    }

} // namespace mt_kahypar::ds