//
// Created by mlaupichler on 07.05.21.
//
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/async/async_common.h"

namespace mt_kahypar::ds {

    TEST(AContractionGroup, EqualityTest) {
        Contraction c1 = {1,2};
        Contraction c2 = {1,4};
        Contraction c3 = {1,3};
        ContractionGroup group1 = {c1,c2};
        ContractionGroup group2 = {c1,c3};
        ContractionGroup group3 = {c2,c1};

        ASSERT_FALSE(group1 == group2);
        ASSERT_TRUE(group1 == group3);

        ContractionGroup group4 = {c1};
        ASSERT_FALSE(group1 == group4);
        ASSERT_FALSE(group4 == group1);
    }

    TEST(AContractionGroupDeathTest, ConstructFailsIfEmptyDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";

        EXPECT_DEATH(ContractionGroup group = {},"");
        std::vector<Contraction> emptyVec;
        EXPECT_DEATH(ContractionGroup group(std::move(emptyVec)),"");
    }

    TEST(AContractionGroupDeathTest, ConstructionFailsIfGroupHasDifferentRepsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        Contraction c1 = {1,2};
        Contraction c2 = {2,3};
        ASSERT_DEATH(ContractionGroup({c1,c2}),"");
    }

    TEST(AContractionToNodeIteratorAdapter, EmptyIterate) {
        auto vec = std::vector<Contraction>();
        auto begin = ContractionToNodeIDIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIDIteratorAdaptor(vec.end());

        ASSERT(std::distance<ContractionToNodeIDIteratorAdaptor>(begin, end) == 0);
        ASSERT_EQ(begin, end);
    }

    TEST(AContractionToNodeIteratorAdapter, IterateAllSameValue) {

        Contraction con = {0,1};
        auto vec = std::vector<Contraction>(5,con);
        auto begin = ContractionToNodeIDIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIDIteratorAdaptor(vec.end());

        int count = 0;
        for (auto cur = begin; cur != end; ++cur) {
            ++count;
            ASSERT_EQ(*cur,1);
        }
        ASSERT_EQ(count,5);
    }

    TEST(AContractionToNodeIteratorAdapter, IterateDifferentValues) {

        auto vec = std::vector<Contraction>(5);
        for (HypernodeID i = 0; i < vec.size(); ++i) {
            vec[i] = {0, i};
        }
        auto begin = ContractionToNodeIDIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIDIteratorAdaptor(vec.end());

        HypernodeID count = 0;
        for (auto cur = begin; cur != end; ++cur) {
            ASSERT_EQ(*cur,count);
            ++count;
        }
        ASSERT_EQ(count,5);
    }

    TEST(AContractionToNodeIteratorAdapter, IterateOnGroup) {
        ContractionGroup group = {{0,1},{0,2},{0,3},{0,4}};
        auto begin = ContractionToNodeIDIteratorAdaptor(group.begin());
        auto end = ContractionToNodeIDIteratorAdaptor(group.end());

        bool seen[4] = {false,false,false,false};
        for (auto cur = begin; cur != end; ++cur) {
            HypernodeID val = *cur;
            ASSERT(!seen[val-1]);
            seen[val-1] = true;
        }

        for (size_t i = 0; i < 4; ++i) {
            ASSERT_TRUE(seen[i]);
        }
    }

    TEST(AGroupNodeIDIterator, EqualityAndInequality) {

        ContractionGroup group1 = {{0,1},{0,2},{0,3},{0,4}};
        auto begin1 = GroupNodeIDIterator::getAtBegin(group1);
        auto end1 = GroupNodeIDIterator::getAtEnd(group1);
        ASSERT_TRUE(begin1 != end1);
        auto it1 = begin1;
        ASSERT_TRUE(it1 == begin1);
        std::advance(it1,5);
        ASSERT_TRUE(it1 == end1);

        // Test that for identical other group, iterators are not seen as equal
        ContractionGroup group2 = {{0,1},{0,2},{0,3},{0,4}};
        auto begin2 = GroupNodeIDIterator::getAtBegin(group2);
        auto end2 = GroupNodeIDIterator::getAtEnd(group2);
        ASSERT_TRUE(begin2 != begin1);
        ASSERT_TRUE(end2 != end1);

    }

    TEST(AGroupNodeIDIterator, IterateSingleContractionInGroup) {

        ContractionGroup group = {{0,1}};
        // Iteration should hit 0 and 1
        bool seen0, seen1 = false;

        auto begin = GroupNodeIDIterator::getAtBegin(group);
        auto end = GroupNodeIDIterator::getAtEnd(group);

        for (auto it = begin; it != end; ++it) {
            ASSERT_TRUE((*it == 0 && !seen0) || (*it == 1 && !seen1));
            if (*it == 0) {
                seen0 = true;
            } else if (*it == 1) {
                seen1 = true;
            }
        }
        ASSERT_TRUE(seen0 && seen1);
    }

    TEST(AGroupNodeIDIterator, IterateMultipleContractionsInGroup) {

        ContractionGroup group = {{0,1},{0,2},{0,3},{0,4}};
        auto begin = GroupNodeIDIterator::getAtBegin(group);
        auto end = GroupNodeIDIterator::getAtEnd(group);

        bool seen[5] = {false,false,false,false,false};
        for (auto cur = begin; cur != end; ++cur) {
            HypernodeID val = *cur;
            ASSERT(!seen[val]);
            seen[val] = true;
        }

        for (size_t i = 0; i < 5; ++i) {
            ASSERT_TRUE(seen[i]);
        }
    }

    TEST(APinSnapshotIterator, StitchingWorksWithNonEmptySequences) {

      Array<HypernodeID> array_stable(5, 0);
      for (HypernodeID i = 0; i < array_stable.size(); ++i) {
        array_stable[i] = i;
      }

      HypernodeID array_volatile[5] = {0, 1, 2, 3, 4};

      auto stable_range = IteratorRange<Array<HypernodeID>::const_iterator>(array_stable.begin(), array_stable.begin() + 3);
      auto volatile_range = IteratorRange<HypernodeID*>(array_volatile + 3, array_volatile + 5);
      IteratorRange<PinSnapshotIterator> stitched = PinSnapshotIterator::stitchPinIterators(stable_range, volatile_range);

      HypernodeID j = 0;
      for (auto it = stitched.begin(); it != stitched.end(); ++it) {
        ASSERT_EQ(*it, j);
        ++j;
      }

      auto it = stitched.begin();
      it++; it++;
      ASSERT_EQ(*it, 2);
    }

    TEST(APinSnapshotIterator, StitchingWorksWithEmptyStable) {

      Array<HypernodeID> array_stable(0, 0);

      HypernodeID array_volatile[5] = {0, 1, 2, 3, 4};

      auto stable_range = IteratorRange<Array<HypernodeID>::const_iterator>(array_stable.begin(), array_stable.end());
      auto volatile_range = IteratorRange<HypernodeID*>(array_volatile, array_volatile + 5);
      IteratorRange<PinSnapshotIterator> stitched = PinSnapshotIterator::stitchPinIterators(stable_range, volatile_range);

      HypernodeID j = 0;
      for (auto it = stitched.begin(); it != stitched.end(); ++it) {
        ASSERT_EQ(*it, j);
        ++j;
      }

      auto it = stitched.begin();
      it++; it++;
      ASSERT_EQ(*it, 2);
    }

    TEST(APinSnapshotIterator, StitchingWorksWithEmptyVolatile) {

      Array<HypernodeID> array_stable(5, 0);
      for (HypernodeID i = 0; i < array_stable.size(); ++i) {
        array_stable[i] = i;
      }

      HypernodeID array_volatile[1] = {0};

      auto stable_range = IteratorRange<Array<HypernodeID>::const_iterator>(array_stable.begin(), array_stable.end());
      auto volatile_range = IteratorRange<HypernodeID*>(array_volatile, array_volatile);
      IteratorRange<PinSnapshotIterator> stitched = PinSnapshotIterator::stitchPinIterators(stable_range, volatile_range);

      HypernodeID j = 0;
      for (auto it = stitched.begin(); it != stitched.end(); ++it) {
        ASSERT_EQ(*it, j);
        ++j;
      }

      auto it = stitched.begin();
      it++; it++;
      ASSERT_EQ(*it, 2);
    }

    TEST(APinSnapshotIterator, StitchingWorksWithBothEmpty) {

      Array<HypernodeID> array_stable(0, 0);

      HypernodeID array_volatile[1] = {0};

      auto stable_range = IteratorRange<Array<HypernodeID>::const_iterator>(array_stable.begin(), array_stable.end());
      auto volatile_range = IteratorRange<HypernodeID*>(array_volatile, array_volatile);
      IteratorRange<PinSnapshotIterator> stitched = PinSnapshotIterator::stitchPinIterators(stable_range, volatile_range);

      ASSERT_EQ(stitched.begin(), stitched.end());
    }

} // end namespace

