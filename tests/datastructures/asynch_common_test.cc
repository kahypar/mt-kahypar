//
// Created by mlaupichler on 07.05.21.
//
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch/asynch_common.h"

namespace mt_kahypar::ds {

    TEST(AContractionToNodeIteratorAdapter, EmptyIterate) {
        auto vec = std::vector<Contraction>();
        auto begin = ContractionToNodeIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIteratorAdaptor(vec.end());

        ASSERT(std::distance<ContractionToNodeIteratorAdaptor>(begin, end) == 0);
        ASSERT_EQ(begin, end);
    }

    TEST(AContractionToNodeIteratorAdapter, IterateAllSameValue) {

        Contraction con = {0,1};
        auto vec = std::vector<Contraction>(5,con);
        auto begin = ContractionToNodeIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIteratorAdaptor(vec.end());

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
        auto begin = ContractionToNodeIteratorAdaptor(vec.begin());
        auto end = ContractionToNodeIteratorAdaptor(vec.end());

        HypernodeID count = 0;
        for (auto cur = begin; cur != end; ++cur) {
            ASSERT_EQ(*cur,count);
            ++count;
        }
        ASSERT_EQ(count,5);
    }

} // end namespace

