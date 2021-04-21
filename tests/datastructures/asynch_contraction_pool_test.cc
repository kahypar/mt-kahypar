//
// Created by mlaupichler on 19.04.21.
//

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch_contraction_pool.h"

#include "tests/datastructures/hypergraph_fixtures.h"

namespace mt_kahypar::ds {

TEST(AAsynchContractionPool,InsertGroup) {

    auto pool = AsynchContractionPool();
    Contraction c1 = {1,2};
    Contraction c2 = {1,4};
    Contraction c3 = {1,12567};
    ContractionGroup group = {c1,c2,c3};

    pool.insertContractionGroup(group);

    ASSERT_FALSE(pool.empty());
    ASSERT_EQ(pool.unsafe_size(),1);

}

TEST(AAsynchContractionPool,InsertSingleContraction) {

    auto pool = AsynchContractionPool();
    Contraction c1 = {1,2};

    pool.insertContraction(c1);

    ASSERT_FALSE(pool.empty());
    ASSERT_EQ(pool.unsafe_size(),1);

}

TEST(AAsynchContractionPool, ContainsTest){

        auto pool = AsynchContractionPool();
        Contraction c1 = {1,2};
        Contraction c2 = {1,4};
        Contraction c3 = {1,3};
        ContractionGroup group1 = {c1,c2};
        ContractionGroup group2 = {c3};

        pool.insertContractionGroup(group1);
        pool.insertContractionGroup(group2);

        ASSERT_TRUE(pool.contains(group1));
        ASSERT_TRUE(pool.contains(group2));

        ASSERT_TRUE(pool.contains(c1));
        ASSERT_TRUE(pool.contains(c2));
        ASSERT_TRUE(pool.contains(c3));

        //Expect containing reordered groups
        ContractionGroup reorderOfGroup1 = {c2,c1};
        ASSERT_TRUE(pool.contains(reorderOfGroup1));

        // Expect not containing completely different contraction
        Contraction c4 = {3,5};
        ASSERT_FALSE(pool.contains(c4));
        ContractionGroup group3 = {c4};
        ASSERT_FALSE(pool.contains(group3));

        // Expect not containing group that is different combination of known contractions
        ContractionGroup unknownGroup = {c1,c3};
        pool.debugPrint();
        std::cout << "Unknown group is: \n";
        unknownGroup.debugPrint();
        ASSERT_FALSE(pool.contains(unknownGroup));



}

TEST(AAsynchContractionPool, PickAnyTest) {

    auto pool = AsynchContractionPool();
    Contraction c1 = {1,2};
    Contraction c2 = {1,4};
    Contraction c3 = {2,3};
    Contraction c4 = {3,5};
    ContractionGroup group1 = {c1,c2};
    ContractionGroup group2 = {c3};

    pool.insertContractionGroup(group1);
    pool.insertContractionGroup(group2);
    pool.insertContraction(c4);

    auto poolCopy = pool;

    auto picked = pool.pickAnyGroup();

    ASSERT_FALSE(pool.contains(picked));
    ASSERT_TRUE(poolCopy.contains(picked));
}

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

TEST(AAsynchContractionPoolDeathTest, InsertFailsIfEmptyDeathTest) {
    testing::FLAGS_gtest_death_test_style="threadsafe";
    auto pool = AsynchContractionPool();
    ContractionGroup group;

    EXPECT_DEATH(pool.insertContractionGroup(group),"");
}

TEST(AContractionGroupDeathTest, ConstructionFailsIfGroupHasDifferentRepsDeathTest) {
    testing::FLAGS_gtest_death_test_style="threadsafe";
    auto pool = AsynchContractionPool();
    Contraction c1 = {1,2};
    Contraction c2 = {2,3};
    ASSERT_DEATH(ContractionGroup({c1,c2}),"");
}

} // namespace mt_kahypar
