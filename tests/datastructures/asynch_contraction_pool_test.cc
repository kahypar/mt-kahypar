//
// Created by mlaupichler on 19.04.21.
//

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch/asynch_contraction_pool.h"
#include "mt-kahypar/datastructures/asynch/mock_group_hierarchy.h"

#include "tests/datastructures/hypergraph_fixtures.h"

namespace mt_kahypar::ds {

    using ::testing::Return;
    using ::testing::_;

    TEST(ASequentialContractionGroupPool, RootsAreActivatedCorrectly) {

        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0,1});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(2));

        auto mockEmpty = parallel::scalable_vector<ContractionGroupID>();
        auto mockEmptyRange = ContractionGroupIDIteratorRange(mockEmpty.begin(),mockEmpty.end());
        EXPECT_CALL(*mockGroupHierarchy,successors(_)).Times(2).WillRepeatedly(Return(mockEmptyRange));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));

        ASSERT(pool.getNumTotal() == 2);
        ASSERT(pool.getNumActive() == 2);

        bool seen0 = false;
        bool seen1 = false;

        auto picked1 = pool.pickAnyActiveID();
        ASSERT(picked1 == 0 || picked1 == 1);
        if (picked1 == 0)
            seen0 = true;
        else
            seen1 = true;

        ASSERT(pool.getNumTotal() == 2);
        ASSERT(pool.getNumActive() == 1);

        pool.activateSuccessors(picked1);
        ASSERT(pool.getNumActive() == 1);

        auto picked2 = pool.pickAnyActiveID();
        ASSERT(picked2 == 0 || picked2 == 1);
        if (picked2 == 0)
            seen0 = true;
        else
            seen1 = true;

        ASSERT(pool.getNumTotal() == 2);
        ASSERT(pool.getNumActive() == 0);

        pool.activateSuccessors(picked2);
        ASSERT(pool.getNumActive() == 0);

        ASSERT(seen0 && seen1);

    }

    TEST(ASequentialContractionGroupPool, ActivatingSuccessorsWorks) {
        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0,1});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(5));

        auto mockSuccessorsOf0 = parallel::scalable_vector<ContractionGroupID>({2,3});
        auto mockSuccessorsOf0Range = ContractionGroupIDIteratorRange(mockSuccessorsOf0.begin(),mockSuccessorsOf0.end());
        EXPECT_CALL(*mockGroupHierarchy,successors(0)).Times(1).WillOnce(Return(mockSuccessorsOf0Range));
        auto mockSuccessorsOf1 = parallel::scalable_vector<ContractionGroupID>({4});
        auto mockSuccessorsOf1Range = ContractionGroupIDIteratorRange(mockSuccessorsOf1.begin(),mockSuccessorsOf1.end());
        EXPECT_CALL(*mockGroupHierarchy,successors(1)).Times(1).WillOnce(Return(mockSuccessorsOf1Range));
        auto mockEmpty = parallel::scalable_vector<ContractionGroupID>();
        auto mockEmptyRange = ContractionGroupIDIteratorRange(mockEmpty.begin(),mockEmpty.end());
        EXPECT_CALL(*mockGroupHierarchy,successors(2)).Times(1).WillOnce(Return(mockEmptyRange));
        EXPECT_CALL(*mockGroupHierarchy,successors(3)).Times(1).WillOnce(Return(mockEmptyRange));
        EXPECT_CALL(*mockGroupHierarchy,successors(4)).Times(1).WillOnce(Return(mockEmptyRange));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));

        ASSERT_EQ(pool.getNumTotal(), 5);
        bool seen[5];
        for (bool & i : seen) {
            i = false;
        }

        auto expectedNumActive = 2;
        while (pool.hasActive()) {
            ASSERT_EQ(pool.getNumActive(), expectedNumActive);
            auto picked = pool.pickAnyActiveID();
            seen[picked] = true;
            pool.activateSuccessors(picked);
            // 0, has 2 more successors, 1 has one more successor, others don't have any
            if (picked == 0) {
                expectedNumActive += 2;
            } else if (picked == 1) {
                expectedNumActive += 1;
            }
            --expectedNumActive;
        }

        ASSERT_EQ(pool.getNumActive(),0);
        for (bool & i : seen) {
            ASSERT(i);
        }

    }

    TEST(ASequentialContractionGroupPool, ActivatingSuccessorsWhenParentStillActiveDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";

        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(1));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));
        ASSERT_EQ(pool.getNumTotal(), 1);
        ASSERT_TRUE(pool.hasActive());
        ASSERT_DEATH(pool.activateSuccessors(0), "");
    }

    TEST(ASequentialContractionGroupPool, ReactivatingTest) {
        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(1));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));
        ASSERT_EQ(pool.getNumTotal(), 1);

        ASSERT_TRUE(pool.hasActive());
        auto picked = pool.pickAnyActiveID();
        ASSERT_FALSE(pool.hasActive());
        pool.reactivate(picked);
        ASSERT_TRUE(pool.hasActive());
    }

    TEST(ASequentialContractionGroupPool, ReactivatingActiveDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";

        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(1));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));
        ASSERT_EQ(pool.getNumTotal(), 1);

        auto picked = pool.pickAnyActiveID();
        pool.reactivate(picked);
        ASSERT_DEATH(pool.reactivate(picked), "");

    }

    TEST(ASequentialContractionGroupPool, ReactivatingWhenSuccessorsAlreadyActivatedDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";

        auto mockGroupHierarchy = std::make_unique<MockGroupHierarchy>();
        auto mockRoots = parallel::scalable_vector<ContractionGroupID>({0});
        auto mockRootRange = ContractionGroupIDIteratorRange(mockRoots.begin(),mockRoots.end());
        EXPECT_CALL(*mockGroupHierarchy,roots()).Times(1).WillOnce(Return(mockRootRange));
        EXPECT_CALL(*mockGroupHierarchy,getNumGroups()).WillRepeatedly(Return(1));
        auto mockEmpty = parallel::scalable_vector<ContractionGroupID>();
        auto mockEmptyRange = ContractionGroupIDIteratorRange(mockEmpty.begin(),mockEmpty.end());
        EXPECT_CALL(*mockGroupHierarchy,successors(0)).Times(1).WillOnce(Return(mockEmptyRange));

        auto pool = SequentialContractionGroupPool(std::move(mockGroupHierarchy));
        ASSERT_EQ(pool.getNumTotal(), 1);

        auto picked = pool.pickAnyActiveID();
        ASSERT_EQ(picked,0);
        pool.activateSuccessors(picked);
        ASSERT_DEATH(pool.reactivate(picked), "");
    }

} // namespace mt_kahypar
