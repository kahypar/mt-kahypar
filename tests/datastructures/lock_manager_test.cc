//
// Created by mlaupichler on 04.05.21.
//


#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/async/array_lock_manager.h"

namespace mt_kahypar::ds {

    using LockedID = uint32_t;
    using OwnerID = uint32_t;
    OwnerID defaultInvalid = std::numeric_limits<OwnerID>::max();

    using LockedIDIterator = std::vector<LockedID>::const_iterator;

    IteratorRange<LockedIDIterator> getRange(LockedID first, LockedID firstInvalid) {

        ASSERT(firstInvalid >= first);
        auto diff = firstInvalid - first;
        auto vec = new std::vector<LockedID>(diff);
        for (LockedID i = 0; i < diff; ++i) {
            (*vec)[i] = i + first;
        }

        return IteratorRange<LockedIDIterator>(vec->begin(),vec->end());
    }

    TEST(AArrayLockManager,InitializeAndEmptyAccessWithSameType) {

        LockedID size = 10;

        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);

        OwnerID testOwner = 0;
        for (LockedID i = 0; i < size; ++i) {
            ASSERT_FALSE(lockManager.isLocked(i));
            ASSERT_TRUE(lockManager.isHeldBy(i,defaultInvalid));
            ASSERT_FALSE(lockManager.tryToReleaseLock(i,testOwner));
        }
    }

    TEST(AArrayLockManager,InitializeAndEmptyAccessWithDifferingTypes) {

        using DiffLockedID = uint16_t;
        using DiffOwnerID = uint8_t;

        DiffLockedID size = std::numeric_limits<DiffLockedID>::max();
        DiffOwnerID invalid = std::numeric_limits<DiffOwnerID>::max();

        auto lockManager = ArrayLockManager<DiffLockedID, DiffOwnerID>(size, invalid);

        DiffOwnerID testOwner = 0;
        ASSERT_NE(testOwner,invalid);

        for (DiffLockedID i = 0; i < size; ++i) {
            ASSERT_FALSE(lockManager.isLocked(i));
            ASSERT_TRUE(lockManager.isHeldBy(i,invalid));
            ASSERT_FALSE(lockManager.tryToReleaseLock(i,testOwner));
        }
    }

    TEST(AArrayLockManager, SequentialAcquireReleaseTest) {
        LockedID size = 10;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size,defaultInvalid);

        // Owner 0 acquires lock 0, owner 1 acquires lock 5, owner 10000 acquires lock 9
        auto acquire1 = lockManager.tryToAcquireLock(0,0);
        auto acquire2 = lockManager.tryToAcquireLock(5,1);
        auto acquire3 = lockManager.tryToAcquireLock(9,10000);

        ASSERT_TRUE(acquire1 && acquire2 && acquire3);
        ASSERT_TRUE(lockManager.isLocked(0) && lockManager.isLocked(5) && lockManager.isLocked(9));
        ASSERT_TRUE(lockManager.isHeldBy(0,0));
        ASSERT_TRUE(lockManager.isHeldBy(5,1));
        ASSERT_TRUE(lockManager.isHeldBy(9,10000));

        // Owner 2 holds no locks so it cannot release any of the three locks and also not a lock that is held by no one (lock 1)
        ASSERT_FALSE(lockManager.tryToReleaseLock(0,2));
        ASSERT_FALSE(lockManager.tryToReleaseLock(5,2));
        ASSERT_FALSE(lockManager.tryToReleaseLock(9,2));
        ASSERT_FALSE(lockManager.tryToReleaseLock(1,2));

        // Lock 5 can be released by owner 1 and then acquired by owner 2
        auto release1 = lockManager.tryToReleaseLock(5,1);
        ASSERT_TRUE(release1);
        ASSERT_FALSE(lockManager.isLocked(5));
        auto acquire4 = lockManager.tryToAcquireLock(5,2);
        ASSERT_TRUE(acquire4);
        ASSERT_TRUE(lockManager.isLocked(5));

        // Release all locks (in no specific order)
        auto release2 = lockManager.tryToReleaseLock(0,0);
        auto release3 = lockManager.tryToReleaseLock(5,2);
        auto release4 = lockManager.tryToReleaseLock(9,10000);
        ASSERT_TRUE(release2 && release3 && release4);
        ASSERT_FALSE(lockManager.isLocked(0) || lockManager.isLocked(5) || lockManager.isLocked(9));

    }

    TEST(AArrayLockManager, AcquireWhenAlreadyHoldingLockTest) {
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size,defaultInvalid);

        // Acquire lock that is not being held
        ASSERT_TRUE(lockManager.tryToAcquireLock(0,0));
        // Acquire lock again with same owner while already holding it, expect the return value to be false
        ASSERT_FALSE(lockManager.tryToAcquireLock(0,0));
    }

    TEST(AArrayLockManager, ReleaseWhenNotHoldingLockTest) {
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size,defaultInvalid);

        ASSERT_FALSE(lockManager.isLocked(0));
        // Release lock while not holding it, expect the return value to be false
        ASSERT_FALSE(lockManager.tryToReleaseLock(0,0));
    }

    TEST(AArrayLockManager, AcquireMultipleTest) {
        LockedID size = 10;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size,defaultInvalid);

        // Owner 0 tries to acquire 0 to 4 at once, expect it to work and that owner 0 owns all of them after
        auto range1 = getRange(0,5);
        bool acquire1 = lockManager.tryToAcquireMultipleLocks(range1, 0);
        ASSERT_TRUE(acquire1);
        for (LockedID i : range1) {
            ASSERT_TRUE(lockManager.isLocked(i));
            ASSERT_TRUE(lockManager.isHeldBy(i,0));
        }

        // Owner 1 tries to acquire 8 to 9 at once, expect it work and that owner 1 owns 8 and 9 after
        auto range4 = getRange(8,10);
        bool acquire4 = lockManager.tryToAcquireMultipleLocks(range4, 1);
        ASSERT_TRUE(acquire4);
        for (auto i : range4) {
            ASSERT_TRUE(lockManager.isLocked(i));
            ASSERT_TRUE(lockManager.isHeldBy(i,1));
        }

        // Owner 0 tries to acquire 4 to 6 at once, expect it to fail and that owner 0 owns 4 still but that 5 and 6 are still free after
        auto range2 = getRange(4,7);
        bool acquire2 = lockManager.tryToAcquireMultipleLocks(range2, 0);
        ASSERT_FALSE(acquire2);
        ASSERT_TRUE(lockManager.isHeldBy(4,0));
        ASSERT_FALSE(lockManager.isLocked(5) || lockManager.isLocked(6));

        // Owner 1 tries to acquire 4 to 6 at once, expect it to fail and that owner 0 owns 4 still but that 5 and 6 are still free after
        auto range3 = getRange(4,7);
        bool acquire3 = lockManager.tryToAcquireMultipleLocks(range3, 1);
        ASSERT_FALSE(acquire3);
        ASSERT_TRUE(lockManager.isHeldBy(4,0));
        ASSERT_FALSE(lockManager.isLocked(5) || lockManager.isLocked(6));

        // Owner 2 tries to acquire 5 to 9 at once, expect it to fail and that owner 1 owns 8 to 9 still, and 5 to 7 are still free
        auto range5 = getRange(5,10);
        bool acquire5 = lockManager.tryToAcquireMultipleLocks(range5, 2);
        ASSERT_FALSE(acquire5);
        auto freeRange = getRange(5,8);
        for (auto i : freeRange) {
            ASSERT_FALSE(lockManager.isLocked(i));
        }
        ASSERT_TRUE(lockManager.isHeldBy(8,1));
        ASSERT_TRUE(lockManager.isHeldBy(9,1));

    }

    TEST(AArrayLockManager, ReleaseMultipleTest) {

        LockedID size = 10;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size,defaultInvalid);

        // Try to release multiple fails on empty
        auto range1 = getRange(0,5);
        auto range2 = getRange(7,8);
        ASSERT_FALSE(lockManager.tryToReleaseMultipleLocks(range1, 0));
        ASSERT_FALSE(lockManager.tryToReleaseMultipleLocks(range2, 4));

        // Individual acquire of 0 to 4 by owner 0
        auto range3 = getRange(0,5);
        for (auto i : range3) {
            bool locked = lockManager.tryToAcquireLock(i,0);
            ASSERT_TRUE(locked);
            ASSERT_TRUE(lockManager.isHeldBy(i,0));
        }

        // Releasing multiple with different owner fails
        bool release1 = lockManager.tryToReleaseMultipleLocks(range3,1);
        ASSERT_FALSE(release1);
        for (auto i : range3) {
            ASSERT_TRUE(lockManager.isHeldBy(i,0));
        }

        // Releasing when not all are owned fails and no locks are changed
        auto wrongRange = getRange(3,8);
        bool release2 = lockManager.tryToReleaseMultipleLocks(wrongRange,0);
        ASSERT_FALSE(release2);
        for (auto i : range3) {
            ASSERT_TRUE(lockManager.isHeldBy(i,0));
        }

        // Releasing all previously locked works and afterwards they are all free
        bool release3 = lockManager.tryToReleaseMultipleLocks(range3,0);
        ASSERT_TRUE(release3);
        for (auto i : range3) {
            ASSERT_FALSE(lockManager.isLocked(i));
        }

        // Individual acquire of 7 to 9 by owner 3
        auto range4 = getRange(7,10);
        for (auto i : range4) {
            bool locked = lockManager.tryToAcquireLock(i,3);
            ASSERT_TRUE(locked);
            ASSERT_TRUE(lockManager.isHeldBy(i,3));
        }

        // Releasing part of previously acquired works and keeps rest untouched
        auto partRange = getRange(7,9);
        bool release4 = lockManager.tryToReleaseMultipleLocks(partRange,3);
        ASSERT_TRUE(release4);
        ASSERT_TRUE(lockManager.isHeldBy(9,3));
        ASSERT_FALSE(lockManager.isLocked(7) || lockManager.isLocked(8));

        // Releasing single lock works
        auto partRange2 = getRange(9,10);
        bool release5 = lockManager.tryToReleaseMultipleLocks(partRange2,3);
        ASSERT_TRUE(release5);
        ASSERT_FALSE(lockManager.isLocked(9));

    }

    TEST(AArrayLockManager, UnqualifiedIndexIsHeldByDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.isHeldBy(1,0), "");
    }

    TEST(AArrayLockManager, UnqualifiedIndexIsLockedDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.isLocked(2), "");
    }

    TEST(AArrayLockManager, UnqualifiedIndexTryToAcquireDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        OwnerID testOwner = 0;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.tryToAcquireLock(100,testOwner), "");
    }

    TEST(AArrayLockManager, UnqualifiedIndexTryToReleaseDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        OwnerID testOwner = 0;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.tryToReleaseLock(100,testOwner), "");
    }

    TEST(AArrayLockManager, EmptyLockRangeTryToAcquireMultipleDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        OwnerID testOwner = 0;
        auto emptyRange = getRange(0,0);
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.tryToAcquireMultipleLocks(emptyRange,testOwner), "");
    }

    TEST(AArrayLockManager, EmptyLockRangeTryToReleaseMultipleDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        OwnerID testOwner = 0;
        auto emptyRange = getRange(0,0);
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.tryToReleaseMultipleLocks(emptyRange,testOwner), "");
    }

//    TEST(ILockManager, ReleaseMultipleWithReleaseFailAfterOwnerCheckDeathTest) {
//
//        using ::testing::_;
//        using ::testing::Return;
//        testing::FLAGS_gtest_death_test_style="threadsafe";
//
//        OwnerID ownerID = 0;
//
//        // Build mock manager who will return true to any isHeldBy query with owner 0 but will return false on any attempt to release
//        // a lock held by 0. Simulates situation where a owner holds a lock but releasing it in tryToReleaseMultiple
//        // still fails (the program is supposed to fail then).
//        auto mockLockManager = MockLockManager<LockedID, OwnerID>();
//        EXPECT_CALL(mockLockManager,isHeldBy(_,ownerID)).WillRepeatedly(Return(true));
//        EXPECT_CALL(mockLockManager,tryToReleaseLock(_,ownerID)).WillRepeatedly(Return(false));
//
//        auto range = getRange(0,5);
//        ASSERT_DEATH(mockLockManager.tryToReleaseMultipleLocks(range,ownerID),"");
//
//    }


} // end namespace mt_kahypar::ds