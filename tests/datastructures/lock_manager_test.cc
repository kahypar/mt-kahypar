//
// Created by mlaupichler on 04.05.21.
//


#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch/array_lock_manager.h"

namespace mt_kahypar::ds {

    using LockedID = uint32_t;
    using OwnerID = uint32_t;
    OwnerID defaultInvalid = std::numeric_limits<OwnerID>::max();

    TEST(AArrayLockManager,InitializeAndEmptyAccessWithSameType) {

        LockedID size = 10;

        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);

        OwnerID testOwner = 0;
        for (LockedID i = 0; i < size; ++i) {
            ASSERT_EQ(lockManager.owner(i),defaultInvalid);
            ASSERT_FALSE(lockManager.isLocked(i));
            ASSERT_EQ(lockManager[i],defaultInvalid);
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
            ASSERT_EQ(lockManager.owner(i),invalid);
            ASSERT_FALSE(lockManager.isLocked(i));
            ASSERT_EQ(lockManager[i],invalid);
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
        ASSERT_EQ(lockManager.owner(0), 0);
        ASSERT_EQ(lockManager.owner(5), 1);
        ASSERT_EQ(lockManager.owner(9), 10000);

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

    TEST(AArrayLockManager, UnqualifiedIndexOwnerDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager.owner(1), "");
    }

    TEST(AArrayLockManager, UnqualifiedIndexRAOperatorDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        LockedID size = 1;
        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, defaultInvalid);
        ASSERT_DEATH(lockManager[2], "");
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


} // end namespace mt_kahypar::ds