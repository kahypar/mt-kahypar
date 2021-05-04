//
// Created by mlaupichler on 04.05.21.
//


#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch/array_lock_manager.h"

namespace mt_kahypar::ds {

    TEST(AArrayLockManager,InitializesToRightSizeAndAllInvalidWithUInt32_T) {

        using LockedID = uint32_t;
        using OwnerID = uint32_t;

        LockedID size = 10;
        OwnerID invalid = std::numeric_limits<OwnerID>::max();

        auto lockManager = ArrayLockManager<LockedID, OwnerID>(size, invalid);

    }

} // end namespace mt_kahypar::ds