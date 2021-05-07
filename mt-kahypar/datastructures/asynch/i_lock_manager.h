//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_I_LOCK_MANAGER_H
#define KAHYPAR_I_LOCK_MANAGER_H

#include <kahypar/meta/mandatory.h>
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/utils/range.h"

#include "gmock/gmock.h"

namespace mt_kahypar::ds {

    /// Interface for managing locks. Locked objects are marked by LockedIDs and mapped to OwnerIDs describing
    /// who holds the lock. Requires LockedID and OwnerID to be unsigned integral number types.
    template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
    class ILockManager {

    public:

        /// Attempt to acquire a lock for the id lockedID with the owner ownerID. Will return true if the lock owner has
        /// changed and false if not.
        virtual bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) = 0;

        /// Release the lock for lockedID held by ownerID. Will return true if the lock has been successfully released or
        /// false if the given ownerID does not hold the lock (the owner does not change in the latter case).
        virtual bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) = 0;

        /// Returns the OwnerID of the current owner of lockedID or the invalid OwnerID if it is not locked.
        virtual OwnerID owner(LockedID lockedID) const = 0;

        /// Returns true if any owner holds the lock for lockedID.
        virtual bool isLocked(LockedID lockedID) const = 0;

        OwnerID operator[] (const LockedID lockedID) const {
            return owner(lockedID);
        };

        /// Attempt to acquire multiple locks for the ids in the range lockedIDs with the owner ownerID.
        /// Will return true if all requested locks were changed and acquired by the owner with ownerID. If any
        /// lock cannot be acquired, all acquired locks will be released again and false will be returned.
        template<typename LockedIDIteratorT>
        bool tryToAcquireMultipleLocks(IteratorRange<LockedIDIteratorT> lockedIDs, OwnerID ownerID) {
            ASSERT(std::distance(lockedIDs.begin(),lockedIDs.end()) > 0 && "Range of ids to lock cannot be empty!");

            auto cur = lockedIDs.begin();
            auto success = true;
            while (cur != lockedIDs.end() && success) {
                success &= tryToAcquireLock(*cur, ownerID);
                ++cur;
            }

            if (success) {
                return true;
            } else {
                // Release locks acquired so far
                auto rev = lockedIDs.begin();
                auto endOfAcquired = --cur;
                while (rev != endOfAcquired) {
                    bool revSuccess = tryToReleaseLock(*rev, ownerID);
                    ASSERT(revSuccess);
                    ++rev;
                }
                return false;
            }
        };

        /// Attempt to release multiple locks for the ids in the range lockedIDs with the owner ownerID.
        /// Will return true if all requested locks were changed and released. If the owner does not hold all of the given locks,
        /// no locks will be released and false will be returned. If any lock cannot be released for a different reason
        /// false will be returned and no guarantees can be given on which locks are still held by the owner and which are not.
        template<typename LockedIDIteratorT>
        bool tryToReleaseMultipleLocks(IteratorRange<LockedIDIteratorT> lockedIDs, OwnerID ownerID) {
            ASSERT(std::distance(lockedIDs.begin(),lockedIDs.end()) > 0 && "Range of locked ids cannot be empty!");

            // Check if any given lock is not held by the given owner. If so, safe return false.
            for (auto i : lockedIDs) {
                if (owner(i) != ownerID) return false;
            }

            auto cur = lockedIDs.begin();
            while (cur != lockedIDs.end()) {
                bool released = tryToReleaseLock(*cur, ownerID);
                ASSERT(released);
                // If release fails at this point it has a different reason than ownerID not being the right owner.
                // If asserts are disabled, in that case one can only return false as already released
                // locks may not be able to be re-acquired (unsafe return false).
                if (!released) return false;
                ++cur;
            }
            return true;
        };

    protected:

        ILockManager(OwnerID invalidOwnerID) :
                _invalid_owner_id(invalidOwnerID) {};

        virtual ~ILockManager() = default;

        OwnerID _invalid_owner_id;

    };

    template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
    class MockLockManager : public ILockManager<LockedID,OwnerID> {
    public:
        MockLockManager() : ILockManager<LockedID, OwnerID>(0) {}

        MOCK_METHOD(bool, tryToAcquireLock,(LockedID lockedID, OwnerID ownerID),(override));
        MOCK_METHOD(bool, tryToReleaseLock,(LockedID lockedID, OwnerID ownerID),(override));
        MOCK_METHOD(OwnerID, owner,(LockedID lockedID), (const, override));
        MOCK_METHOD(bool, isLocked,(LockedID lockedID), (const, override));
    };


} // namespace mt_kahypar::ds




#endif //KAHYPAR_I_LOCK_MANAGER_H
