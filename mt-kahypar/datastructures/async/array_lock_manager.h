//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_ARRAY_LOCK_MANAGER_H
#define KAHYPAR_ARRAY_LOCK_MANAGER_H


#include <mt-kahypar/datastructures/hypergraph_common.h>
#include <mt-kahypar/utils/memory_tree.h>
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"
#include "kahypar/meta/mandatory.h"
#include "async_common.h"

namespace mt_kahypar {

    //Forward declaration
    class AsyncNodeTracker;

namespace ds {

    /// Class template for a LockManager based on an array of atomic wrappers for the given OwnerID type using the
    /// given LockedID type as indices.
    template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
    class ArrayLockManager {

        using UnderlyingType = CAtomic<OwnerID>;

    public:

        /// Constructs a new ArrayLockManager of given size with a given distinct OwnerID to be considered invalid.
        /// \param size the size of the array of locks. All LockedIDs used with an ArrayLockManager instantiated with
        /// this constructor are expected to be in [0,size).
        /// \param invalidOwnerID a distinct OwnerID that an ArrayLockManager instantiated with this constructor
        /// will consider an unassigned or invalid owner.
        ArrayLockManager(LockedID size, OwnerID invalidOwnerID)
            : _size(size), _v(std::make_unique<UnderlyingType[]>(size)), _invalid_owner_id(invalidOwnerID) {
            init();
        }

        ArrayLockManager(const ArrayLockManager&) = delete;
        ArrayLockManager& operator= (const ArrayLockManager&) = delete;

        ArrayLockManager(ArrayLockManager&&)  noexcept = default;
        ArrayLockManager& operator= (ArrayLockManager&&)  noexcept = default;

        ~ArrayLockManager() = default;

        /// Attempt to acquire a lock for the id lockedID with the owner ownerID. Will return true if the lock owner has
        /// changed and false if not.
        bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) {

            if (lockedID >= _size) {
                ERROR("Cannot acquire lock for ID larger than size of lock array.");
            }

            OwnerID expected = owner(lockedID);
            OwnerID desired = ownerID;
            if (expected != _invalid_owner_id) {
                // If the lock is already held by someone, return false
                return false;
            }
            // If the lock is currently not held by anyone, attempt to acquire it; returns true if successful or false if not
            bool acquired = _v[lockedID].compare_exchange_strong(expected, desired);

            return acquired;
        };

        /// Try to release the lock for lockedID held by ownerID. Will return true if the lock has been successfully released or
        /// false if the given ownerID does not hold the lock (the owner does not change in the latter case).
        bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) {
            ASSERT(lockedID < _size);
            if (lockedID >= _size) {
                ERROR("Cannot release lock for ID larger than size of lock array.");
            }

            OwnerID expected = ownerID;
            OwnerID desired = _invalid_owner_id;
            bool released = _v[lockedID].compare_exchange_strong(expected, desired);
            return released;
        };

        /// Returns true if any owner holds the lock for lockedID and false if not.
        bool isLocked(LockedID lockedID) const {
            return owner(lockedID) != _invalid_owner_id;
        };

        /// Checks whether a given lock is currently held by a given owner.
        bool isHeldBy(LockedID lockedID, OwnerID ownerID) const {
            return owner(lockedID) == ownerID;
        }

        /// Attempt to acquire multiple locks for the ids in the range lockedIDs with the owner ownerID.
        /// Will return true if all requested locks were changed and acquired by the owner with ownerID. If any
        /// lock cannot be acquired, all acquired locks will be released again and false will be returned.
        template<typename LockedIDIteratorT>
        bool tryToAcquireMultipleLocks(IteratorRange<LockedIDIteratorT> lockedIDs, OwnerID ownerID) {
            ASSERT(std::distance(lockedIDs.begin(),lockedIDs.end()) > 0 && "Range of ids to lock cannot be empty!");

            LockedIDIteratorT cur = lockedIDs.begin();
            LockedIDIteratorT lastAcquired = lockedIDs.begin();
            auto success = true;
            while (cur != lockedIDs.end() && success) {
                lastAcquired = cur;
                success &= tryToAcquireLock(*cur, ownerID);
                ++cur;
            }

            if (success) {
                return true;
            } else {
                // Release locks acquired so far
                auto rev = lockedIDs.begin();
                while (rev != lastAcquired) {
                    bool revSuccess = tryToReleaseLock(*rev, ownerID);
                    unused(revSuccess);
                    ASSERT(revSuccess);
                    ++rev;
                }
                return false;
            }
        }

        /// Attempt to release multiple locks for the ids in the range lockedIDs with the owner ownerID.
        /// Will return true if all requested locks were changed and released. If the owner does not hold all of the given locks,
        /// no locks will be released and false will be returned. If any lock cannot be released for a different reason
        /// false will be returned and no guarantees can be given on which locks are still held by the owner and which are not.
        template<typename LockedIDIteratorT>
        bool tryToReleaseMultipleLocks(IteratorRange<LockedIDIteratorT> lockedIDs, OwnerID ownerID) {
            ASSERT(std::distance(lockedIDs.begin(),lockedIDs.end()) > 0 && "Range of locked ids cannot be empty!");

            // Check if any given lock is not held by the given owner. If so, safe return false.
            for (auto i : lockedIDs) {
                if (!isHeldBy(i,ownerID)) return false;
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
        }

        /// Acquires lock with given lockedID and owner ownerID with the strong expectation that the operation will work.
        /// This asserts that the lock is actually acquired, i.e. it halts the program if the acquire fails.
        void strongAcquireLock(LockedID lockedID, OwnerID ownerID) {
            bool acquired = tryToAcquireLock(lockedID, ownerID);
            ASSERT(acquired && "Strong acquire of lock failed.");
            if (!acquired) {
              ERROR("strongAcquireLock() failed! Aborting.");
            }
        }

        /// Releases lock with given lockedID and owner ownerID with the strong expectation that the operation will work.
        /// This asserts that the lock is actually released, i.e. it halts the program if the release fails.
        void strongReleaseLock(LockedID lockedID, OwnerID ownerID) {
            bool released = tryToReleaseLock(lockedID, ownerID);
            ASSERT(released && "Strong release of lock failed.");
            if (!released) {
              ERROR("strongReleaseLock() failed! Aborting.");
            }
        }

        /// Releases locks with given lockedIDs and owner ownerID with the strong expectation that the operation will work.
        /// This asserts that the locks are actually all released, i.e. it halts the program if the release fails.
        template<typename LockedIDIteratorT>
        void strongReleaseMultipleLocks(IteratorRange<LockedIDIteratorT> lockedIDs, OwnerID ownerID) {
            bool released = tryToReleaseMultipleLocks(lockedIDs, ownerID);
            ASSERT(released && "Strong release of multiple locks failed.");
            if (!released) {
              ERROR("strongReleaseMultipleLocks() failed! Aborting.");
            }
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* lock_manager_node = parent->addChild("Array Lock Manager");
          lock_manager_node->updateSize(_size * sizeof(OwnerID));
        }

        // ! Only for testing
        bool checkNoneLocked() {
            CAtomic<uint8_t> none_locked(uint8_t(true));
            tbb::parallel_for(LockedID(0), _size, [&](const LockedID& id){
                if (isLocked(id)) none_locked.fetch_and(uint8_t(false));
            });
            return bool(none_locked.load(std::memory_order_relaxed));
        }

    private:

        friend class mt_kahypar::AsyncNodeTracker;

        void init() {
            const OwnerID initializer = _invalid_owner_id;
            for ( LockedID i = 0; i < _size; ++i ) {
                _v[i].store(initializer, std::memory_order_acq_rel);
            }
        }

        OwnerID owner(LockedID lockedID) const {
            ASSERT(lockedID < _size);
            return _v[lockedID].load(std::memory_order_acquire);
        }

        const LockedID _size;
        std::unique_ptr<UnderlyingType[]> _v;
        OwnerID _invalid_owner_id;

    };

    using GroupLockManager = ArrayLockManager<HypernodeID, ContractionGroupID>;

}
}// namespace mt_kahypar::ds





#endif //KAHYPAR_ARRAY_LOCK_MANAGER_H
