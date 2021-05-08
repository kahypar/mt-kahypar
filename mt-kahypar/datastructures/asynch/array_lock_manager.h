//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_ARRAY_LOCK_MANAGER_H
#define KAHYPAR_ARRAY_LOCK_MANAGER_H


#include "i_lock_manager.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar::ds {

    /// Class template for a ILockManager based on an array of atomic wrappers for the given OwnerID type using the
    /// given LockedID type as indices.
    template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
    class ArrayLockManager : public ILockManager<LockedID, OwnerID> {
        using Base = ILockManager<LockedID, OwnerID>;
        friend Base;

        using UnderlyingType = CAtomic<OwnerID>;

    public:

        /// Constructs a new ArrayLockManager of given size with a given distinct OwnerID to be considered invalid.
        /// \param size the size of the array of locks. All LockedIDs used with an ArrayLockManager instantiated with
        /// this constructor are expected to be in [0,size).
        /// \param invalidOwnerID a distinct OwnerID that an ArrayLockManager instantiated with this constructor
        /// will consider an unassigned or invalid owner.
        ArrayLockManager(LockedID size, OwnerID invalidOwnerID)
            : ILockManager<LockedID, OwnerID>(invalidOwnerID), _size(size), _v(std::make_unique<UnderlyingType[]>(size)) {
            init();
        }

        ArrayLockManager(const ArrayLockManager&) = delete;
        ArrayLockManager& operator= (const ArrayLockManager&) = delete;

        ArrayLockManager(ArrayLockManager&&)  noexcept = default;
        ArrayLockManager& operator= (ArrayLockManager&&)  noexcept = default;

        ~ArrayLockManager() = default;

        bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) override {
            ASSERT(lockedID < _size);
            OwnerID expected = _v[lockedID].load(std::memory_order_relaxed);
            OwnerID desired = ownerID;
            if (expected != _invalid_owner_id) {
                // If the lock is already held by someone, return false
                return false;
            }
            // If the lock is currently not held by anyone, attempt to acquire it; returns true if successful or false if not
            return _v[lockedID].compare_exchange_strong(expected, desired);
        };

        bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) override {
            ASSERT(lockedID < _size);
            OwnerID expected = ownerID;
            OwnerID desired = _invalid_owner_id;
            return _v[lockedID].compare_exchange_strong(expected, desired);
        };

        bool isLocked(LockedID lockedID) const override {
            return owner(lockedID) != _invalid_owner_id;
        };

        bool isHeldBy(LockedID lockedID, OwnerID ownerID) const override {
            return owner(lockedID) == ownerID;
        }

    private:

        void init() {
            const OwnerID initializer = _invalid_owner_id;
            for ( LockedID i = 0; i < _size; ++i ) {
                _v[i].store(initializer, std::memory_order_relaxed);
            }
        }

        OwnerID owner(LockedID lockedID) const {
            ASSERT(lockedID < _size);
            return _v[lockedID].load(std::memory_order_relaxed);
        };

        const LockedID _size;
        std::unique_ptr<UnderlyingType[]> _v;

        using Base::_invalid_owner_id;

    };

} // namespace mt_kahypar::ds





#endif //KAHYPAR_ARRAY_LOCK_MANAGER_H
