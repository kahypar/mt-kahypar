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
            : ILockManager<LockedID, OwnerID>(invalidOwnerID), _size(size), _v(std::make_unique<UnderlyingType[]>(size)), _num_locked(0) {
            init();
        }

        ArrayLockManager(const ArrayLockManager&) = delete;
        ArrayLockManager& operator= (const ArrayLockManager&) = delete;

        ArrayLockManager(ArrayLockManager&&)  noexcept = default;
        ArrayLockManager& operator= (ArrayLockManager&&)  noexcept = default;

        ~ArrayLockManager() = default;

        bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) override {
            OwnerID expected = owner(lockedID);
            OwnerID desired = ownerID;
            if (expected != _invalid_owner_id) {
                // If the lock is already held by someone, return false
                return false;
            }
            // If the lock is currently not held by anyone, attempt to acquire it; returns true if successful or false if not
            bool acquired = _v[lockedID].compare_exchange_strong(expected, desired);
            if (acquired) _num_locked.add_fetch(1, std::memory_order_relaxed);
            return acquired;
        };

        bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) override {
            ASSERT(lockedID < _size);
            OwnerID expected = ownerID;
            OwnerID desired = _invalid_owner_id;
            bool released = _v[lockedID].compare_exchange_strong(expected, desired);
            if (released) _num_locked.sub_fetch(1, std::memory_order_relaxed);
            return released;
        };

        bool isLocked(LockedID lockedID) const override {
            return owner(lockedID) != _invalid_owner_id;
        };

        bool isHeldBy(LockedID lockedID, OwnerID ownerID) const override {
            return owner(lockedID) == ownerID;
        }

        int numLocked() const override {
            return _num_locked.load(std::memory_order_acq_rel);
        }

    private:

        void init() {
            const OwnerID initializer = _invalid_owner_id;
            for ( LockedID i = 0; i < _size; ++i ) {
                _v[i].store(initializer, std::memory_order_acq_rel);
            }
        }

        OwnerID owner(LockedID lockedID) const {
            ASSERT(lockedID < _size);
            return _v[lockedID].load(std::memory_order_acq_rel);
        };

        const LockedID _size;
        std::unique_ptr<UnderlyingType[]> _v;
        CAtomic<int> _num_locked;

        using Base::_invalid_owner_id;

    };

} // namespace mt_kahypar::ds





#endif //KAHYPAR_ARRAY_LOCK_MANAGER_H
