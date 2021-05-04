//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_ARRAY_LOCK_MANAGER_H
#define KAHYPAR_ARRAY_LOCK_MANAGER_H


#include "i_lock_manager.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar::ds {

    template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
    class ArrayLockManager : public ILockManager<LockedID, OwnerID> {
        using Base = ILockManager<LockedID, OwnerID>;
        friend Base;

        using UnderlyingType = CAtomic<OwnerID>;

    public:

        ArrayLockManager(LockedID size, OwnerID invalidOwnerID)
            : ILockManager<LockedID, OwnerID>(invalidOwnerID), _size(size), _v(std::make_unique<UnderlyingType[]>(size)) {
            init();
        }

        bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) override {
            OwnerID expected = _v[lockedID].load(std::memory_order_relaxed);
            OwnerID desired = ownerID;
            if ( expected != _invalid_owner_id && _v[lockedID].compare_exchange_strong(expected, desired) ) {
                // Value was successfully set from invalid OwnerID to the given OwnerID
                return true;
            } else {
                // Either expected == _invalid_owner which means it is already locked
                // or compare_exchange_strong failed, which means that another owner holds the lock on lockedID.
                return false;
            }
        };

        bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) override {
            OwnerID expected = ownerID;
            OwnerID desired = _invalid_owner_id;
            return _v[lockedID].compare_exchange_strong(expected, desired);
        };

        OwnerID owner(LockedID lockedID) override {
            return _v[lockedID].load(std::memory_order_relaxed);
        };

        bool isLocked(LockedID lockedID) override {
            return _v[lockedID].load(std::memory_order_relaxed) == _invalid_owner_id;
        };

    private:

        void init() {
            for ( LockedID i = 0; i < _size; ++i ) {
                _v[i].store(_invalid_owner_id, std::memory_order_relaxed);
            }
        }

        LockedID _size;
        std::unique_ptr<UnderlyingType[]> _v;

        using Base::_invalid_owner_id;

    };

} // namespace mt_kahypar::ds





#endif //KAHYPAR_ARRAY_LOCK_MANAGER_H
