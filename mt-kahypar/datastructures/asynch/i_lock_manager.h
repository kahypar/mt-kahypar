//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_I_LOCK_MANAGER_H
#define KAHYPAR_I_LOCK_MANAGER_H

/// Pure virtual interface for managing locks. Locked objects are marked by LockedIDs and mapped to OwnerIDs describing
/// who holds the lock. Requires LockedID and OwnerID to be integral numbers.
template<typename LockedID = Mandatory, typename OwnerID = Mandatory>
class ILockManager {

public:

    /// Attempt to acquire a lock for the id lockedID with the owner ownerID. Will return true if afterwards the lock is
    /// owned by the owner with ownerID and false if not.
    virtual bool tryToAcquireLock(LockedID lockedID, OwnerID ownerID) = 0;

    /// Release the lock for lockedID held by ownerID. Will return true if the lock has been successfully released or
    /// false if the given ownerID did not hold the lock prior to the call.
    virtual bool tryToReleaseLock(LockedID lockedID, OwnerID ownerID) = 0;

    /// Returns the OwnerID of the current owner of lockedID or the invalid OwnerID if it is not locked.
    virtual OwnerID owner(LockedID lockedID) = 0;

    /// Returns true if any owner holds the lock for lockedID.
    virtual bool isLocked(LockedID lockedID) = 0;

    OwnerID operator[] (const LockedID lockedID) const {
        return owner(lockedID);
    };

protected:

    ILockManager(OwnerID invalidOwnerID) :
        _invalid_owner_id(invalidOwnerID) {};

    ~ILockManager() = default;

    OwnerID _invalid_owner_id;

};


#endif //KAHYPAR_I_LOCK_MANAGER_H
