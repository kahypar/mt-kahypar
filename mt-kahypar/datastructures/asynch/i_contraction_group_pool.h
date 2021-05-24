//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_I_CONTRACTION_GROUP_POOL_H
#define KAHYPAR_I_CONTRACTION_GROUP_POOL_H

#include "gmock/gmock.h"

#include "asynch_common.h"

namespace mt_kahypar::ds {

    //Forward
    class IContractionGroupPool;


    /// Pure virtual interface to define functions for a pool that manages active/inactive uncontraction groups.
    class IContractionGroupPool {
    public:

        /// Get the number of currently active contraction groups in the pool.
        virtual uint32_t getNumActive() const = 0;

        /// Get the total number of contraction groups that may be activated in the pool through activateAllSuccessors
        /// (i.e. the size of the reservoir of groups).
        virtual uint32_t getNumTotal() const = 0;

        /// Get the version of the underlying hypergraph that the Contractions in the ContractionGroups in this pool
        /// have been contracted in.
        virtual size_t getVersion() const = 0;

        /// Returns whether there are any active groups in the pool.
        virtual bool hasActive() const = 0;

        /// Retrieves the group associated with a specific group id.
        virtual const ContractionGroup& group(const ContractionGroupID id) const = 0;

        /// Picks any active group from the pool and deactivates it. This does not at all restrict the order in which
        /// groups are picked. (More specifically, it is not necessarily random).
        virtual ContractionGroupID pickAnyActiveID() = 0;

        /// Reactivates the group with the given id. A group can only be reactivated if it is not active and as long
        /// as its successors have not been activated yet.
        virtual void reactivate(ContractionGroupID id) = 0;

        /// Activates all successors of the group with the given id. The definition of a successor depends on the
        /// underlying hierarchy of groups and is implementation specific.
        virtual void activateSuccessors(ContractionGroupID id) = 0;


        template<typename F> void doForAllGroupIDsInParallel(const F& f) {
            tbb::parallel_for(all(), [&](const tbb::blocked_range<ContractionGroupID>& r) {
                for (ContractionGroupID i = r.begin(); i != r.end(); ++i) {
                    f(i);
                }
            });
        };


        template<typename F>
        void doForAllGroupsInParallel(const F& f) {
            static_cast<const IContractionGroupPool&>(*this).template doForAllGroupsInParallel(f);
        };

        template<typename F>
        void doForAllGroupsInParallel(const F& f) const {
            tbb::parallel_for(all(), [&](const tbb::blocked_range<ContractionGroupID>& r) {
                for (ContractionGroupID i = r.begin(); i != r.end(); ++i) {
                    f(group(i));
                }
            });
        };

    protected:
        virtual BlockedGroupIDIterator all() const = 0;
    };


}

#endif //KAHYPAR_I_CONTRACTION_GROUP_POOL_H
