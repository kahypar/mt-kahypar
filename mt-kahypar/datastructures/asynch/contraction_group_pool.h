//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_CONTRACTION_GROUP_POOL_H
#define KAHYPAR_CONTRACTION_GROUP_POOL_H

#include "gmock/gmock.h"

#include "asynch_common.h"

namespace mt_kahypar::ds {

    //Forward
    class IContractionGroupPool;

    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<IContractionGroupPool>>;

    /// Pure virtual interface to define functions for a pool that manages active/inactive uncontraction groups.
    class IContractionGroupPool {
    public:
        virtual uint32_t getNumActive() const = 0;
        virtual uint32_t getNumTotal() const = 0;
        virtual size_t getVersion() const = 0;
        virtual bool hasActive() const = 0;
        virtual const ContractionGroup& group(const ContractionGroupID id) const = 0;
        virtual ContractionGroupID pickAnyActiveID() = 0;
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

#endif //KAHYPAR_CONTRACTION_GROUP_POOL_H
