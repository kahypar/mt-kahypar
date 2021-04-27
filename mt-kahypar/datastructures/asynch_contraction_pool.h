//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_queue.h>
#include "gmock/gmock.h"

#include <utility>
#include "hypergraph_common.h"
#include "uncontraction_group_tree.h"
#include <list>

namespace mt_kahypar::ds
{
    //Forward
    class IContractionGroupPool;

    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<IContractionGroupPool>>;

    /// Pure virtual interface to define functions for a pool that manages active/inactive uncontraction groups.
    class IContractionGroupPool {
    public:
        virtual uint32_t getNumActive() const = 0;
        virtual bool hasActive() const = 0;
        virtual const ContractionGroup& group(const ContractionGroupID id) const = 0;
        virtual ContractionGroupID pickAnyActiveID() = 0;
        virtual void activateSuccessors(ContractionGroupID id) = 0;

        template<typename F> void doForAllGroupIDsInParallel(const F& f) {
            tbb::parallel_for(all(), [&](const ContractionGroupID& id) {
                f(id);
            });
        };

        template<typename F> void doForAllGroupsInParallel(const F& f) {
            tbb::parallel_for(all(), [&](const ContractionGroupID& id) {
                f(group(id));
            });
        };

    protected:
        virtual BlockedGroupIDIterator all() const = 0;
    };

    class SequentialContractionGroupPool : public IContractionGroupPool {

    private:
        parallel::scalable_vector<ContractionGroupID> _active;
        std::unique_ptr<IUncontractionGroupHierarchy> _hierarchy;

    public:
        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit SequentialContractionGroupPool(std::unique_ptr<IUncontractionGroupHierarchy> hierarchy) : _hierarchy(hierarchy.release()) {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                _active.push_back(r);
            }
        }

        uint32_t getNumActive() const override;

        const ContractionGroup &group(ContractionGroupID id) const override;

        ContractionGroupID pickAnyActiveID() override;

        void activateSuccessors(ContractionGroupID id) override;

        bool hasActive() const override;

    protected:
        BlockedGroupIDIterator all() const override;

    };

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
