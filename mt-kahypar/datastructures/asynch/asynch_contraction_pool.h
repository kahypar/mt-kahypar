//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/asynch/asynch_common.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace mt_kahypar::ds
{

    template<typename GroupHierarchy>
    class SequentialContractionGroupPool {

    private:

        std::unique_ptr<GroupHierarchy> _hierarchy;
        Array<bool> _active;
        Array<bool> _successors_activated;

//        parallel::scalable_vector<ContractionGroupID> _active_ids;
        std::queue<ContractionGroupID> _active_ids;

    public:
        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit SequentialContractionGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
            : _hierarchy(hierarchy.release()),
            _active(_hierarchy->getNumGroups(),false),
            _successors_activated(_hierarchy->getNumGroups(),false) {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        uint32_t getNumActive() const;

        uint32_t getNumTotal() const;

        size_t getVersion() const;

        const ContractionGroup &group(ContractionGroupID id) const;

        ContractionGroupID pickAnyActiveID();

        void activateSuccessors(ContractionGroupID id);

        bool hasActive() const;

        void reactivate(ContractionGroupID id);

    protected:
        BlockedGroupIDIterator all() const;

    private:

        bool isActive(ContractionGroupID id) const;
        void activate(ContractionGroupID id);



    };

    using TreeGroupPool = SequentialContractionGroupPool<UncontractionGroupTree>;
    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
