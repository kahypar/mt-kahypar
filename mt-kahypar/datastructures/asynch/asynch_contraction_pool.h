//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/asynch/asynch_common.h"
#include "mt-kahypar/datastructures/asynch/i_contraction_group_pool.h"

namespace mt_kahypar::ds
{

    class SequentialContractionGroupPool : public IContractionGroupPool {

    private:

        std::unique_ptr<IUncontractionGroupHierarchy> _hierarchy;
        Array<bool> _active;
        Array<bool> _successors_activated;

        parallel::scalable_vector<ContractionGroupID> _active_ids;

    public:
        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit SequentialContractionGroupPool(std::unique_ptr<IUncontractionGroupHierarchy> hierarchy)
            : _hierarchy(hierarchy.release()),
            _active(_hierarchy->getNumGroups(),false),
            _successors_activated(_hierarchy->getNumGroups(),false) {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        uint32_t getNumActive() const override;

        uint32_t getNumTotal() const override;

        size_t getVersion() const override;

        const ContractionGroup &group(ContractionGroupID id) const override;

        ContractionGroupID pickAnyActiveID() override;

        void activateSuccessors(ContractionGroupID id) override;

        bool hasActive() const override;

        void reactivate(ContractionGroupID id) override;

    protected:
        BlockedGroupIDIterator all() const override;

    private:

        bool isActive(ContractionGroupID id) const;
        void activate(ContractionGroupID id);

    };

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
