//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/asynch/asynch_common.h"
#include "mt-kahypar/datastructures/asynch/contraction_group_pool.h"

namespace mt_kahypar::ds
{

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

        uint32_t getNumTotal() const override;

        const ContractionGroup &group(ContractionGroupID id) const override;

        ContractionGroupID pickAnyActiveID() override;

        void activateSuccessors(ContractionGroupID id) override;

        bool hasActive() const override;

    protected:
        BlockedGroupIDIterator all() const override;

    };

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
