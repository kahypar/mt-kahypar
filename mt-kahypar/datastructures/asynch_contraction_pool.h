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
    /// Pure virtual interface to define functions for a pool that manages active/inactive uncontraction groups.
    class IContractionGroupPool {
    public:
        virtual uint32_t getNumActive() const = 0;
        virtual bool hasActive() const = 0;
        virtual const ContractionGroup& group(ContractionGroupID id) const = 0;
        virtual ContractionGroupID pickAnyActiveID() = 0;
        virtual void activateSuccessors(ContractionGroupID id) = 0;
    };

    class SequentialContractionGroupPool : public IContractionGroupPool {

    private:
        std::vector<ContractionGroupID> _active;
        IUncontractionGroupHierarchy& _hierarchy;

    public:

        explicit SequentialContractionGroupPool(IUncontractionGroupHierarchy& hierarchy) : _hierarchy(hierarchy) {
            auto roots = _hierarchy.roots();
            for(auto r: roots) {
                _active.push_back(r);
            }
        }

        uint32_t getNumActive() const override;

        const ContractionGroup &group(ContractionGroupID id) const override;

        ContractionGroupID pickAnyActiveID() override;

        void activateSuccessors(ContractionGroupID id) override;

        bool hasActive() const override;

    };

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
