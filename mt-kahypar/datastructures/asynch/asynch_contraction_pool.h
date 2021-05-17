//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include <tbb/concurrent_queue.h>
#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/asynch/asynch_common.h"

namespace mt_kahypar::ds
{

    template<typename GroupHierarchy>
    class SequentialContractionGroupPool {

    private:

        std::unique_ptr<GroupHierarchy> _hierarchy;

        std::queue<ContractionGroupID> _active_ids;

    public:
        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit SequentialContractionGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
            : _hierarchy(hierarchy.release()) {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        uint32_t getNumActive() const {
            return _active_ids.size();
        }

        uint32_t getNumTotal() const {
            return _hierarchy->getNumGroups();
        }

        size_t getVersion() const {
            return _hierarchy->getVersion();
        }

        const ContractionGroup &group(ContractionGroupID id) const {
            return _hierarchy->group(id);
        }

        bool pickAnyActiveID(ContractionGroupID& destination) {
            ASSERT(!_active_ids.empty());
            destination = _active_ids.front();
            _active_ids.pop();
            return true;
        }

        void activateSuccessors(ContractionGroupID id) {
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s);
            }
        }

        bool hasActive() const {
            return !(_active_ids.empty());
        }

        void reactivate(ContractionGroupID id) {

            activate(id);
        }

    private:

        void activate(ContractionGroupID id) {
            _active_ids.push(id);
        }

    };

    template<typename GroupHierarchy>
    class ConcurrentQueueGroupPool {

    private:

        std::unique_ptr<GroupHierarchy> _hierarchy;

        tbb::concurrent_queue<ContractionGroupID> _active_ids;

    public:

        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit ConcurrentQueueGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
                : _hierarchy(hierarchy.release()) {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        uint32_t getNumActive() const {
            return _active_ids.unsafe_size();
        }

        uint32_t getNumTotal() const {
            return _hierarchy->getNumGroups();
        }

        size_t getVersion() const {
            return _hierarchy->getVersion();
        }

        const ContractionGroup &group(ContractionGroupID id) const {
            return _hierarchy->group(id);
        }

        bool pickAnyActiveID(ContractionGroupID& destination) {
            bool picked = _active_ids.try_pop(destination);
            return picked;
        }

        void activateSuccessors(ContractionGroupID id) {

            // todo mlaupichler are iterators like successors thread safe on parallel::scalable_vector?
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s);
            }
        }

        bool hasActive() const {
            return !(_active_ids.empty());
        }

        void reactivate(ContractionGroupID id) {
            activate(id);
        }

        // todo mlaupichler: This is used to entirely circumvent the actual pool for parallel_do (which manages
        //  the queue itself via the feeder). If parallel_do is supposed to be used long term, then there should not be
        //  a pool anymore at all and only hierarchies should be calculated
        GroupHierarchy& hierarchy() {
            return *_hierarchy.get();
        }

    private:

        void activate(ContractionGroupID id) {
            _active_ids.push(id);
        }

    };

    using TreeGroupPool = ConcurrentQueueGroupPool<UncontractionGroupTree>;
    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H
