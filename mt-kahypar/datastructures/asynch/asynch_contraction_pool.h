//
// Created by mlaupichler on 19.04.21.
//

#pragma once

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

        ContractionGroupID _num_groups;
        Array<bool> _active;
        Array<bool> _successors_activated;

    public:
        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit SequentialContractionGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
            : _hierarchy(hierarchy.release()),
            _num_groups(_hierarchy->getNumGroups()),
            _active(_num_groups, false),
            _successors_activated(_num_groups, false) {
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
            ASSERT(isActive(destination));
            ASSERT(!_successors_activated[destination]);
            _active_ids.pop();
            _active[destination] = false;
            return true;
        }

        void activateSuccessors(ContractionGroupID id) {
            ASSERT(!isActive(id));
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s);
            }
            _successors_activated[id] = true;
        }

        bool hasActive() const {
            return !(_active_ids.empty());
        }

        void reactivate(ContractionGroupID id) {
            ASSERT(!isActive(id));
            ASSERT(!_successors_activated[id]);
            activate(id);
        }

    private:

        bool isActive(ContractionGroupID id) const {
            ASSERT(id < _num_groups);
            return _active[id];
        }

        void activate(ContractionGroupID id) {
            ASSERT(!isActive(id));
            _active_ids.push(id);
            _active[id] = true;
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
