//
// Created by mlaupichler on 19.04.21.
//

#pragma once

#include <queue>
#include <tbb/concurrent_queue.h>
#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/async/async_common.h"

namespace mt_kahypar::ds
{
    using DoParallelForAllGroupsFunction = std::function<void (const ContractionGroup&)>;
    using DoParallelForAllGroupIDsFunction = std::function<void (const ContractionGroupID&)>;

    template<typename GroupHierarchy>
    class SequentialGroupPool {

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
        explicit SequentialGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
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

        const size_t& depth(ContractionGroupID id) const {
            return _hierarchy->depth(id);
        }

//        const size_t& node_depth(HypernodeID id) const {
//            return _hierarchy->node_depth(id);
//        }

        bool tryToPickActiveID(ContractionGroupID& destination) {
            ASSERT(!_active_ids.empty());
            destination = _active_ids.front();
            ASSERT(isActive(destination));
            ASSERT(!_successors_activated[destination]);
            _active_ids.pop();
            _active[destination] = false;
            return true;
        }

        void pickActiveID(ContractionGroupID& destination) {
            tryToPickActiveID(destination);
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

        void doParallelForAllGroups(const DoParallelForAllGroupsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(group(id));
                }
            });
        }

        void doParallelForAllGroupIDs(const DoParallelForAllGroupIDsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(id);
                }
            });
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

        BlockedGroupIDIterator all() const {
            return _hierarchy->all();
        }

    };

    template<typename GroupHierarchy>
    class ConcurrentQueueGroupPool {

    public:

        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit ConcurrentQueueGroupPool(std::unique_ptr<GroupHierarchy> hierarchy)
                : _hierarchy(hierarchy.release()),
//                _queue_lock(),
//                _active_ids(_cmp)
                  _active_ids()
                {
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        uint32_t getNumActive() {
            return unsafeNumActive();
        }

        uint32_t getNumTotal() {
            return _hierarchy->getNumGroups();
        }

        size_t getVersion() const {
            return _hierarchy->getVersion();
        }

        const ContractionGroup &group(ContractionGroupID id) const {
            return _hierarchy->group(id);
        }

        const size_t& depth(ContractionGroupID id) const {
            return _hierarchy->depth(id);
        }

        bool tryToPickActiveID(ContractionGroupID& destination) {
//            ASSERT(hasActive());
            bool picked = tryPopActive(destination);
            return picked;
        }

        void pickActiveID(ContractionGroupID& destination) {
            while (!tryToPickActiveID(destination)) {}
        }

        /// Convenience method to call activate() on all successors
        void activateAllSuccessors(ContractionGroupID id) {
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s);
            }
        }

        ContractionGroupIDIteratorRange successors(ContractionGroupID id) {
            return _hierarchy->successors(id);
        }

        ContractionGroupID numSuccessors(ContractionGroupID id) {
            return _hierarchy->numSuccessors(id);
        }

        bool hasActive() {
            return !unsafeEmpty();
        }

        void activate(ContractionGroupID id) {
            insertActive(id);
        }

//        void reactivate(ContractionGroupID id) {
//            activate(id);
//        }

        void doParallelForAllGroups(const DoParallelForAllGroupsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(group(id));
                }
            });
        }

        void doParallelForAllGroupIDs(const DoParallelForAllGroupIDsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(id);
                }
            });
        }

        size_t getTotalNumUncontractions() {
            size_t num = 0;
            auto range = all();
            for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                num += group(id).size();
            }
            return num;
        }

    private:

        BlockedGroupIDIterator all() const {
            return _hierarchy->all();
        }

        bool tryInsertActive(ContractionGroupID id) {
//            bool locked = _queue_lock.tryLock();
//            if (locked) {
//                _active_ids.push(id);
//                _queue_lock.unlock();
//            }
//            return locked;
              _active_ids.push(id);
              return true;
        }


        bool tryPopActive(ContractionGroupID& destination) {
//            bool locked = _queue_lock.tryLock();
//            if (locked) {
//                destination = _active_ids.top();
//                _active_ids.pop();
//                _queue_lock.unlock();
//            }
//            return locked;
              return _active_ids.try_pop(destination);
        }

        void insertActive(ContractionGroupID id) {
//            _queue_lock.lock();
            _active_ids.push(id);
//            _queue_lock.unlock();
        }

        void popActive(ContractionGroupID& destination) {
//            _queue_lock.lock();
//            destination = _active_ids.top();
//            _active_ids.pop();
//            _queue_lock.unlock();
            while (!_active_ids.try_pop(destination)) {/* continue trying to pop */}
        }

        ContractionGroupID unsafeNumActive() {
//            _queue_lock.lock();
            ContractionGroupID num_active = _active_ids.unsafe_size();
//            _queue_lock.unlock();
            return num_active;
        }

        bool unsafeEmpty() {
//            _queue_lock.lock();
            bool empty = _active_ids.empty();
//            _queue_lock.unlock();
            return empty;
        }

//        using DepthCompare = std::function<bool (ContractionGroupID, ContractionGroupID)>;

        // Used to compare groups by their depth in the GroupHierarchy. The comparator for std::priority queue has to be
        // true if id1 < id2 in the queue but std::priority_queue emits the largest elements first, so here we reverse
        // the ordering in order to emit the smallest depth group first.
//        DepthCompare _cmp = [&](ContractionGroupID id1, ContractionGroupID id2) {
//            return _hierarchy->depth(id1) > _hierarchy->depth(id2);
//        };

        std::unique_ptr<GroupHierarchy> _hierarchy;

//        SpinLock _queue_lock;
//        std::priority_queue<ContractionGroupID, std::vector<ContractionGroupID>, DepthCompare> _active_ids;
        tbb::concurrent_queue<ContractionGroupID> _active_ids;

    };

    using TreeGroupPool = ConcurrentQueueGroupPool<UncontractionGroupTree>;
    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar
