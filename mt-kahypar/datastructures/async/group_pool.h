//
// Created by mlaupichler on 19.04.21.
//

#pragma once

//#include <queue>
//#include <tbb/concurrent_queue.h>
#include <mt-kahypar/partition/context.h>
#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/async/async_common.h"
#include "depth_priority_queue.h"
#include "node_region_comparator.h"
#include "mt-kahypar/definitions.h"

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

    template<typename GroupHierarchy = Mandatory,
        typename Hypergraph = Mandatory>
    class ConcurrentQueueGroupPool {

    public:

        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit ConcurrentQueueGroupPool(std::unique_ptr<GroupHierarchy> hierarchy, const Context& context)
                : _hierarchy(hierarchy.release()),
//                _cmp([&](ContractionGroupID id1, ContractionGroupID id2) {
//                    return _hierarchy->depth(id1) > _hierarchy->depth(id2);
//                  }),
//                _queue_lock(),
//                _active_ids(_cmp)
                  _active_ids(_hierarchy->getNumberOfDepths(), std::move(_hierarchy->getNumberOfGroupsPerDepth())),
                  _node_region_comparator(nullptr),
                  _region_similarity_retries(context.uncoarsening.node_region_similarity_retries),
                  _region_similarity_threshold(context.uncoarsening.node_region_similarity_threshold)
                {
            ASSERT(_hierarchy);
            auto roots = _hierarchy->roots();
            for(auto r: roots) {
                activate(r);
            }
        }

        ConcurrentQueueGroupPool(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool& operator=(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool(ConcurrentQueueGroupPool&& other) = delete;
        ConcurrentQueueGroupPool& operator=(ConcurrentQueueGroupPool&& other) = delete;

        ~ConcurrentQueueGroupPool() = default;

        uint32_t getNumActive() {
            return unsafeNumActive();
        }

        const GroupHierarchy* getPtrToHierarchyForQueries() {
          ASSERT(_hierarchy);
          return _hierarchy.get();
        }

        void setNodeRegionComparator(const NodeRegionComparator<Hypergraph>* comparator) {
          ASSERT(comparator);
          _node_region_comparator = comparator;
        }

        // ===== Hierarchy forwards =====

        uint32_t getNumTotal() {
          ASSERT(_hierarchy);
            return _hierarchy->getNumGroups();
        }

        HypernodeID getTotalNumUncontractions() const {
          ASSERT(_hierarchy);
          return _hierarchy->getNumContainedContractions();
        }

        size_t getVersion() const {
          ASSERT(_hierarchy);
            return _hierarchy->getVersion();
        }

        const ContractionGroup &group(const ContractionGroupID id) const {
          ASSERT(_hierarchy);
            return _hierarchy->group(id);
        }

        const HypernodeID& depth(const ContractionGroupID id) const {
          ASSERT(_hierarchy);
            return _hierarchy->depth(id);
        }

        const HypernodeID getNumberOfUncontractionsInDepth(const HypernodeID depth) const {
          ASSERT(_hierarchy);
          return _hierarchy->getNumberOfUncontractionsInDepth(depth);
        }

        bool isLastContractionGroupOfRepresentative(const ContractionGroupID groupID) const {
          ASSERT(_hierarchy);
          return _hierarchy->isLastContractionGroupOfRepresentative(groupID);
        }

        // ! Returns whether a node is initially stable for this version, i.e. whether it is not the representative in
        // ! any uncontraction in the underlying hierarchy
        bool isInitiallyStableNode(const HypernodeID hn) const {
          ASSERT(_hierarchy);
          return _hierarchy->isInitiallyStableNode(hn);
        }

        ContractionGroupIDIteratorRange successors(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          return _hierarchy->successors(id);
        }

        ContractionGroupID numSuccessors(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          return _hierarchy->numSuccessors(id);
        }

        // ===== Activation/Deactivation of Hypernodes =====

        bool tryToPickActiveID(ContractionGroupID& destination) {
            bool picked = tryPopActive(destination);
            return picked;
        }

        void pickActiveID(ContractionGroupID& destination) {

            destination = invalidGroupID;
            if (!_node_region_comparator) {
              while (!tryToPickActiveID(destination)) {}
              ASSERT(!_hierarchy->isInitiallyStableNode(group(destination).getRepresentative()));
            } else {
              std::vector<ContractionGroupID> pickedIDs(_region_similarity_retries, invalidGroupID);
              double cur_min_similarity = 1.0;
              ContractionGroupID cur_min_idx = 0;
              for (uint32_t i = 0; i < _region_similarity_retries; ++i) {
                ContractionGroupID pickedID = invalidGroupID;
                bool picked = tryToPickActiveID(pickedID);
                if (!picked) {
                  // If none found in PQ and none have been previously, wait for an element and return it immediately
                  if (i == 0) {
                    while(!tryToPickActiveID(destination)) {}
                    _last_picked_ets.local() = destination;
                    return;
                  }
                  // If no more found in PQ but previously some ID was found, reinsert all but best and return best
                  destination = pickedIDs[cur_min_idx];
                  for (uint32_t j = 0; j < i; ++j) {
                    if (j != cur_min_idx) {
                      insertActive(pickedIDs[j]);
                    }
                  }
                  _last_picked_ets.local() = destination;
                  return;
                }
                ASSERT(pickedID != invalidGroupID);
                double sim = maxRegionSimilarityToLastPicked(group(pickedID).getRepresentative());
                if (sim <= _region_similarity_threshold) {
                  // If good candidate found, return it right away and reinsert others that were picked
                  destination = pickedID;
                  for (uint32_t j = 0; j < i; ++j) {
                      insertActive(pickedIDs[j]);
                  }
                  _last_picked_ets.local() = destination;
                  return;
                } else {
                  // If candidate is not good, insert it into pickedIDs and possibly update min
                  pickedIDs[i] = pickedID;
                  if (sim < cur_min_similarity) {
                    cur_min_similarity = sim;
                    cur_min_idx = i;
                  }
                }
              }
              // If tried as often as possible, reinsert all but best and return best
              destination = pickedIDs[cur_min_idx];
              for (uint32_t j = 0; j < _region_similarity_retries; ++j) {
                if (j != cur_min_idx) {
                  insertActive(pickedIDs[j]);
                }
              }
            }
            ASSERT(destination != invalidGroupID);
            _last_picked_ets.local() = destination;
        }

        /// Convenience method to call activate() on all successors
        void activateAllSuccessors(const ContractionGroupID id) {
          ASSERT(_hierarchy);
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s);
            }
        }

        bool hasActive() {
            return !unsafeEmpty();
        }

        void activate(const ContractionGroupID id) {
            insertActive(id);
        }

        // ! Mark a given id as finished, i.e. mark that it will not be reinserted into the pool again.
        void finalize(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          _active_ids.increment_finished(_hierarchy->depth(id));
        }

        // ! Mark a given id as finished, i.e. mark that it will not be reinserted into the pool again.
        // ! Sets completed_depth_of_id to true if this group is the last to be uncontracted in its depth.
        void finalize(const ContractionGroupID id, bool& completed_depth_of_id) {
          ASSERT(_hierarchy);
          _active_ids.increment_finished(_hierarchy->depth(id), completed_depth_of_id);
        }

        // ! Only for testing/debugging
        bool checkAllFinalized() const {
          return _active_ids.allDepthsCompleted();
//            return true;
        }

        // ===== Parallel Iteration Convenience Methods =====

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

        BlockedGroupIDIterator all() const {
          ASSERT(_hierarchy);
          return _hierarchy->all();
        }

//        bool tryInsertActive(ContractionGroupID id) {
//          ASSERT(_hierarchy);
////          _active_ids.push(id, _hierarchy->depth(id));
////          return true;
//          bool locked = _queue_lock.tryLock();
//          if (!locked) return false;
//          _active_ids.push(id);
//          _queue_lock.unlock();
//          return true;
//        }


        bool tryPopActive(ContractionGroupID& destination) {
//                _queue_lock.lock();
//                destination = _active_ids.top();
//                _active_ids.pop();
//                _queue_lock.unlock();
//                ASSERT(destination != invalidGroupID);
//                return true;
              return _active_ids.try_pop(destination);
        }

        void insertActive(ContractionGroupID id) {
            ASSERT(_hierarchy);
            _active_ids.push(id, _hierarchy->depth(id));
//            _queue_lock.lock();
//            _active_ids.push(id);
//            _queue_lock.unlock();
        }

//        void popActive(ContractionGroupID& destination) {
//            while (!_active_ids.try_pop(destination)) {/* continue trying to pop */}
//        }

        ContractionGroupID unsafeNumActive() {
            ContractionGroupID num_active = _active_ids.unsafe_size();
            return num_active;
//            return _active_ids.size();
        }

        bool unsafeEmpty() {
            bool empty = _active_ids.unsafe_empty();
//            bool empty = _active_ids.empty();
            return empty;
        }

        double maxRegionSimilarityToLastPicked(const HypernodeID hn) {
          ASSERT(_node_region_comparator);
          double cur_max = 0.0;
          for (const ContractionGroupID last_picked_id : _last_picked_ets) {
            double sim = _node_region_comparator->regionSimilarity(hn, group(last_picked_id).getRepresentative());
            if (sim > cur_max) {
              cur_max = sim;
            }
          }
          return cur_max;
        }

//        using DepthCompare = std::function<bool (ContractionGroupID, ContractionGroupID)>;

        // Used to compare groups by their depth in the GroupHierarchy. The comparator for std::priority queue has to be
        // true if id1 < id2 in the queue but std::priority_queue emits the largest elements first, so here we reverse
        // the ordering in order to emit the smallest depth group first.
//        DepthCompare _cmp; // = [&](ContractionGroupID id1, ContractionGroupID id2) {
//            return _hierarchy->depth(id1) > _hierarchy->depth(id2);
//        };

        std::unique_ptr<GroupHierarchy> _hierarchy;

//        SpinLock _queue_lock;
//        std::priority_queue<ContractionGroupID, std::vector<ContractionGroupID>, DepthCompare> _active_ids;
//        tbb::concurrent_queue<ContractionGroupID> _active_ids;
          DepthPriorityQueue _active_ids;

          const NodeRegionComparator<Hypergraph>* _node_region_comparator;

          tbb::enumerable_thread_specific<ContractionGroupID> _last_picked_ets;

          const size_t _region_similarity_retries;
          const double _region_similarity_threshold;

    };

//    using TreeGroupPool = ConcurrentQueueGroupPool<UncontractionGroupTree, Hypergraph>;
//    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar::ds
