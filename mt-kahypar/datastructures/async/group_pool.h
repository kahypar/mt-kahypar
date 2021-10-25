//
// Created by mlaupichler on 19.04.21.
//

#pragma once

#include <queue>
//#include <tbb/concurrent_queue.h>
#include <mt-kahypar/partition/context.h>
#include <multiqueue/int_multiqueue.hpp>
#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/async/async_common.h"
#include "depth_priority_queue.h"
#include "node_region_comparator.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

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
        typename RegionComparator = Mandatory>
    class ConcurrentQueueGroupPool {

        struct StickyMQConfig : multiqueue::configuration::FullBuffering {

            static constexpr unsigned int K = 8;
            static constexpr unsigned int C = 8;

        };

        using Multiqueue = multiqueue::int_multiqueue<uint32_t, ContractionGroupID, StickyMQConfig>;

    public:

        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit ConcurrentQueueGroupPool(std::unique_ptr<GroupHierarchy> hierarchy, const Context &context,
                                          const bool use_multiqueue)
                : _hierarchy(hierarchy.release()),
                  _node_region_comparator(nullptr),
                  _region_similarity_retries(context.uncoarsening.node_region_similarity_retries)
                  , _picked_ids_ets(std::vector<ContractionGroupID>(context.uncoarsening.node_region_similarity_retries, invalidGroupID)),
                  _pick_candidates_ets(std::vector<pick_candidate>(context.uncoarsening.node_region_similarity_retries)),
                  _total_calls_to_pick(0),
                  _calls_to_pick_that_reached_max_retries(0),
                  _calls_to_pick_with_empty_pq(0),
                  _num_accepted_uncontractions(0),
                  _calculate_full_similarities(context.uncoarsening.region_comparison_with_full_similarities),
                  _use_multiqueue(use_multiqueue),
                  _similarity_to_others_threshold(context.uncoarsening.node_region_similarity_threshold)
                {
            ASSERT(_hierarchy);

            // Ugly code duplication because interfaces can't be used (slow vtable lookups)
            if (_use_multiqueue) {
              _mq_active_ids = std::make_unique<Multiqueue>(context.shared_memory.num_threads);
            } else {
              _dpq_active_ids = std::make_unique<DepthPriorityQueue>(_hierarchy->getNumberOfDepths(), std::move(_hierarchy->getNumberOfGroupsPerDepth()));
            }

            auto roots = _hierarchy->roots();
            size_t i = 0;
            for(auto r: roots) {
                activate(r, i++);
                if (i == context.shared_memory.num_threads) {
                  i = 0;
                }
            }
        }

        ConcurrentQueueGroupPool(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool& operator=(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool(ConcurrentQueueGroupPool&& other) = delete;
        ConcurrentQueueGroupPool& operator=(ConcurrentQueueGroupPool&& other) = delete;

        ~ConcurrentQueueGroupPool() = default;

        const GroupHierarchy* getPtrToHierarchyForQueries() {
          ASSERT(_hierarchy);
          return _hierarchy.get();
        }

        void setNodeRegionComparator(RegionComparator* comparator) {
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

        bool tryToPickActiveID(ContractionGroupID& destination, const size_t task_id) {
          if (_calculate_full_similarities) {
            return tryToPickActiveIDWithFullSimilarities(destination, task_id);
          } else {
            return tryToPickActiveIDWithEarlyBreak(destination, task_id);
          }
        }

        void pickActiveID(ContractionGroupID& destination, const size_t task_id) {
          if (_calculate_full_similarities) {
            while (!tryToPickActiveIDWithFullSimilarities(destination, task_id)) {/* keep trying */};
          } else {
            while (!tryToPickActiveIDWithEarlyBreak(destination, task_id)) {/* keep trying */};
          }
        }

        bool tryToPickActiveIDWithFullSimilarities(ContractionGroupID& destination, const size_t task_id, const bool accept_based_on_similarity_to_others = true) {
          _total_calls_to_pick.add_fetch(1, std::memory_order_relaxed);

            destination = invalidGroupID;
            if (!_node_region_comparator || (_region_similarity_retries == 0)) {
              return tryPopActive(destination, task_id);
            } else {
              const size_t num_edges_active_in_this_and_other_tasks = _node_region_comparator->numberOfEdgesActiveInThisTaskAndAnyOtherTask(task_id);
              std::vector<pick_candidate>& candidates = _pick_candidates_ets.local();
              ContractionGroupID best_idx = 0;
              for (uint32_t i = 0; i < _region_similarity_retries; ++i) {
                pick_candidate next_candidate;
                bool picked = tryToGetNewCandidate(next_candidate, task_id, num_edges_active_in_this_and_other_tasks);
                if (!picked) {
                  // If none found in PQ and none have been previously, return false (no elements available)
                  if (i == 0) {
                    _calls_to_pick_with_empty_pq.fetch_add(1, std::memory_order_relaxed);
                    return false;
                  }
                  // If no more found in PQ but previously some ID was found, reinsert all but best and return best
                  destination = candidates[best_idx].group_id;
                  for (uint32_t j = 0; j < i; ++j) {
                    if (j != best_idx) {
                      insertActive(candidates[j].group_id, task_id);
                    }
                  }
                  if (destination == invalidGroupID) {
                    ERROR("Destination invalid after picking best from so far picked when PQ popping failed!" << V(destination));
                  }
                  return true;
                }
                ASSERT(next_candidate.group_id != invalidGroupID);
                if (next_candidate.group_id == invalidGroupID) {
                  ERROR("No correct group has been picked!" << V(next_candidate.group_id));
                }
                if (isGoodCandidate(next_candidate, accept_based_on_similarity_to_others)) {
                  // If good candidate found, return it right away and reinsert others that were picked
                  destination = next_candidate.group_id;
                  for (uint32_t j = 0; j < i; ++j) {
                      insertActive(candidates[j].group_id, task_id);
                  }
                  return true;
                } else {
                  // If candidate is not good, insert it into candidates and possibly update best
                  candidates[i] = next_candidate;
                  if (isBetterThanOtherCandidate(next_candidate, candidates[best_idx])) {
                    best_idx = i;
                  }
                }
              }
              // If tried as often as possible, reinsert all but best and return best
              _calls_to_pick_that_reached_max_retries.fetch_add(1, std::memory_order_relaxed);
              destination = candidates[best_idx].group_id;
              for (uint32_t j = 0; j < _region_similarity_retries; ++j) {
                if (j != best_idx) {
                  insertActive(candidates[j].group_id, task_id);
                }
              }
              ASSERT(destination != invalidGroupID);
              return true;
            }
        }

        bool tryToPickActiveIDWithEarlyBreak(ContractionGroupID& destination, const size_t task_id) {
          _total_calls_to_pick.add_fetch(1, std::memory_order_relaxed);

          destination = invalidGroupID;
          if (!_node_region_comparator || (_region_similarity_retries == 0)) {
            return tryPopActive(destination, task_id);
          } else {
            std::vector<ContractionGroupID>& pickedIDs = _picked_ids_ets.local();
            size_t best_idx = 0;
            double best_part_of_edges_before_fail = 0.0;
            for (uint32_t i = 0; i < _region_similarity_retries; ++i) {
              ContractionGroupID pickedID = invalidGroupID;
              bool picked = tryPopActive(pickedID, task_id);
              if (!picked) {
                if (i == 0) {
                  _calls_to_pick_with_empty_pq.fetch_add(1, std::memory_order_relaxed);
                  return false;
                } else {
                  // If no more found in PQ and no good one was found just return the one that failed at the latest (heuristic) and reinsert the rest
                  destination = pickedIDs[best_idx];
                  for (uint32_t j = 0; j < i; ++j) {
                    if (j != best_idx) {
                      insertActive(pickedIDs[j], task_id);
                    }
                  }
                }
                ASSERT(destination != invalidGroupID);
                return true;
              }
              ASSERT(pickedID != invalidGroupID);
              HypernodeID repr = group(pickedID).getRepresentative();
              double part_of_edges_seen_before_failure;
              if (_node_region_comparator->regionIsNotTooSimilarToActiveNodesWithEarlyBreak(repr, task_id, part_of_edges_seen_before_failure)) {
                // If good candidate found return it and reinsert all previously found
                destination = pickedID;
                for (uint32_t j = 0; j < i; ++j) {
                  insertActive(pickedIDs[j], task_id);
                }
                return true;
              } else {
                pickedIDs[i] = pickedID;
                if (part_of_edges_seen_before_failure > best_part_of_edges_before_fail) {
                  best_idx = i;
                  best_part_of_edges_before_fail = part_of_edges_seen_before_failure;
                }
              }
            }
            // If tried as often as possible, reinsert all but the one that failed at the latest (heuristic)
            _calls_to_pick_that_reached_max_retries.fetch_add(1, std::memory_order_relaxed);
            destination = pickedIDs[best_idx];
            for (uint32_t j = 0; j < _region_similarity_retries; ++j) {
              if (j != best_idx) {
                insertActive(pickedIDs[j], task_id);
              }
            }
            ASSERT(destination != invalidGroupID);
            return true;
          }
        }

        /// Convenience method to call activate() on all successors
        void activateAllSuccessors(const ContractionGroupID id, const size_t task_id) {
          ASSERT(_hierarchy);
            auto succs = _hierarchy->successors(id);
            for (auto s : succs) {
                activate(s, task_id);
            }
        }

        void activate(const ContractionGroupID id, const size_t task_id) {
            insertActive(id, task_id);
        }

        // ! Mark a given id as accepted, i.e. mark that a thread will commence working on uncoarsening the group and
        // ! that it will not be reinserted into the pool again.
        void markAccepted(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          if (!_use_multiqueue) {
            ASSERT(_dpq_active_ids);
            _dpq_active_ids->increment_finished(_hierarchy->depth(id));
          }
          _num_accepted_uncontractions.fetch_add(group(id).size(), std::memory_order_relaxed);
        }

        bool allAccepted() const {
          return _num_accepted_uncontractions.load(std::memory_order_relaxed) == getTotalNumUncontractions();
//          return _active_ids.allDepthsCompleted();
        }

        size_t getNumAccepted() const {
          return _num_accepted_uncontractions.load(std::memory_order_relaxed);
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

        size_t getTotalCallsToPick() const {
          return _total_calls_to_pick.load(std::memory_order_relaxed);
        }

        size_t getCallsToPickThatReachedMaxRetries() const {
          return _calls_to_pick_that_reached_max_retries.load(std::memory_order_relaxed);
        }

        size_t getCallsToPickWithEmptyPQ() const {
          return _calls_to_pick_with_empty_pq.load(std::memory_order_relaxed);
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* hierarchy_node = parent->addChild("Uncontraction Hierarchy");
          _hierarchy->memoryConsumption(hierarchy_node);
          if (!_use_multiqueue) {
            ASSERT(_dpq_active_ids);
            utils::MemoryTreeNode* depth_pq_node = parent->addChild("Depth Priority Queue");
            _dpq_active_ids->memoryConsumption(depth_pq_node);
          }
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


        bool tryPopActive(ContractionGroupID& destination, const size_t task_id) {
//                _queue_lock.lock();
//                destination = _active_ids.top();
//                _active_ids.pop();
//                _queue_lock.unlock();
//                ASSERT(destination != invalidGroupID);
//                return true;

          if (_use_multiqueue) {
            ASSERT(_mq_active_ids);
            auto handle = _mq_active_ids->get_handle(task_id);
            typename Multiqueue::value_type ret;
            bool success = _mq_active_ids->extract_top(handle, ret);
            if (success) {
              destination = ret.second;
            }
            return success;
          } else {
            ASSERT(_dpq_active_ids);
            unused(task_id);
            return _dpq_active_ids->try_pop(destination);
          }
        }

        void insertActive(ContractionGroupID id, const size_t task_id) {
            ASSERT(_hierarchy);
          if (_use_multiqueue) {
            ASSERT(_mq_active_ids);
            auto handle = _mq_active_ids->get_handle(task_id);
            _mq_active_ids->push(handle, {_hierarchy->depth(id), id});
          } else {
            ASSERT(_dpq_active_ids);
            unused(task_id);
            _dpq_active_ids->push(id, _hierarchy->depth(id));
          }
        }

        struct pick_candidate {
            ContractionGroupID group_id = invalidGroupID;
            double similarity_to_calling_task = 0.0;
            double similarity_to_other_tasks = 1.0;
        };

        bool tryToGetNewCandidate(pick_candidate& candidate, const size_t task_id, const HyperedgeID num_active_in_this_and_other_tasks) {
          bool picked = tryPopActive(candidate.group_id, task_id);
          if (!picked) return false;

          HypernodeID repr = group(candidate.group_id).getRepresentative();
          auto [sim_to_this_task, sim_to_other_tasks] = _node_region_comparator->regionSimilarityToThisTaskAndOtherTasks(repr, task_id, num_active_in_this_and_other_tasks);
          candidate.similarity_to_calling_task = sim_to_this_task;
          candidate.similarity_to_other_tasks = sim_to_other_tasks;
          return true;
        }

        bool isGoodCandidate(const pick_candidate& candidate, const bool accept_based_only_on_similarity_to_others) const {
          return (candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON
            && (accept_based_only_on_similarity_to_others || candidate.similarity_to_calling_task > _similarity_to_calling_task_threshold - CMP_EPSILON));
        }

        // ! Returns true iff this_candidate is better than other_candidate
        bool isBetterThanOtherCandidate(const pick_candidate& this_candidate, const pick_candidate& other_candidate) const {

          // If both are below threshold for other tasks, the one with greater similarity to this task is better
          if (this_candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON
              && other_candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON) {
            return this_candidate.similarity_to_calling_task > other_candidate.similarity_to_calling_task;
          }

          // Otherwise the one with the lower similarity to other tasks is better
          return (this_candidate.similarity_to_other_tasks < other_candidate.similarity_to_other_tasks
            || (this_candidate.similarity_to_other_tasks <= other_candidate.similarity_to_other_tasks + CMP_EPSILON
              && this_candidate.similarity_to_calling_task > other_candidate.similarity_to_calling_task));
        }

        std::unique_ptr<GroupHierarchy> _hierarchy;

        // todo: Idea: PQs by depth are somewhat irrelevant because threads explicitly work on different regions now
        // => shuffle roots, divide roots into p parts, give each task its own part of the roots; each task grows its own pool, later work stealing
        std::unique_ptr<Multiqueue> _mq_active_ids;
        std::unique_ptr<DepthPriorityQueue> _dpq_active_ids;

        RegionComparator* _node_region_comparator;
        const size_t _region_similarity_retries;
        tbb::enumerable_thread_specific<std::vector<ContractionGroupID>> _picked_ids_ets;
        tbb::enumerable_thread_specific<std::vector<pick_candidate>> _pick_candidates_ets;

        CAtomic<size_t> _total_calls_to_pick;
        CAtomic<size_t> _calls_to_pick_that_reached_max_retries;
        CAtomic<size_t> _calls_to_pick_with_empty_pq;

        CAtomic<size_t> _num_accepted_uncontractions;

        const bool _calculate_full_similarities;
        const bool _use_multiqueue;

        const double _similarity_to_others_threshold;
        const double _similarity_to_calling_task_threshold = 0.05;

        static constexpr double CMP_EPSILON = 1.0e-30;


    };

//    using TreeGroupPool = ConcurrentQueueGroupPool<UncontractionGroupTree, Hypergraph>;
//    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar::ds
