//
// Created by mlaupichler on 20.06.21.
//

#pragma once

#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/datastructures/async/array_lock_manager.h>

namespace mt_kahypar {


    class AsyncNodeTracker {

    public:

        explicit AsyncNodeTracker(const HypernodeID node_freeze_time) :
            _node_freeze_time(node_freeze_time),
            _underlying_locks(),
            _time(0),
            _last_tracked_timestamps() {}

        void resize(const HypernodeID num_nodes) {
          _underlying_locks = std::make_unique<ds::GroupLockManager>(num_nodes, ds::invalidGroupID);
          _time = CAtomic<HypernodeID>(_node_freeze_time + 1);
          _last_tracked_timestamps.resize(num_nodes, CAtomic<HypernodeID>(0));
        }

        bool tryAcquireNode(const HypernodeID u, const ds::ContractionGroupID search) {
          ASSERT(_underlying_locks);
          bool acquired = !isFrozen(u) && _underlying_locks->tryToAcquireLock(u, search);
          return acquired;
        }

        void releaseNode(const HypernodeID u, const ds::ContractionGroupID search) {
          ASSERT(_underlying_locks);
          _underlying_locks->strongReleaseLock(u,search);
        }

        // Signature without ContractionGroupID to conform to mt_kahypar::NodeTracker method used in FM Strategies
        void releaseNode(const HypernodeID u) {
          ASSERT(_underlying_locks);
          releaseNode(u, _underlying_locks->owner(u));
        }

        ds::ContractionGroupID owner(const HypernodeID u) const {
          return _underlying_locks->owner(u);
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* node_tracker_node = parent->addChild("Async Node Tracker");
          _underlying_locks->memoryConsumption(node_tracker_node);
          node_tracker_node->addChild("Freeze Timestamps", sizeof(CAtomic<HypernodeID>) * _last_tracked_timestamps.size());
        }

        bool claimedFirstTime(const HypernodeID, const SearchID) {
          ERROR("AsyncNodeTracker does not keep track of which nodes have been claimed before.");
          return false;
        }

        // ! Only for testing
        bool checkNoneLocked() {
          return _underlying_locks->checkNoneLocked();
        }

        void incrementTime(const HypernodeID delta) {
          ASSERT(delta != 0);
          HypernodeID new_time = _time.add_fetch(delta, std::memory_order_relaxed);
          HypernodeID old_time = new_time - delta;
          if (new_time < old_time) {
            ERROR("Counting uncontractions in _time in AsyncNodeTracker experienced overflow!");
          }
        }

        void setNodeFreezeTime(const HypernodeID freeze_time) {
          _node_freeze_time = freeze_time;
        }

        void freeze(const HypernodeID hn) {
          ASSERT(hn < _last_tracked_timestamps.size());
          HypernodeID current_time = _time.load(std::memory_order_relaxed);
          _last_tracked_timestamps[hn].store(current_time, std::memory_order_relaxed);
        }

    private:

        bool isFrozen(const HypernodeID hn) const {
          ASSERT(hn < _last_tracked_timestamps.size());
          return _time.load(std::memory_order_relaxed) - _last_tracked_timestamps[hn].load(std::memory_order_relaxed) <= _node_freeze_time;
        }

        HypernodeID _node_freeze_time;

        std::unique_ptr<ds::GroupLockManager> _underlying_locks;

        CAtomic<HypernodeID> _time;
        ds::Array<CAtomic<HypernodeID>> _last_tracked_timestamps;

    };

    struct AsyncFMSharedData {

        using ConcreteSearchID = ds::ContractionGroupID;

        // ! PQ handles shared by all threads (each vertex is only held by one thread)
        vec<PosT> vertexPQHandles;

        // ! num parts
        PartitionID numParts;

        // ! Tracks the current search of a node, and if a node can still be added to an active search
        AsyncNodeTracker nodeTracker;

        // ! Stores the designated target part of a vertex, i.e. the part with the highest gain to which moving is feasible
        vec<PartitionID> targetPart;

        // ! Stat counters for testing/debugging
        CAtomic<size_t> total_moves;
        CAtomic<size_t> total_reverts;
        CAtomic<size_t> total_find_moves_calls;
        CAtomic<size_t> find_moves_calls_with_good_prefix;
        CAtomic<size_t> find_move_retries;
        CAtomic<size_t> total_pushes_pos_gain;
        CAtomic<size_t> total_pushes_non_pos_gain;
        CAtomic<size_t> num_pins_touched_by_delta_gain_cache_updates;
        CAtomic<size_t> num_case_from_zero_gc_updates;
        CAtomic<size_t> num_case_from_one_gc_updates;
        CAtomic<size_t> num_case_to_one_gc_updates;
        CAtomic<size_t> num_case_to_two_gc_updates;

        const bool release_nodes = true;

        explicit AsyncFMSharedData(size_t numNodes = 0, PartitionID numParts = 0, size_t numPQHandles = 0, HypernodeID node_freeze_time = 0) :
          vertexPQHandles(),
          numParts(numParts),
          nodeTracker(node_freeze_time),
          targetPart(),
          total_moves(0),
          total_reverts(0),
          total_find_moves_calls(0),
          find_moves_calls_with_good_prefix(0),
          find_move_retries(0),
          total_pushes_pos_gain(0),
          total_pushes_non_pos_gain(0),
          num_pins_touched_by_delta_gain_cache_updates(0),
          num_case_from_zero_gc_updates(0),
          num_case_from_one_gc_updates(0),
          num_case_to_one_gc_updates(0),
          num_case_to_two_gc_updates(0) {

          tbb::parallel_invoke( [&] {
              nodeTracker.resize(numNodes);
          }, [&] {
              vertexPQHandles.resize(numPQHandles, invalid_position);
          }, [&] {
              targetPart.resize(numNodes, kInvalidPartition);
          });
        }

        AsyncFMSharedData(size_t numNodes, const Context& context) :
            AsyncFMSharedData(
            numNodes,
            context.partition.k,
            getNumberOfPQHandles(context, numNodes),
            context.refinement.fm.async_node_freeze_time)  {}


        static size_t getNumberOfPQHandles(const Context& context, size_t numNodes) {
          if (context.refinement.fm.algorithm == FMAlgorithm::fm_gain_delta) {
            return numNodes * context.partition.k;
          } else {
            return numNodes;
          }
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* async_shared_fm_data_node = parent->addChild("Async Shared FM Data");

          utils::MemoryTreeNode* pq_handles_node = async_shared_fm_data_node->addChild("PQ Handles");
          pq_handles_node->updateSize(vertexPQHandles.capacity() * sizeof(PosT));

          nodeTracker.memoryConsumption(async_shared_fm_data_node);
        }

        double getFractionOfRevertedMoves() const {
          return (double) total_reverts.load(std::memory_order_relaxed) / (double) total_moves.load(std::memory_order_relaxed);
        }

        double getFractionOfFMCallsWithGoodPrefix() const {
          return (double) find_moves_calls_with_good_prefix / (double) total_find_moves_calls;
        }

        double getFractionOfPosGainPushes() const {
          return (double) total_pushes_pos_gain / (double) (total_pushes_pos_gain + total_pushes_non_pos_gain);
        }

    };

} // namespace mt_kahypar
