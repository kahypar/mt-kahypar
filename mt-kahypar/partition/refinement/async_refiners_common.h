//
// Created by mlaupichler on 20.06.21.
//

#pragma once

#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/datastructures/async/array_lock_manager.h>

namespace mt_kahypar {


    class AsyncNodeTracker {

    public:

        AsyncNodeTracker() = default;

        void resize(HypernodeID num_nodes) {
          _underlying_locks = std::make_unique<ds::GroupLockManager>(num_nodes, ds::invalidGroupID);
        }

        bool tryAcquireNode(HypernodeID u, ds::ContractionGroupID search) {
          ASSERT(_underlying_locks);
          return _underlying_locks->tryToAcquireLock(u, search);
        }

        void releaseNode(HypernodeID u, ds::ContractionGroupID search) {
          ASSERT(_underlying_locks);
          _underlying_locks->strongReleaseLock(u,search);
        }

        // Signature without ContractionGroupID to conform to mt_kahypar::NodeTracker method used in FM Strategies
        void releaseNode(HypernodeID u) {
          ASSERT(_underlying_locks);
          _underlying_locks->strongReleaseLock(u,_underlying_locks->owner(u));
        }

        ds::ContractionGroupID owner(HypernodeID u) const {
          return _underlying_locks->owner(u);
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* node_tracker_node = parent->addChild("Async Node Tracker");
          _underlying_locks->memoryConsumption(node_tracker_node);
        }

        bool claimedFirstTime(HypernodeID, SearchID) {
          ERROR("AsyncNodeTracker does not keep track of which nodes have been claimed before.");
          return false;
        }

        // ! Only for testing
        bool checkNoneLocked() {
          return _underlying_locks->checkNoneLocked();
        }

    private:

        std::unique_ptr<ds::GroupLockManager> _underlying_locks;

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

        const bool release_nodes = true;

        explicit AsyncFMSharedData(size_t numNodes = 0, PartitionID numParts = 0, size_t numPQHandles = 0) :
          vertexPQHandles(), //numPQHandles, invalid_position),
          numParts(numParts),
          targetPart() {

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
            getNumberOfPQHandles(context, numNodes))  { }


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

    };

} // namespace mt_kahypar
