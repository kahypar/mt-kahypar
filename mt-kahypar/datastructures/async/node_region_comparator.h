//
// Created by mlaupichler on 24.07.21.
//

#pragma once

#include <tbb/enumerable_thread_specific.h>
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/incident_net_array.h"
#include "mt-kahypar/datastructures/async/thread_wise_flag_array.h"

namespace mt_kahypar::ds {

    /// Data structure that is supposed to approximate a coefficient that describes the similarity of the HG region of
    /// a node against a set of active nodes, i.e. nodes that are currently being uncontracted or used as seeds for refinement.
    /// The coefficient is in [0,1] where 0 means no similarity and 1 means all incident edges of the new node are
    /// subsumed by the union of incident edges of active nodes.
    template<typename Hypergraph = Mandatory>
    class NodeRegionComparator {

    public:

        explicit NodeRegionComparator(const Hypergraph& hypergraph, const double similarity_threshold, const size_t num_threads) :
            _hg(hypergraph),
            _region_similarity_threshold(similarity_threshold),
            _active_nodes_signature_union_size(0),
            _active_nodes_combined_signatures(_hg.initialNumEdges(), num_threads),
            _active_edges_per_task(num_threads) {
          for (size_t tid = 0; tid < num_threads; ++tid) {
            _active_edges_per_task[tid] = std::make_unique<std::vector<HyperedgeID>>();
          }
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Active Nodes Combined Signatures", _active_nodes_combined_signatures.size_in_bytes());
        }

        /*void resetHyperedgeSizesParallel() {
          tbb::parallel_for(ID(0), ID(_hg.initialNumEdges()), [&](const HyperedgeID he) {
              if (!_hg.edgeIsEnabled(he)) {
                _hyperedge_sizes_at_beginning_of_version[he] = 0;
                return;
              }
              HypernodeID edge_size = _hg.edgeSize(he);
              if (edge_size == 0) {
                _hyperedge_sizes_at_beginning_of_version[he] = 0;
              } else {
//                double log_size = std::ceil(std::log(_hg.edgeSize(he)));
                _hyperedge_sizes_at_beginning_of_version[he] = _hg.edgeSize(he);
              }
          });
        }*/

        template<typename IncidentEdgeIteratorT>
        void markActive(const IteratorRange<IncidentEdgeIteratorT> incident_edges, const size_t task_id) {
          ASSERT(task_id < _active_edges_per_task.size());
          HyperedgeID num_edges_activated = 0;
          for (const HypernodeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool activated_edge = false;
            bool changed_bit = _active_nodes_combined_signatures.set_true(he, task_id, activated_edge);
            if (changed_bit) {
              _active_edges_per_task[task_id]->push_back(he);
            }
            if (activated_edge) {
//              _active_nodes_signature_union_size.fetch_add(1, std::memory_order_relaxed);
                ++num_edges_activated;
            }
          }
          HyperedgeID num_active = _active_nodes_signature_union_size.fetch_add(num_edges_activated, std::memory_order_relaxed);
          unused(num_active);
          ASSERT(num_active <= num_active + num_edges_activated);
        }

        void markAllEdgesForTaskInactive(const size_t task_id) {
          ASSERT(task_id < _active_edges_per_task.size());
          std::vector<HyperedgeID>& active_in_task = *_active_edges_per_task[task_id];
          HyperedgeID num_edges_deactivated = 0;
          while (!active_in_task.empty()) {
            const HyperedgeID he = active_in_task.back();
            active_in_task.pop_back();
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool deactivated_edge = false;
            bool changed_bit = _active_nodes_combined_signatures.set_false(he, task_id, deactivated_edge);
            ASSERT(changed_bit); unused(changed_bit);
            if (deactivated_edge) {
              ++num_edges_deactivated;
            }
          }
          ASSERT(_active_edges_per_task[task_id]->empty());
          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
          unused(union_size);
          ASSERT(union_size >= num_edges_deactivated);
        }

//        template<typename IncidentEdgeIteratorT1>
//        void markInactive(const IteratorRange<IncidentEdgeIteratorT1> incident_edges, const size_t task_id) {
//          HyperedgeID num_edges_deactivated = 0;
//          markInactiveImpl(incident_edges, task_id, num_edges_deactivated);
//          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
//          unused(union_size);
//          ASSERT(union_size >= num_edges_deactivated);
//        }
//
//        template<typename IncidentEdgeIteratorT1, typename IncidentEdgeIteratorT2>
//        void markInactive(const IteratorRange<IncidentEdgeIteratorT1> current_incident_edges, const IteratorRange<IncidentEdgeIteratorT2> dropped_incident_edges, const size_t task_id) {
//          HyperedgeID num_edges_deactivated = 0;
//          markInactiveImpl(current_incident_edges, task_id, num_edges_deactivated);
//          markInactiveImpl(dropped_incident_edges, task_id, num_edges_deactivated);
//          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
//          unused(union_size);
//          ASSERT(union_size >= num_edges_deactivated);
//        }

        bool regionIsNotTooSimilarToActiveNodesWithFullSimilarity(const HypernodeID hn, const size_t task_id,
                                                                  double &similarity) {
          similarity = regionSimilarityToActiveNodes(hn, task_id);
          return similarity < _region_similarity_threshold + CMP_EPSILON;
        }

        double regionSimilarityToActiveNodes(const HypernodeID hn, const size_t task_id) {
          ASSERT(_hg.nodeIsEnabled(hn));
          size_t union_size_lower_bound = _hg.nodeDegree(hn) + _active_nodes_signature_union_size.load(std::memory_order_relaxed) - _active_edges_per_task[task_id]->size();
          size_t intersection_size_with_other_tasks = 0;
          size_t intersection_size_with_all_tasks = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            auto [intersects_with_any_task, intersects_with_other_tasks] = _active_nodes_combined_signatures.any_set_with_and_without_thread(he, task_id);
            if (intersects_with_any_task) {
              ++intersection_size_with_all_tasks;
            }
            if (intersects_with_other_tasks) {
              ++intersection_size_with_other_tasks;
            }
          }
          // |A \cup B| = |A| + |B| - |A \cap B|
          ASSERT(union_size_lower_bound >= intersection_size_with_all_tasks);
          union_size_lower_bound -= intersection_size_with_all_tasks;

          if (union_size_lower_bound == 0) return 0.0;
          return (double) intersection_size_with_other_tasks / (double) union_size_lower_bound;
        }

        bool regionIsNotTooSimilarToActiveNodesWithEarlyBreak(const HypernodeID hn, const size_t task_id) {
          ASSERT(_hg.nodeIsEnabled(hn));
          if (_region_similarity_threshold <= 0.0 + CMP_EPSILON) {
            return regionDoesNotIntersectsWithActiveNodes(hn, task_id);
          }
          size_t union_size_lower_bound = _hg.nodeDegree(hn) + _active_nodes_signature_union_size.load(std::memory_order_relaxed) - _active_edges_per_task[task_id]->size();
          size_t intersection_size_with_other_tasks = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            auto [intersects_with_any_task, intersects_with_other_tasks] = _active_nodes_combined_signatures.any_set_with_and_without_thread(he, task_id);
            if (intersects_with_other_tasks) {
              ++intersection_size_with_other_tasks;
            }
            if (intersects_with_any_task) {
              // |A \cup B| = |A| + |B| - |A \cap B|
              ASSERT(union_size_lower_bound > 0);
              --union_size_lower_bound;
            }

            if (intersects_with_any_task || intersects_with_other_tasks) {
              // Break early if threshold has been surpassed already
              double lower_bound_sim = (double) intersection_size_with_other_tasks / (double) union_size_lower_bound;
              if (lower_bound_sim > _region_similarity_threshold + CMP_EPSILON) return false;
            }
          }
          ASSERT(union_size_lower_bound >= intersection_size_with_other_tasks);
          return true;
        }

        // Simplified version of regionIsNotTooSimilarToActiveNodesWithEarlyBreak in case the similarity threshold is 0.0.
        // Simply returns true exactly if the given node shares no incident edges to the active nodes.
        bool regionDoesNotIntersectsWithActiveNodes(const HypernodeID hn, const size_t task_id) {
          ASSERT(_region_similarity_threshold <= 0.0 + CMP_EPSILON);
          ASSERT(_hg.nodeIsEnabled(hn));
          auto incident_edges = _hg.incidentEdges(hn);
          return std::none_of(incident_edges.begin(), incident_edges.end(), [&](const HyperedgeID he) {
            return _active_nodes_combined_signatures.any_set_except_thread(he, task_id);
          });
        }

        // ! Only for testing/debugging
        bool checkNoneActiveParallel() {
          if (_active_nodes_signature_union_size.load(std::memory_order_relaxed) != 0) return false;
          bool noneActive = true;
          tbb::parallel_for(ID(0), ID(_active_nodes_combined_signatures.numElements()), [&](const HyperedgeID he) {
            bool active = _active_nodes_combined_signatures.any_set(he);
            if (active) {
              noneActive = false;
            }
          });
          return noneActive;
        }

    private:

        template<typename IncidentEdgeIteratorT>
        void markInactiveImpl(const IteratorRange<IncidentEdgeIteratorT> incident_edges, const size_t task_id, HyperedgeID& num_edges_deactivated) {
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool deactivated_edge = false;
            _active_nodes_combined_signatures.set_false(he, task_id, deactivated_edge);
            if (deactivated_edge) {
              ++num_edges_deactivated;
            }
          }
        }

        const Hypergraph& _hg;
        const double _region_similarity_threshold;

        CAtomic<HyperedgeID> _active_nodes_signature_union_size;
        ThreadWiseFlagArray<HyperedgeID> _active_nodes_combined_signatures;

        std::vector<std::unique_ptr<std::vector<HyperedgeID>>> _active_edges_per_task;

        static constexpr double CMP_EPSILON = 1.0e-100;

    };

    template<class Hypergraph = Mandatory>
    class DummyRegionComparator {

    public:
        explicit DummyRegionComparator(const Hypergraph&) {};

        void memoryConsumption(utils::MemoryTreeNode*) const {}

        void markActive(const HypernodeID) {}

        void markInactive(const HypernodeID) {}

        double regionSimilarityToActiveNodes(const HypernodeID) const {
          return 0.0;
        }

        bool regionIsNotTooSimilarToActiveNodesWithFullSimilarity(const HypernodeID, double& similarity) const {
          similarity = 0.0;
          return true;
        }

        bool regionIsNotTooSimilarToActiveNodesWithEarlyBreak(const HypernodeID) const {
          return true;
        }

    };

} // end namespace mt_kahypar::ds
