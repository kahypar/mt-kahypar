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

    private:

        class ActiveEdges {

        public:

            explicit ActiveEdges(const size_t initial_capacity) :
              _num_expired_threshold(std::max(initial_capacity / 4, size_t(1))),
              _first_valid(0),
              _first_invalid(0),
              _edges(initial_capacity, kInvalidHyperedge) {}

            ActiveEdges(const ActiveEdges& other) = default;
            ActiveEdges& operator=(const ActiveEdges& other) = default;
            ActiveEdges(ActiveEdges&& other)  noexcept = default;
            ActiveEdges& operator=(ActiveEdges&& other)  noexcept = default;

            size_t size() const {
              return _first_invalid - _first_valid;
            }

            bool empty() const {
              return _first_invalid == _first_valid;
            }

            size_t size_in_bytes() const {
              return (size_t) _edges.size() * sizeof(HyperedgeID);
            };

            void push_back(const HyperedgeID he) {
              if (_first_invalid == _edges.size()) {
                if (_first_valid >= _num_expired_threshold) {
                  // If at least as many elements as the threshold are expired, remove expired part
                  remove_expired();
                } else {
                  // Else have std::vector deal with the push_back on full vector and return
                  _edges.push_back(he);
                  _num_expired_threshold = _edges.size() / 4;
                  ++_first_invalid;
                  return;
                }
              }
              ASSERT(_first_invalid < _edges.size());
              _edges[_first_invalid] = he;
              ++_first_invalid;
            }

            HyperedgeID pop_back() {
              ASSERT(_first_invalid > _first_valid);
              const HyperedgeID back = _edges[--_first_invalid];
              if (_first_invalid == _first_valid) {
                // Remove expired in O(1) if valid edges became empty
                remove_expired();
              }
              return back;
            }

            HyperedgeID pop_front() {
              ASSERT(_first_invalid > _first_valid);
              const HyperedgeID front = _edges[_first_valid++];
              if (_first_invalid == _first_valid) {
                // Remove expired in O(1) if valid edges became empty
                remove_expired();
              }
              return front;
            }

            auto begin() {ASSERT(_first_valid <= _edges.size()); return _edges.begin() + _first_valid;};
            auto end() {ASSERT(_first_invalid <= _edges.size()); return _edges.begin() + _first_invalid;};
            auto cbegin() const {ASSERT(_first_valid <= _edges.size()); return _edges.cbegin() + _first_valid;};
            auto cend() const {ASSERT(_first_invalid <= _edges.size()); return _edges.cbegin() + _first_invalid;};
            auto begin() const {ASSERT(_first_valid <= _edges.size()); return _edges.begin() + _first_valid;};
            auto end() const {ASSERT(_first_invalid <= _edges.size()); return _edges.begin() + _first_invalid;};

        private:

            void remove_expired() {
              for (size_t i = 0; i < size(); ++i) {
                _edges[i] = _edges[i + _first_valid];
              }
              _first_invalid -= _first_valid;
              _first_valid = 0;
            }

            size_t _num_expired_threshold;
            // Marks first index of the currently valid elements (everything before this is considered expired)
            size_t _first_valid;
            // Marks first index past the currently valid elements
            size_t _first_invalid;
            vec<HyperedgeID> _edges;
        };

    public:

        explicit NodeRegionComparator(const Hypergraph& hypergraph, const double similarity_threshold, const size_t num_threads) :
            _hg(hypergraph),
            _region_similarity_threshold(similarity_threshold),
            _active_nodes_signature_union_size(0),
            _active_nodes_combined_signatures(_hg.initialNumEdges(), num_threads),
            _active_edges_per_task(num_threads) {
          for (size_t tid = 0; tid < num_threads; ++tid) {
            _active_edges_per_task[tid] = std::make_unique<ActiveEdges>(std::max<size_t>(static_cast<size_t>(_hg.initialNumEdges() / num_threads), size_t(16)));
          }
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Active Nodes Combined Signatures", _active_nodes_combined_signatures.size_in_bytes());
          size_t combined_size_of_edge_buffers = 0;
          for (const auto& active_edges_in_task : _active_edges_per_task) {
            combined_size_of_edge_buffers += active_edges_in_task->size_in_bytes();
          }
          parent->addChild("Active Edges Per Task", combined_size_of_edge_buffers);
        }

        template<typename IncidentEdgeIteratorT>
        void markActive(const IteratorRange<IncidentEdgeIteratorT> incident_edges, const size_t task_id, HyperedgeID& num_edges_activated_for_this_task) {
          ASSERT(task_id < _active_edges_per_task.size());
          HyperedgeID num_edges_activated = 0;
          for (const HypernodeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool activated_edge = false;
            bool changed_bit = _active_nodes_combined_signatures.set_true(he, task_id, activated_edge);
            if (changed_bit) {
              _active_edges_per_task[task_id]->push_back(he);
              ++num_edges_activated_for_this_task;
            }
            if (activated_edge) {
                ++num_edges_activated;
            }
          }
          HyperedgeID num_active = _active_nodes_signature_union_size.fetch_add(num_edges_activated, std::memory_order_relaxed);
          unused(num_active);
          ASSERT(num_active <= num_active + num_edges_activated);
        }

        void markLastActivatedEdgesForTaskInactive(const size_t task_id, const HyperedgeID num_edges_to_deactivate_for_task) {
          ASSERT(task_id < _active_edges_per_task.size());
          ActiveEdges& active_in_task = *_active_edges_per_task[task_id];
          HyperedgeID num_edges_deactivated = 0;
          for (HyperedgeID i = 0; i < num_edges_to_deactivate_for_task; ++i) {
            const HyperedgeID he = active_in_task.pop_back();
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool deactivated_edge = false;
            bool changed_bit = _active_nodes_combined_signatures.set_false(he, task_id, deactivated_edge);
            ASSERT(changed_bit); unused(changed_bit);
            if (deactivated_edge) {
              ++num_edges_deactivated;
            }
          }
          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
          unused(union_size);
          ASSERT(union_size >= num_edges_deactivated);
        }

        void markFirstActivatedEdgesForTaskInactive(const size_t task_id, const HyperedgeID num_edges_to_deactivate_for_task) {
          ASSERT(task_id < _active_edges_per_task.size());
          ActiveEdges& active_in_task = *_active_edges_per_task[task_id];
          HyperedgeID num_edges_deactivated = 0;
          for (HyperedgeID i = 0; i < num_edges_to_deactivate_for_task; ++i) {
            const HyperedgeID he = active_in_task.pop_front();
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            bool deactivated_edge = false;
            bool changed_bit = _active_nodes_combined_signatures.set_false(he, task_id, deactivated_edge);
            ASSERT(changed_bit); unused(changed_bit);
            if (deactivated_edge) {
              ++num_edges_deactivated;
            }
          }
          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
          unused(union_size);
          ASSERT(union_size >= num_edges_deactivated);
        }

        void markAllEdgesForTaskInactive(const size_t task_id) {
          ASSERT(task_id < _active_edges_per_task.size());
          const HyperedgeID total_num_active_in_task = _active_edges_per_task[task_id]->size();
          markLastActivatedEdgesForTaskInactive(task_id, total_num_active_in_task);
          ASSERT(_active_edges_per_task[task_id]->empty());
        }

        bool isInGoodRegionWithFullSimilarity(const HypernodeID hn, const size_t task_id, const size_t num_active_in_this_and_other_tasks,
                                              double &similarity_to_this_task, double& similarity_to_other_tasks) {
          auto [sim_to_this_task, sim_to_other_tasks] = regionSimilarityToThisTaskAndOtherTasks(hn, task_id, num_active_in_this_and_other_tasks);
          similarity_to_this_task = sim_to_this_task;
          similarity_to_other_tasks = sim_to_other_tasks;
          return (similarity_to_other_tasks < _region_similarity_threshold + CMP_EPSILON) && (similarity_to_this_task > (1.0 - _region_similarity_threshold) - CMP_EPSILON);
        }

        bool regionIsNotTooSimilarToActiveNodesWithEarlyBreak(const HypernodeID hn, const size_t task_id, double& part_of_edges_seen_before_failure) {
          ASSERT(_hg.nodeIsEnabled(hn));
          if (_region_similarity_threshold <= 0.0 + CMP_EPSILON) {
            return regionDoesNotIntersectsWithActiveNodes(hn, task_id, part_of_edges_seen_before_failure);
          }
          size_t union_size_lower_bound = _hg.nodeDegree(hn) + _active_nodes_signature_union_size.load(std::memory_order_relaxed) - _active_edges_per_task[task_id]->size();
          size_t intersection_size_with_other_tasks = 0;
          HyperedgeID num_edges_seen = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            auto set_state = _active_nodes_combined_signatures.set_state_for_task(he, task_id);
            if (set_state.is_set_for_any_except_task) {
              ++intersection_size_with_other_tasks;
            }
            if (set_state.is_set_for_any_task()) {
              // |A \cup B| = |A| + |B| - |A \cap B|
              ASSERT(union_size_lower_bound > 0);
              --union_size_lower_bound;
            }

            if (set_state.is_set_for_any_task()) {
              // Break early if threshold has been surpassed already
              double lower_bound_sim = (double) intersection_size_with_other_tasks / (double) union_size_lower_bound;
              if (lower_bound_sim > _region_similarity_threshold - CMP_EPSILON) {
                part_of_edges_seen_before_failure = (double) num_edges_seen / (double) _hg.nodeDegree(hn);
                return false;
              }
            }

            ++num_edges_seen;
          }
          ASSERT(union_size_lower_bound >= intersection_size_with_other_tasks);
          ASSERT(num_edges_seen == _hg.nodeDegree(hn));
          part_of_edges_seen_before_failure = 1.0;
          return true;
        }

        // Simplified version of regionIsNotTooSimilarToActiveNodesWithEarlyBreak in case the similarity threshold is 0.0.
        // Simply returns true exactly if the given node shares no incident edges to the active nodes.
        bool regionDoesNotIntersectsWithActiveNodes(const HypernodeID hn, const size_t task_id, double& part_of_edges_seen_before_failure) {
          ASSERT(_region_similarity_threshold <= 0.0 + CMP_EPSILON);
          ASSERT(_hg.nodeIsEnabled(hn));
          HyperedgeID num_edges_seen = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const auto he : incident_edges) {
            if (_active_nodes_combined_signatures.any_set_except_thread(he, task_id)) {
              part_of_edges_seen_before_failure = (double) num_edges_seen / (double) _hg.nodeDegree(hn);
              return false;
            }
            ++num_edges_seen;
          }

          ASSERT(num_edges_seen == _hg.nodeDegree(hn));
          part_of_edges_seen_before_failure = 1.0;
          return true;
        }

        HyperedgeID numberOfEdgesActiveInThisTaskAndAnyOtherTask(const size_t task_id) {
          ASSERT(task_id < _active_edges_per_task.size());
          ActiveEdges& active_in_task = *_active_edges_per_task[task_id];
          size_t intersection_size = 0;
          for (const auto& he : active_in_task) {
            if (_active_nodes_combined_signatures.any_set_except_thread(he, task_id)) {
              ++intersection_size;
            }
          }
          return intersection_size;
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

        std::pair<double, double> regionSimilarityToThisTaskAndOtherTasks(const HypernodeID hn, const size_t task_id, const size_t num_active_in_this_and_other_tasks) {
          ASSERT(_hg.nodeIsEnabled(hn));
          size_t intersection_size_with_this_task = 0;
          size_t intersection_size_with_other_tasks = 0;
          size_t intersection_size_with_all_tasks = 0;
          size_t intersection_size_with_this_and_any_other_tasks = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.numElements());
            auto set_state = _active_nodes_combined_signatures.set_state_for_task(he, task_id);
            if (set_state.is_set_for_any_task()) {
              ++intersection_size_with_all_tasks;
            }
            if (set_state.is_set_for_task) {
              ++intersection_size_with_this_task;
            }
            if (set_state.is_set_for_any_except_task) {
              ++intersection_size_with_other_tasks;
            }
            if (set_state.is_set_for_task && set_state.is_set_for_any_except_task) {
              ++intersection_size_with_this_and_any_other_tasks;
            }
          }
          // calculate size of union between active edges in other threads and edges incident to given node without
          // explicitly knowing the former set
          const size_t union_size_with_others = _active_nodes_signature_union_size.load(std::memory_order_relaxed) // total number active edges
                  + _hg.nodeDegree(hn) // total number incident edges
                  - _active_edges_per_task[task_id]->size() // number active in this task
                  - intersection_size_with_all_tasks;
                  + intersection_size_with_this_task
                  + num_active_in_this_and_other_tasks
                  - intersection_size_with_this_and_any_other_tasks;
          ASSERT(union_size_with_others >= intersection_size_with_other_tasks);

          double similarity_to_others = 0.0;
          if (union_size_with_others > 0) {
            similarity_to_others = (double) intersection_size_with_other_tasks / (double) union_size_with_others;
          }

          const size_t union_size_with_this_task = _active_edges_per_task[task_id]->size() + _hg.nodeDegree(hn) - intersection_size_with_this_task;

          double similarity_to_this = 0.0;
          if (union_size_with_this_task > 0) {
            similarity_to_this = (double) intersection_size_with_this_task / (double) union_size_with_this_task;
          }

          return std::make_pair(similarity_to_this, similarity_to_others);
        }

    private:

        const Hypergraph& _hg;
        const double _region_similarity_threshold;

        CAtomic<HyperedgeID> _active_nodes_signature_union_size;
        ThreadWiseFlagArray<HyperedgeID> _active_nodes_combined_signatures;

        std::vector<std::unique_ptr<ActiveEdges>> _active_edges_per_task;

        static constexpr double CMP_EPSILON = 1.0e-100;

    };

    template<class Hypergraph = Mandatory>
    class DummyRegionComparator {

    public:
        explicit DummyRegionComparator(const Hypergraph&) {};

        void memoryConsumption(utils::MemoryTreeNode*) const {}

        template<typename IncidentEdgeIteratorT>
        void markActive(const IteratorRange<IncidentEdgeIteratorT>, const size_t, HyperedgeID&) {}

        void markLastActivatedEdgesForTaskInactive(const size_t, const HyperedgeID) {}

        void markAllEdgesForTaskInactive(const size_t) {}

        bool regionIsNotTooSimilarToActiveNodesWithFullSimilarity(const HypernodeID, double& similarity) const {
          similarity = 0.0;
          return true;
        }

        bool regionIsNotTooSimilarToActiveNodesWithEarlyBreak(const HypernodeID) const {
          return true;
        }

    };

} // end namespace mt_kahypar::ds
