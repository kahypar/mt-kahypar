//
// Created by mlaupichler on 24.07.21.
//

#pragma once

#include <tbb/enumerable_thread_specific.h>
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/incident_net_array.h"

namespace mt_kahypar::ds {

    /// Data structure that is supposed to approximate a coefficient that describes the similarity of the HG region of
    /// a node against a set of active nodes, i.e. nodes that are currently being uncontracted or used as seeds for refinement.
    /// The coefficient is in [0,1] where 0 means no similarity and 1 means all incident edges of the new node are
    /// subsumed by the union of incident edges of active nodes.
    template<typename Hypergraph = Mandatory>
    class NodeRegionComparator {

    public:

        explicit NodeRegionComparator(const Hypergraph& hypergraph, const Context& context) :
            _hg(hypergraph),
            _region_similarity_threshold(context.uncoarsening.node_region_similarity_threshold),
            _active_nodes_signature_union_size(0),
            _active_nodes_combined_signatures(_hg.initialNumEdges(), CAtomic<size_t>(0)) {}

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Active Nodes Combined Signatures", sizeof(CAtomic<size_t>) * _active_nodes_combined_signatures.size());
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
        void markActive(const IteratorRange<IncidentEdgeIteratorT> incident_edges) {
          HyperedgeID num_edges_activated = 0;
          for (const HypernodeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            uint32_t count = _active_nodes_combined_signatures[he].add_fetch(1, std::memory_order_relaxed);
            if (count == 1) {
//              _active_nodes_signature_union_size.fetch_add(1, std::memory_order_relaxed);
                ++num_edges_activated;
            }
          }
          HyperedgeID num_active = _active_nodes_signature_union_size.fetch_add(num_edges_activated, std::memory_order_relaxed);
          unused(num_active);
          ASSERT(num_active < num_active + num_edges_activated);
        }


        template<typename IncidentEdgeIteratorT1, typename IncidentEdgeIteratorT2>
        void markInactive(const IteratorRange<IncidentEdgeIteratorT1> current_incident_edges, const IteratorRange<IncidentEdgeIteratorT2> dropped_incident_edges) {
          HyperedgeID num_edges_deactivated = 0;
          markInactive(current_incident_edges, num_edges_deactivated);
          markInactive(dropped_incident_edges, num_edges_deactivated);
          uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(num_edges_deactivated, std::memory_order_relaxed);
          unused(union_size);
          ASSERT(union_size >= num_edges_deactivated);
        }

        bool regionIsNotTooSimilarToActiveNodesWithFullSimilarity(const HypernodeID hn, double& similarity) const {
          similarity = regionSimilarityToActiveNodes(hn);
          return similarity < _region_similarity_threshold + CMP_EPSILON;
        }

        double regionSimilarityToActiveNodes(const HypernodeID hn) const {
          ASSERT(_hg.nodeIsEnabled(hn));
          size_t union_size = _hg.nodeDegree(hn) + _active_nodes_signature_union_size.load(std::memory_order_relaxed);
          size_t intersection_size = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            size_t num_active_incident_to_he = _active_nodes_combined_signatures[he].load(std::memory_order_relaxed);
            if (num_active_incident_to_he > 0) {
              ++intersection_size;
            }
          }
          // |A \cup B| = |A| + |B| - |A \cap B|
          ASSERT(union_size >= intersection_size);
          union_size -= intersection_size;

          if (union_size == 0) return 0.0;
          return (double) intersection_size / (double) union_size;
        }

        bool regionIsNotTooSimilarToActiveNodesWithEarlyBreak(const HypernodeID hn) const {
          ASSERT(_hg.nodeIsEnabled(hn));
          size_t union_size = _hg.nodeDegree(hn) + _active_nodes_signature_union_size.load(std::memory_order_relaxed);
          size_t intersection_size = 0;
          auto incident_edges = _hg.incidentEdges(hn);
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            size_t num_active_incident_to_he = _active_nodes_combined_signatures[he].load(std::memory_order_relaxed);
            if (num_active_incident_to_he > 0) {
              ++intersection_size;
              // |A \cup B| = |A| + |B| - |A \cap B|
              ASSERT(union_size > 0);
              --union_size;
              // Break early if threshold has been surpassed already
              double lower_bound_sim = (double) intersection_size / (double) union_size;
              if (lower_bound_sim > _region_similarity_threshold + CMP_EPSILON) return false;
            }
          }
          ASSERT(union_size >= intersection_size);
          return true;
        }

        // ! Only for testing/debugging
        bool checkNoneActiveParallel() const {
          if (_active_nodes_signature_union_size.load(std::memory_order_relaxed) != 0) return false;
          bool noneActive = true;
          tbb::parallel_for(ID(0), ID(_active_nodes_combined_signatures.size()), [&](const HyperedgeID he) {
            bool active = _active_nodes_combined_signatures[he].load(std::memory_order_relaxed) > 0;
            if (active) {
              noneActive = false;
            }
          });
          return noneActive;
        }

    private:

        template<typename IncidentEdgeIteratorT>
        void markInactive(const IteratorRange<IncidentEdgeIteratorT> incident_edges, HyperedgeID& num_edges_deactivated) {
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            uint32_t count = _active_nodes_combined_signatures[he].sub_fetch(1, std::memory_order_relaxed);
            ASSERT(count < count + 1);
            if (count == 0) {
//              uint32_t union_size = _active_nodes_signature_union_size.fetch_sub(1, std::memory_order_relaxed);
//              unused(union_size);
//              ASSERT(union_size > 0);
              ++num_edges_deactivated;
            }
          }
        }

        const Hypergraph& _hg;
        const double _region_similarity_threshold;

        CAtomic<HyperedgeID> _active_nodes_signature_union_size;
        Array<CAtomic<size_t>> _active_nodes_combined_signatures;

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
