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

        explicit NodeRegionComparator(const Hypergraph& hypergraph) :
            _hg(hypergraph),
            _active_nodes_signature_union_size(0),
            _active_nodes_combined_signatures(_hg.initialNumEdges(), CAtomic<uint32_t>(0)) {}

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Active Nodes Combined Signatures", sizeof(CAtomic<uint32_t>) * _active_nodes_combined_signatures.size());
        }

        template<typename IncidentEdgeIteratorT>
        void markActive(const IteratorRange<IncidentEdgeIteratorT> incident_edges) {
          for (const HypernodeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            uint32_t count = _active_nodes_combined_signatures[he].add_fetch(1, std::memory_order_relaxed);
            if (count == 1) {
              _active_nodes_signature_union_size.fetch_add(1, std::memory_order_relaxed);
            }
          }
        }

        template<typename IncidentEdgeIteratorT>
        void markInactive(const IteratorRange<IncidentEdgeIteratorT> incident_edges) {
          for (const HyperedgeID he : incident_edges) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            uint32_t count = _active_nodes_combined_signatures[he].sub_fetch(1, std::memory_order_relaxed);
            ASSERT(count < count + 1);
            if (count == 0) {
              uint32_t union_size = _active_nodes_signature_union_size.sub_fetch(1, std::memory_order_relaxed);
              unused(union_size);
              ASSERT(union_size < union_size + 1);
            }
          }
        }

        template<typename IncidentEdgeIteratorT1, typename IncidentEdgeIteratorT2>
        void markInactive(const IteratorRange<IncidentEdgeIteratorT1> current_incident_edges, const IteratorRange<IncidentEdgeIteratorT2> dropped_incident_edges) {
          markInactive(current_incident_edges);
          markInactive(dropped_incident_edges);
        }

        double regionSimilarityToActiveNodes(const HypernodeID hn) const {
          HyperedgeID union_size = 0;
          HyperedgeID intersection_size = 0;
          for (const HyperedgeID he : _hg.incidentEdges(hn)) {
            ASSERT(he < _active_nodes_combined_signatures.size());
            ++union_size;
            if (_active_nodes_combined_signatures[he].load(std::memory_order_relaxed) > 0) {
              ++intersection_size;
            }
          }
          ASSERT(union_size >= intersection_size);
          // |A \cup B| = |A| + |B| - |A \cap B|
          HypernodeID active_nodes_union_size = _active_nodes_signature_union_size.load(std::memory_order_relaxed);
          union_size += active_nodes_union_size - intersection_size;

          if (union_size == 0) return 0.0;
          return (double) intersection_size / (double) union_size;
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

        const Hypergraph& _hg;

        CAtomic<HyperedgeID> _active_nodes_signature_union_size;
        Array<CAtomic<uint32_t>> _active_nodes_combined_signatures;


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

    };

} // end namespace mt_kahypar::ds
