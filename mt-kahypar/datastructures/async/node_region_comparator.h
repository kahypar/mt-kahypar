//
// Created by mlaupichler on 24.07.21.
//

#pragma once

#include <tbb/enumerable_thread_specific.h>
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/incident_net_array.h"

namespace mt_kahypar::ds {

    /// Data structure that is supposed to approximate a coefficient that describes the similarity of the HG region of
    /// two nodes. The coefficient is in [0,1] where 0 means no similarity and 1 means equivalent regions. The coefficient
    /// is based on a signature set of a constant number of the largest incident edges of a hypernode and equals the Jaccard-Index
    /// between these signature sets of two nodes. The approximation only guarantees a similarity value of 0 for two nodes that are not
    /// neighbors, i.e. share no incident edges, including degree-zero nodes and a value of 1 for nodes that have the
    /// exact same and more than zero incident edges. Larger signature sets promise more accurate region similarity
    /// values but require more memory, linear in the number of hypernodes.
    template<typename Hypergraph = Mandatory>
    class NodeRegionComparator {

    private:

        using SizeCmpFunc = std::function<bool (const HyperedgeID, const HyperedgeID)>;

    public:

        explicit NodeRegionComparator(const size_t signature_size = DEFAULT_SIGNATURE_SIZE) :
            _hg(nullptr),
            _signatures(),
            _signature_size(signature_size),
            _sort_buffer_ets() {}

        // ! Calculates or recalculates the signature sets for every node of a given hypergraph in parallel.
        void calculateSignaturesParallel(const Hypergraph* hypergraph) {
          ASSERT(hypergraph);
          if (_hg != hypergraph) {
            _hg = hypergraph;
            _signatures = Array<HyperedgeID>(_hg->initialNumNodes() * _signature_size, kInvalidHyperedge);
          }
          ASSERT(_hg);
          const HyperedgeID max_degree = _hg->maxNodeDegree();
          _sort_buffer_ets = tbb::enumerable_thread_specific<std::vector<HyperedgeID>>(max_degree, kInvalidHyperedge);

          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
          };

          tbb::parallel_for(ID(0), _hg->initialNumNodes(), [&](const HypernodeID hn) {
            ASSERT(hn * _signature_size < _signatures.size());
            auto& local_sort_buffer = _sort_buffer_ets.local();
            ASSERT(_hg->nodeDegree(hn) <= local_sort_buffer.size(), V(_hg->nodeDegree(hn)) << V(local_sort_buffer.size()));
            calculateSignatureOfNode(hn, local_sort_buffer, size_cmp);
          });
        }

        // ! Returns the approximate region similarity coefficient for two hypernodes. Runs in O(|signature|).
        double regionSimilarity(const HypernodeID hn1, const HypernodeID hn2) const {
          ASSERT(_hg);
          ASSERT(hn1 * _signature_size < _signatures.size());
          ASSERT(hn2 * _signature_size < _signatures.size());
          if (hn1 == hn2) return 1;
          HyperedgeID intersection_size = 0;
          HyperedgeID union_size = 0;
          HyperedgeID idx1 = 0;
          HyperedgeID idx2 = 0;
          while (idx1 < _signature_size && idx2 < _signature_size && _signatures[sig_entry(hn1, idx1)] != kInvalidHyperedge && _signatures[sig_entry(hn2, idx2)] != kInvalidHyperedge) {
            if (_signatures[sig_entry(hn1, idx1)] == _signatures[sig_entry(hn2, idx2)]) {
              ++intersection_size;
              ++union_size;
              ++idx1;
              ++idx2;
            } else if (_hg->edgeSize(_signatures[sig_entry(hn1, idx1)]) < _hg->edgeSize(_signatures[sig_entry(hn2, idx2)])) {
              ++union_size;
              ++idx1;
            } else {
//              ASSERT(_hg->edgeSize(_signatures[sig_entry(hn1, idx1)]) >= _hg->edgeSize(_signatures[sig_entry(hn2, idx2)]));
              ++union_size;
              ++idx2;
            }
          }

          while (idx1 < _signature_size && _signatures[sig_entry(hn1, idx1)] != kInvalidHyperedge) {
            ++union_size;
            ++idx1;
          }

          while (idx2 < _signature_size && _signatures[sig_entry(hn2, idx2)] != kInvalidHyperedge) {
            ++union_size;
            ++idx2;
          }

          if (union_size == 0) return 0;

          return (double) intersection_size / (double) union_size;
        }

        void adaptForPinRemoval(const HypernodeID removed_pin, const HyperedgeID he) {
          ASSERT(_hg);
          HypernodeID smallestSizeInSignature = _hg->edgeSize(_signatures[sig_entry(removed_pin, 0)]);
          HypernodeID sizeOfHEBeforeRemoval = _hg->edgeSize(he) + 1;
          if (sizeOfHEBeforeRemoval >= smallestSizeInSignature) {
            recalculateSignatureOfNode(removed_pin);
          }
        }

        void adaptForNewPin(const HypernodeID new_pin, const HyperedgeID he) {
          ASSERT(_hg);
          recalculateSignatureOfNode(new_pin);

          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
          };

          for (auto pin : _hg->pins(he)) {
            if (pin == new_pin) continue;
            ASSERT(std::is_sorted(_signatures.begin() + sig_entry(pin, 0), _signatures.begin() + _hg->nodeDegree(pin), size_cmp));
            // If edge is still smaller than smallest edge in signature, no change
            if (size_cmp(he, _signatures[sig_entry(pin, 0)])) continue;
            // Else find slot where edge belongs now in signature of pin
            HyperedgeID i = 1;
            while (!size_cmp(he, _signatures[sig_entry(pin, i)])) {
              ++i;
            }
            ASSERT(i >= 1);
            const HyperedgeID new_slot_of_he = i - 1;
            // Move everything left of new slot of he one slot further left until the old spot of he or until index 0 of
            // the signature in case he was previously not in the signature
            HyperedgeID old = _signatures[sig_entry(pin, new_slot_of_he)];
            HyperedgeID j = new_slot_of_he;
            while (j > 0 && old != he) {
              HyperedgeID overwritten_value = _signatures[sig_entry(pin, j-1)];
              _signatures[sig_entry(pin, j-1)] = old;
              old = overwritten_value;
              --j;
            }
            // Write he to its new slot
            _signatures[sig_entry(pin, new_slot_of_he)] = he;
            ASSERT(std::is_sorted(_signatures.begin() + sig_entry(pin, 0), _signatures.begin() + _hg->nodeDegree(pin), size_cmp));
          }
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Signature Sets", sizeof(HyperedgeID) * _signatures.size());
        }

    private:

        void calculateSignatureOfNode(const HypernodeID hn, std::vector<HypernodeID>& sort_buffer, SizeCmpFunc& size_cmp) {
          ASSERT(_hg);
          ASSERT(hn * _signature_size < _signatures.size());

          auto incident_nets_range = _hg->incidentEdges(hn);
          // Note: Using _hg->nodeDegree() here does not work as that returns the degree in the original graph, not the
          // current one
          const HyperedgeID degree = std::distance(incident_nets_range.begin(), incident_nets_range.end());

          if (degree > sort_buffer.size()) {
            sort_buffer.resize(degree, kInvalidHyperedge);
          }

          if (degree <= _signature_size) {
            // Simple case: Incident nets fit entirely into signature: Copy the edge ids in, sort them and return
            HyperedgeID i = 0;
            for (auto it = incident_nets_range.begin(); it != incident_nets_range.end(); ++it) {
              ASSERT( i < _signature_size);
              _signatures[sig_entry(hn, i)] = *it;
              ++i;
            }
            std::sort(_signatures.begin() + sig_entry(hn, 0), _signatures.begin() + sig_entry(hn , 0) + degree, size_cmp);
            for (HyperedgeID j = degree; j < _signature_size; ++j) {
              _signatures[sig_entry(hn, j)] = kInvalidHyperedge;
            }
          } else {
            // General case: If degree > _signature_size find _signature_size largest edges incident to nodes and use those as signature
            ASSERT(sort_buffer.size() >= degree);
            HypernodeID idx = 0;
            for (IncidentNetIterator it = incident_nets_range.begin(); it != incident_nets_range.end(); ++it) {
              ASSERT(idx < degree);
              sort_buffer[idx] = *it;
              ++idx;
            }
            std::sort(sort_buffer.begin(), sort_buffer.begin() + degree, size_cmp);
            HyperedgeID signature_begin_index = degree - _signature_size;
            ASSERT(signature_begin_index < degree);
            for (HyperedgeID i = 0; i < _signature_size; ++i) {
              _signatures[sig_entry(hn, i)] = sort_buffer[signature_begin_index + i];
            }
          }
        }

        void recalculateSignatureOfNode(const HypernodeID hn) {
          ASSERT(_hg);
          ASSERT(hn * _signature_size < _signatures.size());

          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
          };

          auto& sort_buffer = _sort_buffer_ets.local();
          calculateSignatureOfNode(hn, sort_buffer, size_cmp);
        }

        size_t sig_entry(const HypernodeID hn, const HypernodeID sig_idx) const {
          ASSERT(sig_idx >= 0 && sig_idx < _signature_size);
          return (size_t) hn * _signature_size + sig_idx;
        }

        const Hypergraph* _hg;

        Array<HyperedgeID> _signatures;
        const size_t _signature_size;
        tbb::enumerable_thread_specific<std::vector<HyperedgeID>> _sort_buffer_ets;

        static constexpr size_t DEFAULT_SIGNATURE_SIZE = 50;

    };

} // end namespace mt_kahypar::ds
