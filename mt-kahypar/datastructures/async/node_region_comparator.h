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
    /// neighbors, i.e. share no incident edges, including degree-zero nodes, and a value of 1 for nodes that have the
    /// exact same (and more than zero) incident edges. Larger signature sets may promise more accurate region similarity
    /// values but require more memory, linear in the number of hypernodes.
    template<typename Hypergraph = Mandatory>
    class NodeRegionComparator {

    private:

        using SizeCmpFunc = std::function<bool (const HyperedgeID, const HyperedgeID)>;
        using SigIterator = Array<HyperedgeID>::iterator;

    public:

        using InitiallyStablePredicate = std::function<bool (const HypernodeID)>;

        explicit NodeRegionComparator( const HyperedgeID max_signature_size = DEFAULT_SIGNATURE_SIZE) :
            _hg(nullptr),
            _is_initially_stable([](const HypernodeID) {return false;}),
            _signatures(),
            _signature_offsets(),
            _max_signature_size(max_signature_size),
            _sort_buffer_ets() {}

        // ! Calculates or recalculates the signature sets for every node of a given hypergraph in parallel.
        void calculateSignaturesParallel(const Hypergraph* hypergraph, InitiallyStablePredicate&& initially_stable_predicate = [](const HypernodeID) {return false;}) {
          ASSERT(hypergraph);
          if (_hg != hypergraph) {
            _hg = hypergraph;
          }
          ASSERT(_hg);
          _is_initially_stable = initially_stable_predicate;

          // Calculate signature size for all nodes and store current degrees, too, so they don't have to be recalculated.
          // (Initially store signature size in offsets, prefix sum afterwards.)
          _signature_offsets = vec<size_t>(_hg->initialNumNodes(), 0UL);
          vec<HyperedgeID> degrees(_hg->initialNumNodes(), ID(0));
          tbb::enumerable_thread_specific<HyperedgeID> max_degree_ets(ID(0));
          tbb::parallel_for(ID(0), _hg->initialNumNodes(), [&](const HypernodeID hn) {
              // Nodes that are initially stable will not be uncontracted any further and their region is therefore
              // never compared, so we do not need to calculate a signature for them
              if (_is_initially_stable(hn)) return;
              auto incident_nets_range = _hg->incidentEdges(hn);
              // Note: Using _hg->nodeDegree() here does not work as that returns the degree in the original graph, not the
              // current one
              const HyperedgeID degree = std::distance(incident_nets_range.begin(), incident_nets_range.end());
              degrees[hn] = degree;
              _signature_offsets[hn] = std::min(degree, _max_signature_size);
              max_degree_ets.local() = std::max(max_degree_ets.local(), degree);
          });

          HyperedgeID max_degree = max_degree_ets.combine([](HyperedgeID d1, HypernodeID d2) {
              return std::max(d1, d2);
          });

          // Compute prefix sum over signature sizes
          auto signature_sizes_prefix_sum = parallel::TBBPrefixSum<size_t>(_signature_offsets);
          tbb::parallel_scan(tbb::blocked_range<size_t>(
              0UL, UI64(_hg->initialNumNodes())), signature_sizes_prefix_sum);

          // Resize signatures array to fit sum of signature sizes
          _signatures = Array<HyperedgeID>(signature_sizes_prefix_sum.total_sum(), kInvalidHyperedge);


          _sort_buffer_ets = tbb::enumerable_thread_specific<std::vector<HyperedgeID>>(max_degree, kInvalidHyperedge);
          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
          };

          // Calculate signatures of nodes
          tbb::parallel_for(ID(0), _hg->initialNumNodes(), [&](const HypernodeID hn) {
            auto& local_sort_buffer = _sort_buffer_ets.local();
            ASSERT(sig_size(hn) <= degrees[hn]);
            calculateSignatureOfNode(hn, degrees[hn],local_sort_buffer, size_cmp);
          });
        }

        // ! Returns the approximate region similarity coefficient for two hypernodes. Runs in O(|signature|).
        double regionSimilarity(const HypernodeID hn1, const HypernodeID hn2) const {
          ASSERT(_hg);
          ASSERT(!_is_initially_stable(hn1));
          ASSERT(!_is_initially_stable(hn2));
          if (hn1 == hn2) return 1;
          IteratorRange<SigIterator> sig1 = signature(hn1);
          IteratorRange<SigIterator> sig2 = signature(hn2);
          ASSERT(sig1.begin() >= _signatures.cbegin() && sig1.end() <= _signatures.cend());
          ASSERT(sig2.begin() >= _signatures.cbegin() && sig2.end() <= _signatures.cend());
          HyperedgeID intersection_size = 0;
          HyperedgeID union_size = 0;
          SigIterator it1 = sig1.begin();
          SigIterator it2 = sig2.begin();
          while (it1 < sig1.end() && it2 < sig2.end()) {
            if (*it1 == *it2) {
              ++intersection_size;
              ++union_size;
              ++it1;
              ++it2;
            } else if (_hg->edgeSize(*it1) < _hg->edgeSize(*it2)) {
              ++union_size;
              ++it1;
            } else {
//              ASSERT(_hg->edgeSize(_signatures[sig_entry(hn1, idx1)]) >= _hg->edgeSize(_signatures[sig_entry(hn2, idx2)]));
              ++union_size;
              ++it2;
            }
          }

          while (it1 < sig1.end()) {
            ++union_size;
            ++it1;
          }

          while (it2 < sig2.end()) {
            ++union_size;
            ++it2;
          }

          if (union_size == 0) return 0;

          return (double) intersection_size / (double) union_size;
        }

//        void adaptForPinRemoval(const HypernodeID removed_pin, const HyperedgeID he) {
//          ASSERT(_hg);
//          HypernodeID smallestSizeInSignature = _hg->edgeSize(_signatures[sig_entry(removed_pin, 0)]);
//          HypernodeID sizeOfHEBeforeRemoval = _hg->edgeSize(he) + 1;
//          if (sizeOfHEBeforeRemoval >= smallestSizeInSignature) {
//            recalculateSignatureOfNode(removed_pin);
//          }
//        }
//
//        void adaptForNewPin(const HypernodeID new_pin, const HyperedgeID he) {
//          ASSERT(_hg);
//          recalculateSignatureOfNode(new_pin);
//
//          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
//              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
//              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
//          };
//
//          for (auto pin : _hg->pins(he)) {
//            if (pin == new_pin) continue;
//            ASSERT(std::is_sorted(_signatures.begin() + sig_entry(pin, 0), _signatures.begin() + _hg->nodeDegree(pin), size_cmp));
//            // If edge is still smaller than smallest edge in signature, no change
//            if (size_cmp(he, _signatures[sig_entry(pin, 0)])) continue;
//            // Else find slot where edge belongs now in signature of pin
//            HyperedgeID i = 1;
//            while (!size_cmp(he, _signatures[sig_entry(pin, i)])) {
//              ++i;
//            }
//            ASSERT(i >= 1);
//            const HyperedgeID new_slot_of_he = i - 1;
//            // Move everything left of new slot of he one slot further left until the old spot of he or until index 0 of
//            // the signature in case he was previously not in the signature
//            HyperedgeID old = _signatures[sig_entry(pin, new_slot_of_he)];
//            HyperedgeID j = new_slot_of_he;
//            while (j > 0 && old != he) {
//              HyperedgeID overwritten_value = _signatures[sig_entry(pin, j-1)];
//              _signatures[sig_entry(pin, j-1)] = old;
//              old = overwritten_value;
//              --j;
//            }
//            // Write he to its new slot
//            _signatures[sig_entry(pin, new_slot_of_he)] = he;
//            ASSERT(std::is_sorted(_signatures.begin() + sig_entry(pin, 0), _signatures.begin() + _hg->nodeDegree(pin), size_cmp));
//          }
//        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);
          parent->addChild("Signature Sets", sizeof(HyperedgeID) * _signatures.size());
        }

    private:

        void
        calculateSignatureOfNode(const HypernodeID hn, const HyperedgeID degree, std::vector<HypernodeID> &sort_buffer,
                                 SizeCmpFunc &size_cmp) {
          ASSERT(_hg);
          ASSERT(sig_size(hn) <= degree);

          auto incident_nets_range = _hg->incidentEdges(hn);

          if (degree > sort_buffer.size()) {
            sort_buffer.resize(degree, kInvalidHyperedge);
          }

          if (degree == sig_size(hn)) {
            // Simple case: Incident nets fit entirely into signature: Copy the edge ids in, sort them and return
            auto it = incident_nets_range.begin();
            for (HyperedgeID i = 0; i < sig_size(hn); ++i) {
              ASSERT(it != incident_nets_range.end());
              _signatures[sig_entry(hn, i)] = *it;
              ++it;
            }
            ASSERT(_is_initially_stable(hn) || it == incident_nets_range.end());
            std::sort(signature(hn).begin(), signature(hn).end(), size_cmp);
          } else {
            // General case: If degree > _max_signature_size find _max_signature_size largest edges incident to nodes and use those as signature
            ASSERT(sig_size(hn) == _max_signature_size);
            ASSERT(sort_buffer.size() >= degree);
            HypernodeID idx = 0;
            for (auto it = incident_nets_range.begin(); it != incident_nets_range.end(); ++it) {
              ASSERT(idx < degree);
              sort_buffer[idx] = *it;
              ++idx;
            }
            ASSERT(idx == degree);
            std::sort(sort_buffer.begin(), sort_buffer.begin() + degree, size_cmp);
            HyperedgeID signature_begin_index = degree - sig_size(hn);
            ASSERT(signature_begin_index < degree);
            for (HyperedgeID i = 0; i < sig_size(hn); ++i) {
              _signatures[sig_entry(hn, i)] = sort_buffer[signature_begin_index + i];
            }
          }
        }

//        void recalculateSignatureOfNode(const HypernodeID hn) {
//          ASSERT(_hg);
//          ASSERT(hn * _max_signature_size < _signatures.size());
//
//          SizeCmpFunc size_cmp = [&](const HyperedgeID he1, const HyperedgeID he2) {
//              // Sort primarily by edgeSize and within same sizes by ID so signatures can later be compared in O(signature size)
//              return _hg->edgeSize(he1) < _hg->edgeSize(he2) || (_hg->edgeSize(he1) == _hg->edgeSize(he2) && he1 < he2);
//          };
//
//          auto& sort_buffer = _sort_buffer_ets.local();
//          calculateSignatureOfNode(hn, sort_buffer, size_cmp);
//        }

        size_t sig_entry(const HypernodeID hn, const HypernodeID sig_idx) const {
          ASSERT(_hg);
          // _signature_offsets contains starting index of each hypernode starting at HypernodeID 1 at index 0 (so shift indices by one)
          ASSERT(hn < _signature_offsets.size());
          size_t offset;
          if (hn > 0) {
            ASSERT(sig_idx < _signature_offsets[hn] - _signature_offsets[hn - 1]);
            offset = _signature_offsets[hn - 1];
          } else {
            ASSERT(sig_idx < _signature_offsets[0UL]);
            offset = UI64(0);
          }

          size_t result = offset + sig_idx;
          ASSERT(result < _signatures.size());
          return result;
        }

        HyperedgeID sig_size(const HypernodeID hn) const {
          ASSERT(_hg);
          ASSERT(hn < _signature_offsets.size());
          if (hn > 0) {
            return _signature_offsets[hn] - _signature_offsets[hn - 1];
          } else {
            return _signature_offsets[0];
          }
        }

        IteratorRange<SigIterator> signature(const HypernodeID hn) const {
          ASSERT(_hg);
          ASSERT(hn < _signature_offsets.size());
          size_t begin_offset;
          if (hn > 0) {
            begin_offset = _signature_offsets[hn - 1];
          } else {
            begin_offset = UI64(0);
          }

          if (hn == _signature_offsets.size() - 1) {
            ASSERT(_signature_offsets[hn] == _signatures.size());
            return IteratorRange<SigIterator>(_signatures.cbegin() + begin_offset, _signatures.cend());
          } else {
            size_t end_offset = _signature_offsets[hn];
            return IteratorRange<SigIterator>(_signatures.cbegin() + begin_offset,
                                                                     _signatures.cbegin() + end_offset);
          }
        }

        const Hypergraph* _hg;
        InitiallyStablePredicate _is_initially_stable;

        parallel::scalable_vector<size_t> _signature_offsets;
        Array<HyperedgeID> _signatures;
        const HyperedgeID _max_signature_size;
        tbb::enumerable_thread_specific<std::vector<HyperedgeID>> _sort_buffer_ets;

        static constexpr HyperedgeID DEFAULT_SIGNATURE_SIZE = 10;

    };

} // end namespace mt_kahypar::ds
