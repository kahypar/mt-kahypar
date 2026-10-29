/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren
 * Copyright (C) 2019 Tobias Heuer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <vector>
#include <functional>
#include <cstdint>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <limits>
#include "include/mtkahypartypes.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/datastructures/two_level_vector.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/bit_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

// Forward declarations
class CompressedHypergraphFactory;
template <typename Hypergraph,
          typename ConnectivityInformation>
class PartitionedHypergraph;

// ####################### Compression Utility Functions #######################

inline void encode_varint(uint64_t value, std::vector<uint8_t> &out) {
  while (value >= 0x80) {
    out.push_back(static_cast<uint8_t>(value & 0x7F) | 0x80);
    value >>= 7;
  }
  out.push_back(static_cast<uint8_t>(value & 0x7F));
}

inline uint64_t decode_varint(const std::vector<uint8_t> &in, size_t &pos) {
  uint64_t result = 0;
  int shift = 0;
  while (pos < in.size()) {
    const uint8_t byte = in[pos++];
    result |= (static_cast<uint64_t>(byte & 0x7F) << shift);
    if ((byte & 0x80) == 0) break;
    shift += 7;
  }
  return result;
}

inline uint64_t decode_varint_bounded(const std::vector<uint8_t>& in, size_t& pos, const size_t limit) {
  uint64_t result = 0;
  int shift = 0;
  while (pos < limit) {
    const uint8_t byte = in[pos++];
    result |= (static_cast<uint64_t>(byte & 0x7F) << shift);
    if ((byte & 0x80) == 0) break;   // last byte
    shift += 7;
  }
  return result;
}

// Sort and deduplicate before gap encoding to guarantee a monotone sequence
inline void append_compressed_sequence(const std::vector<HypernodeID> &sequence, std::vector<uint8_t> &target) {
  if (sequence.empty()) return;
  std::vector<HypernodeID> sorted = sequence;
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  HypernodeID prev = 0;
  for (const HypernodeID node : sorted) {
    ASSERT(node >= prev);
    const HypernodeID gap = node - prev;
    encode_varint(gap, target);
    prev = node;
  }
}

inline std::vector<uint8_t> compress_sequence(const std::vector<HypernodeID> &sequence) {
  std::vector<uint8_t> encoded;
  append_compressed_sequence(sequence, encoded);
  return encoded;
}

class CompressedHypergraph {
    static constexpr bool enable_heavy_assert = false;
    static constexpr HyperedgeID HIGH_DEGREE_CONTRACTION_THRESHOLD = ID(500000);
    
    static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
    static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");
    
        using UncontractionFunction = std::function<void(const HypernodeID, const HypernodeID, const HyperedgeID)>;

    // Used for the two-level vectors below
    using WeightBaseT = uint16_t;

#define NOOP_BATCH_FUNC [] (const HypernodeID, const HypernodeID, const HyperedgeID) { }

    // Proxy for a hypernode backed by parent offsets/enabled arrays
    class Hypernode {
    public:
        using IDType = HypernodeID;
        Hypernode() : _hg(nullptr), _id(0) {}
        Hypernode(CompressedHypergraph* hg, HypernodeID id) : _hg(hg), _id(id) {}

    bool isDisabled() const { return !_hg || !_hg->_hypernode_enabled[_id]; }
    void enable() const { ASSERT(_hg); _hg->_hypernode_enabled[_id] = true; }
    void disable() const { ASSERT(_hg); _hg->_hypernode_enabled[_id] = false; }

    // Start of incident list (after degree header)
    size_t firstEntry() const {
        const size_t begin = _hg->_hypernode_offsets[_id];
        const size_t raw_end = nextNodePos();
        size_t pos = begin;
        decode_varint_bounded(_hg->_compressed_incident_nets, pos, raw_end);
        return pos;
    }
    void setFirstEntry(const size_t begin) const {
        ASSERT(!isDisabled());
        _hg->_hypernode_offsets[_id] = begin;
    }

    // End position (one past last data byte), trim 0x00 padding bytes if present (created upon edge removal)
    size_t firstInvalidEntry() const {
        const size_t start = firstEntry();
        size_t end = nextNodePos();
        while (
            end > (start + 1)
            && _hg->_compressed_incident_nets[end - 1] == 0x00u
            && (_hg->_compressed_incident_nets[end - 2] & 0x80) == 0
        ) end--;
        if (end == start + 1 && this->degree() == 0) return start; // edge case: only padding bytes
        return end;
    }
    HyperedgeID degree() const {
        const size_t begin = _hg->_hypernode_offsets[_id];
        size_t pos = begin;
        return static_cast<HyperedgeID>(decode_varint_bounded(_hg->_compressed_incident_nets, pos, nextNodePos()));
    }
        HypernodeWeight weight() const {
            // Two-level lazy vector: default is 1 when not initialized
            return _hg->_hypernode_weights[_id];
        }
        void setWeight(HypernodeWeight w) const {
            if (w == 1 && _hg->_hypernode_weights.empty()) return; // keep lazy default
            if (_hg->_hypernode_weights.empty()) _hg->_hypernode_weights.ensure_initialized(_hg->_num_hypernodes);
            _hg->_hypernode_weights[_id] = w;
        }
        HypernodeID id() const { return _id; }

    private:
        CompressedHypergraph* _hg;
        HypernodeID _id;
        size_t nextNodePos() const {
            return (_id + 1 < _hg->_num_hypernodes)
                ? _hg->_hypernode_offsets[_id + 1]
                : _hg->_compressed_incident_nets.size();
        }
    };

    // Proxy for a hyperedge
    class Hyperedge {
    public:
        using IDType = HyperedgeID;
        Hyperedge() : _hg(nullptr), _id(0) {}
        Hyperedge(CompressedHypergraph* hg, HyperedgeID id) : _hg(hg), _id(id) {}

        bool isDisabled() const { return !_hg || !_hg->_hyperedge_enabled[_id]; }
        void enable() const { ASSERT(_hg); _hg->_hyperedge_enabled[_id] = true; }
        void disable() const { ASSERT(_hg); _hg->_hyperedge_enabled[_id] = false; }

        // Start of pin list (after size header)
        size_t firstEntry() const {
            const size_t begin = _hg->_hyperedge_offsets[_id];
            size_t pos = begin;
            (void)decode_varint_bounded(_hg->_compressed_incidence_array, pos, nextEdgePos());
            return pos;
        }
        void setFirstEntry(const size_t begin) const { ASSERT(!isDisabled()); _hg->_hyperedge_offsets[_id] = begin; }

        // End position (one past last data byte), trim a single trailing 0x00 padding byte if present
        size_t firstInvalidEntry() const {
            const size_t start = firstEntry();
            size_t end = nextEdgePos();
            while (
                end > (start + 1)
                && _hg->_compressed_incidence_array[end - 1] == 0x00u
                && (_hg->_compressed_incidence_array[end - 2] & 0x80) == 0
            ) end--;
            if (end == start + 1 && this->size() == 0) return start; // edge case: only padding bytes
            return end;
        }
        HypernodeID size() const {
            const size_t begin = _hg->_hyperedge_offsets[_id];
            size_t pos = begin;
            return static_cast<HypernodeID>(decode_varint_bounded(_hg->_compressed_incidence_array, pos, nextEdgePos()));
        }
        HyperedgeWeight weight() const { return _hg->_hyperedge_weights[_id]; }
        void setWeight(HyperedgeWeight w) const {
            if (w == 1 && _hg->_hyperedge_weights.empty()) return;
            if (_hg->_hyperedge_weights.empty()) _hg->_hyperedge_weights.ensure_initialized(_hg->_num_hyperedges);
            _hg->_hyperedge_weights[_id] = w;
        }
        HyperedgeID id() const { return _id; }

        bool operator==(const Hyperedge& other) const { return _hg == other._hg && _id == other._id; }
        bool operator!=(const Hyperedge& other) const { return !(*this == other); }

    private:
        CompressedHypergraph* _hg;
        HyperedgeID _id;
        size_t nextEdgePos() const {
            return (_id + 1 < _hg->_num_hyperedges)
                ? _hg->_hyperedge_offsets[_id + 1]
                : _hg->_compressed_incidence_array.size();
        }
    };

    /**
     * Compressed Iterator for incident nets and pins
     * Decompresses data on-the-fly during iteration
     */
    template<typename ElementType>
    class CompressedIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = ElementType;
        using reference = ElementType;
        using pointer = const ElementType*;
        using difference_type = std::ptrdiff_t;

        CompressedIterator(const std::vector<uint8_t>& compressed_data, 
                          size_t start_pos, size_t end_pos) :
            _compressed_data(compressed_data),
            _start_pos(start_pos),
            _end_pos(end_pos),
            _current_pos(start_pos),
            _current_value(0),
            _accumulated_value(0),
            _has_value(false) {
            if (_current_pos < _end_pos) {
                advance();
                _has_value = true;
            }
        }

        // End iterator
        CompressedIterator(const std::vector<uint8_t>& compressed_data, size_t end_pos) :
            _compressed_data(compressed_data),
            _start_pos(end_pos),
            _end_pos(end_pos),
            _current_pos(end_pos),
            _current_value(0),
            _accumulated_value(0),
            _has_value(false) {}

        ElementType operator*() const {
            ASSERT(_has_value);
            return _current_value;
        }

        CompressedIterator& operator++() {
            if (_has_value) advance();
            return *this;
        }

        CompressedIterator operator++(int) {
            CompressedIterator copy = *this;
            ++(*this);
            return copy;
        }

        // Compare iterator state by position/bounds to be robust
        bool operator!=(const CompressedIterator& rhs) const {
            // Both at end
            if (!_has_value && !rhs._has_value) return false;
            return _current_pos != rhs._current_pos || _end_pos != rhs._end_pos || _has_value != rhs._has_value;
        }

        bool operator==(const CompressedIterator& rhs) const {
            return !(*this != rhs);
        }

    private:
        void advance() {
            if (_current_pos >= _end_pos) {
                _has_value = false;
                return;
            };
            const uint64_t gap = decode_varint_bounded(_compressed_data, _current_pos, _end_pos);
            _accumulated_value += static_cast<ElementType>(gap);
            _current_value = _accumulated_value;
        }

        const std::vector<uint8_t>& _compressed_data;
        size_t _start_pos;
        size_t _end_pos;
        size_t _current_pos;
        ElementType _current_value;
        ElementType _accumulated_value;
        bool _has_value;
    };

    // Iterators over IDs using enabled flags
    template <typename ID>
    class IDIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = ID;
        using reference = ID;
        using pointer = void;
        using difference_type = std::ptrdiff_t;

        IDIterator(const CompressedHypergraph* hg, ID id, ID max, bool edge) :
            _hg(hg), _id(id), _max(max), _edge(edge) {
            if (_id != _max && !isEnabled(_id)) ++(*this);
        }

        ID operator*() const { return _id; }

        IDIterator& operator++() {
            do { ++_id; } while (_id < _max && !isEnabled(_id));
            return *this;
        }

        IDIterator operator++(int) { auto c = *this; ++(*this); return c; }
        bool operator!=(const IDIterator& rhs) const { return _id != rhs._id; }
        bool operator==(const IDIterator& rhs) const { return _id == rhs._id; }

    private:
        bool isEnabled(ID x) const {
            return _edge ? _hg->_hyperedge_enabled[x]
                         : _hg->_hypernode_enabled[x];
        }
        const CompressedHypergraph* _hg;
        ID _id;
        ID _max;
        bool _edge;
    };

    using CompressedIncidenceArray = std::vector<uint8_t>;    // Compressed pins storage
    using CompressedIncidentNets = std::vector<uint8_t>;      // Compressed incident nets storage

    // Contains buffers needed during contractions
    struct TmpContractionBuffer {
        explicit TmpContractionBuffer(const HypernodeID num_hypernodes,
                                      const HyperedgeID num_hyperedges,
                                      const HyperedgeID num_pins) {
            // Allocate through memory pool-backed Arrays to reduce peak RSS and fragmentation
            tbb::parallel_invoke([&] {
                mapping.resize("Coarsening", "mapping", num_hypernodes);
            }, [&] {
                tmp_hypernodes.resize("Coarsening", "tmp_hypernodes", num_hypernodes);
            }, [&] {
                tmp_incident_nets.resize("Coarsening", "tmp_incident_nets", num_pins);
            }, [&] {
                tmp_num_incident_nets.resize("Coarsening", "tmp_num_incident_nets", num_hypernodes);
            }, [&] {
                hn_weights.resize("Coarsening", "hn_weights", num_hypernodes);
            }, [&] {
                tmp_hyperedges.resize("Coarsening", "tmp_hyperedges", num_hyperedges);
            }, [&] {
                tmp_incidence_array.resize("Coarsening", "tmp_incidence_array", num_pins);
            }, [&] {
                he_sizes.resize("Coarsening", "he_sizes", num_hyperedges);
            }, [&] {
                valid_hyperedges.resize("Coarsening", "valid_hyperedges", num_hyperedges);
            });
        }

        ds::Array<size_t> mapping;
        ds::Array<Hypernode> tmp_hypernodes;
        // For incident nets/pins during contraction we need uncompressed ID buffers
        ds::Array<HyperedgeID> tmp_incident_nets;
        ds::Array<parallel::IntegralAtomicWrapper<size_t>> tmp_num_incident_nets;
        ds::Array<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
        ds::Array<Hyperedge> tmp_hyperedges;
        ds::Array<HypernodeID> tmp_incidence_array;
        ds::Array<size_t> he_sizes;
        ds::Array<size_t> valid_hyperedges;
    };

public:
    static constexpr bool is_graph = false;
    static constexpr bool is_static_hypergraph = true;
    static constexpr bool is_partitioned = false;
    static constexpr size_t SIZE_OF_HYPERNODE = sizeof(Hypernode);
    static constexpr size_t SIZE_OF_HYPEREDGE = sizeof(Hyperedge);
    static constexpr mt_kahypar_hypergraph_type_t TYPE = COMPRESSED_HYPERGRAPH;

    // Factory
    using Factory = CompressedHypergraphFactory;

    // Iterator types (ID-based over enabled flags)
    // Note: legacy element iterator from other hypergraph variants is not used here
    // Pins of a hyperedge are HypernodeID (degree-bounded, skip zero-gap padding bytes)
    class PinsIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = HypernodeID;
        using reference = HypernodeID;
        using pointer = const HypernodeID*;
        using difference_type = std::ptrdiff_t;

        PinsIterator(const std::vector<uint8_t>& data,
                     size_t start_pos,
                     size_t end_pos,
                     size_t remaining)
            : _data(data), _end(end_pos), _pos(start_pos), _acc(0), _cur(0), _has(false), _left(remaining), _emitted(0) {
            advance_to_valid();
        }

        PinsIterator(const std::vector<uint8_t>& data, size_t end_pos)
            : _data(data), _end(end_pos), _pos(end_pos), _acc(0), _cur(0), _has(false), _left(0), _emitted(0) {}

        HypernodeID operator*() const { ASSERT(_has); return _cur; }
        PinsIterator& operator++() { if (_has) advance_to_valid(); return *this; }
        PinsIterator operator++(int) { auto c = *this; ++(*this); return c; }
        bool operator!=(const PinsIterator& rhs) const { if (!_has && !rhs._has) return false; return _pos != rhs._pos || _end != rhs._end || _has != rhs._has || _left != rhs._left; }
        bool operator==(const PinsIterator& rhs) const { return !(*this != rhs); }

    private:
        void advance_to_valid() {
            _has = false;
            while (_pos < _end && _left > 0) {
                const uint64_t gap = decode_varint_bounded(_data, _pos, _end);
                // skip zero-gap placeholders except possibly the very first emitted element
                if (_emitted > 0 && gap == 0) { continue; }
                _acc += static_cast<HypernodeID>(gap);
                _cur = _acc;
                --_left;
                ++_emitted;
                _has = true;
                return;
            }
            _has = false;
        }

        const std::vector<uint8_t>& _data;
        size_t _end;
        size_t _pos;
        HypernodeID _acc;
        HypernodeID _cur;
        bool _has;
        size_t _left;
        size_t _emitted;
    };
        using IncidenceIterator = PinsIterator;
        // Incident nets of a hypernode are HyperedgeID
        // We wrap decoding with a filtering iterator to skip out-of-range IDs defensively.
    class IncidentNetsFilteringIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = HyperedgeID;
        using reference = HyperedgeID;
        using pointer = const HyperedgeID*;
        using difference_type = std::ptrdiff_t;

        IncidentNetsFilteringIterator(const std::vector<uint8_t>& data,
                                      size_t start_pos,
                                      size_t end_pos,
                                      HyperedgeID max_edge_id,
                                      const BitVector* enabled,
                                      size_t remaining)
            : _data(data), _end(end_pos), _pos(start_pos),
              _acc(0), _cur(0), _has(false), _max(max_edge_id), _enabled(enabled), _left(remaining), _emitted(0) {
            advance_to_valid();
        }

        // End iterator
        IncidentNetsFilteringIterator(const std::vector<uint8_t>& data,
                                      size_t end_pos)
            : _data(data), _end(end_pos), _pos(end_pos),
              _acc(0), _cur(0), _has(false), _max(0), _enabled(nullptr), _left(0), _emitted(0) { }

        HyperedgeID operator*() const { ASSERT(_has); return _cur; }

        IncidentNetsFilteringIterator& operator++() { advance_to_valid(); return *this; }

        IncidentNetsFilteringIterator operator++(int) { auto copy = *this; ++(*this); return copy; }

        bool operator!=(const IncidentNetsFilteringIterator& rhs) const {
            if (!_has && !rhs._has) return false;
            return _pos != rhs._pos || _end != rhs._end || _has != rhs._has;
        }

        bool operator==(const IncidentNetsFilteringIterator& rhs) const { return !(*this != rhs); }

    private:
        void advance_to_valid() {
            _has = false;
            while (_pos < _end && _left > 0) {
                const uint64_t gap = decode_varint_bounded(_data, _pos, _end);
                // skip zero-gap placeholders except possibly the very first emitted element
                if (_emitted > 0 && gap == 0) { continue; }
                _acc += static_cast<HyperedgeID>(gap);
                if (_acc < _max && (_enabled == nullptr || (*_enabled)[_acc])) {
                    _cur = _acc;
                    _has = true;
                    --_left;
                    ++_emitted;
                    return;
                }
                // else: skip invalid ID and continue without consuming logical count
            }
            // exhausted
            _has = false;
        }

        const std::vector<uint8_t>& _data;
        size_t _end;
        size_t _pos;
        HyperedgeID _acc;
        HyperedgeID _cur;
        bool _has;
        HyperedgeID _max;
        const BitVector* _enabled;
        size_t _left;
        size_t _emitted;
    };

        using IncidentNetsIterator = IncidentNetsFilteringIterator;

    struct ParallelHyperedge {
        HyperedgeID removed_hyperedge;
        HyperedgeID representative;
    };

    explicit CompressedHypergraph() :
        _num_hypernodes(0),
        _num_removed_hypernodes(0),
        _removed_degree_zero_hn_weight(0),
        _num_hyperedges(0),
        _num_removed_hyperedges(0),
        _max_edge_size(0),
        _num_pins(0),
        _total_degree(0),
        _total_weight(0),
        _hypernode_offsets(),
        _hypernode_enabled(),
        _compressed_incident_nets(),
        _hyperedge_offsets(),
        _hyperedge_enabled(),
        _compressed_incidence_array(),
        _hypernode_weights(),
        _hyperedge_weights(),
        _community_ids(0),
        _fixed_vertices(),
        _tmp_contraction_buffer(nullptr) { }

    explicit CompressedHypergraph(const std::string& filename,
        bool remove_single_pin_hes = true);
        
    CompressedHypergraph(const CompressedHypergraph&) = delete;
    CompressedHypergraph& operator=(const CompressedHypergraph&) = delete;

    CompressedHypergraph(CompressedHypergraph&& other) :
        _num_hypernodes(other._num_hypernodes),
        _num_removed_hypernodes(other._num_removed_hypernodes),
        _removed_degree_zero_hn_weight(other._removed_degree_zero_hn_weight),
        _num_hyperedges(other._num_hyperedges),
        _num_removed_hyperedges(other._num_removed_hyperedges),
        _max_edge_size(other._max_edge_size),
        _num_pins(other._num_pins),
        _total_degree(other._total_degree),
        _total_weight(other._total_weight),
        _hypernode_offsets(std::move(other._hypernode_offsets)),
        _hypernode_enabled(std::move(other._hypernode_enabled)),
        _compressed_incident_nets(std::move(other._compressed_incident_nets)),
        _hyperedge_offsets(std::move(other._hyperedge_offsets)),
        _hyperedge_enabled(std::move(other._hyperedge_enabled)),
        _compressed_incidence_array(std::move(other._compressed_incidence_array)),
        _hypernode_weights(std::move(other._hypernode_weights)),
        _hyperedge_weights(std::move(other._hyperedge_weights)),
        _community_ids(std::move(other._community_ids)),
        _fixed_vertices(std::move(other._fixed_vertices)),
        _tmp_contraction_buffer(std::move(other._tmp_contraction_buffer)) {
        _fixed_vertices.setHypergraph(this);
        other._tmp_contraction_buffer = nullptr;
    }

    CompressedHypergraph& operator=(CompressedHypergraph&& other) {
        if (this != &other) {
            _num_hypernodes = other._num_hypernodes;
            _num_removed_hypernodes = other._num_removed_hypernodes;
            _removed_degree_zero_hn_weight = other._removed_degree_zero_hn_weight;
            _num_hyperedges = other._num_hyperedges;
            _num_removed_hyperedges = other._num_removed_hyperedges;
            _max_edge_size = other._max_edge_size;
            _num_pins = other._num_pins;
            _total_degree = other._total_degree;
            _total_weight = other._total_weight;
            _hypernode_offsets = std::move(other._hypernode_offsets);
            _hypernode_enabled = std::move(other._hypernode_enabled);
            _compressed_incident_nets = std::move(other._compressed_incident_nets);
            _hyperedge_offsets = std::move(other._hyperedge_offsets);
            _hyperedge_enabled = std::move(other._hyperedge_enabled);
            _compressed_incidence_array = std::move(other._compressed_incidence_array);
            _hypernode_weights = std::move(other._hypernode_weights);
            _hyperedge_weights = std::move(other._hyperedge_weights);
            _community_ids = std::move(other._community_ids);
            _fixed_vertices = std::move(other._fixed_vertices);
            _fixed_vertices.setHypergraph(this);
            _tmp_contraction_buffer = std::move(other._tmp_contraction_buffer);
            other._tmp_contraction_buffer = nullptr;
        }
        return *this;
    }

    ~CompressedHypergraph() {
        freeInternalData();
    }

    // ####################### General Hypergraph Stats #######################
    
    HypernodeID initialNumNodes() const {
        return _num_hypernodes;
    }

    HypernodeID numRemovedHypernodes() const {
        return _num_removed_hypernodes;
    }

    HypernodeWeight weightOfRemovedDegreeZeroVertices() const {
        return _removed_degree_zero_hn_weight;
    }

    HyperedgeID initialNumEdges() const {
        return _num_hyperedges;
    }

    HyperedgeID numRemovedHyperedges() const {
        return _num_removed_hyperedges;
    }

    void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
        _num_removed_hyperedges = num_removed_hyperedges;
    }

    HypernodeID initialNumPins() const {
        return _num_pins;
    }

    HypernodeID initialTotalVertexDegree() const {
        return _total_degree;
    }

    HypernodeWeight totalWeight() const {
        return _total_weight;
    }

    void computeAndSetTotalNodeWeight(parallel_tag_t);


    // ####################### Iterators #######################

    template<typename F>
    void doParallelForAllNodes(const F& f) const {
        // Sequential version (parallel logic omitted as requested)
        for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
            if (nodeIsEnabled(hn)) {
                f(hn);
            }
        }
    }

    template<typename F>
    void doParallelForAllEdges(const F& f) const {
        // Sequential version (parallel logic omitted as requested)
        for (HyperedgeID he = 0; he < _num_hyperedges; ++he) {
            if (edgeIsEnabled(he)) {
                f(he);
            }
        }
    }

    using HypernodeIterator = IDIterator<HypernodeID>;
    using HyperedgeIterator = IDIterator<HyperedgeID>;

    IteratorRange<HypernodeIterator> nodes() const {
        return IteratorRange<HypernodeIterator>(
            HypernodeIterator(this, ID(0), _num_hypernodes, /*edge*/false),
            HypernodeIterator(this, _num_hypernodes, _num_hypernodes, /*edge*/false));
    }

    IteratorRange<HyperedgeIterator> edges() const {
        return IteratorRange<HyperedgeIterator>(
            HyperedgeIterator(this, ID(0), _num_hyperedges, /*edge*/true),
            HyperedgeIterator(this, _num_hyperedges, _num_hyperedges, /*edge*/true));
    }

    IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
        ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
        const auto hn = hypernode(u);
        const size_t begin = hn.firstEntry();
        const size_t end = hn.firstInvalidEntry();
        const size_t deg = static_cast<size_t>(hn.degree());
        return IteratorRange<IncidentNetsIterator>(
            IncidentNetsIterator(_compressed_incident_nets, begin, end, _num_hyperedges, &_hyperedge_enabled, deg),
            IncidentNetsIterator(_compressed_incident_nets, end));
    }

    IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
        ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
        const auto he = hyperedge(e);
        const size_t begin = he.firstEntry();
        const size_t end = he.firstInvalidEntry();
        const size_t esize = static_cast<size_t>(he.size());

        return IteratorRange<IncidenceIterator>(
            IncidenceIterator(_compressed_incidence_array, begin, end, esize),
            IncidenceIterator(_compressed_incidence_array, end));
    }

    // ####################### Hypernode Information #######################

    HypernodeWeight nodeWeight(const HypernodeID u) const { return _hypernode_weights[u]; }

    void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
        ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
        // Read old value first (works with lazy default = 1)
    const HypernodeWeight old_w = _hypernode_weights[u];
        if (_hypernode_weights.empty() && weight == 1) return; // keep lazy default
        if (_hypernode_weights.empty()) _hypernode_weights.ensure_initialized(_num_hypernodes);
        _hypernode_weights[u] = weight;
        // Maintain total weight incrementally
        _total_weight += (weight - old_w);
    }

    HyperedgeID nodeDegree(const HypernodeID u) const {
        ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
        return hypernode(u).degree();
    }

    bool nodeIsEnabled(const HypernodeID u) const { return _hypernode_enabled[u]; }
    void enableHypernode(const HypernodeID u) { _hypernode_enabled[u] = true; }
    void disableHypernode(const HypernodeID u) { _hypernode_enabled[u] = false; }

    void removeHypernode(const HypernodeID u) {
        hypernode(u).disable();
        ++_num_removed_hypernodes;
    }

    void removeDegreeZeroHypernode(const HypernodeID u) {
        ASSERT(nodeDegree(u) == 0);
        _removed_degree_zero_hn_weight += nodeWeight(u);
        removeHypernode(u);
    }

    void restoreDegreeZeroHypernode(const HypernodeID u) {
        hypernode(u).enable();
        ASSERT(nodeDegree(u) == 0);
        _removed_degree_zero_hn_weight -= nodeWeight(u);
    }

    // ####################### Hyperedge Information #######################

    HypernodeWeight edgeWeight(const HyperedgeID e) const {
        ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    return _hyperedge_weights[e];
    }

    void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
        ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    if (_hyperedge_weights.empty() && weight == 1) return;
    if (_hyperedge_weights.empty()) _hyperedge_weights.ensure_initialized(_num_hyperedges);
    _hyperedge_weights[e] = weight;
    }

    HypernodeID edgeSize(const HyperedgeID e) const {
        ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
        return hyperedge(e).size();
    }

    HypernodeID maxEdgeSize() const {
        return _max_edge_size;
    }

    bool edgeIsEnabled(const HyperedgeID e) const { return _hyperedge_enabled[e]; }
    void enableHyperedge(const HyperedgeID e) { _hyperedge_enabled[e] = true; }
    void disableHyperedge(const HyperedgeID e) { _hyperedge_enabled[e] = false; }

    PartitionID communityID(const HypernodeID u) const {
        return _community_ids[u];
    }

    void setCommunityID(const HypernodeID u, const PartitionID community_id) {
        _community_ids[u] = community_id;
    }

    // ####################### Fixed Vertex Support #######################

    void addFixedVertexSupport(FixedVertexSupport<CompressedHypergraph>&& fixed_vertices) {
        _fixed_vertices = std::move(fixed_vertices);
        _fixed_vertices.setHypergraph(this);
    }

    bool hasFixedVertices() const {
        return _fixed_vertices.hasFixedVertices();
    }

    HypernodeWeight totalFixedVertexWeight() const {
        return _fixed_vertices.totalFixedVertexWeight();
    }

    HypernodeWeight fixedVertexBlockWeight(const PartitionID block) const {
        return _fixed_vertices.fixedVertexBlockWeight(block);
    }

    bool isFixed(const HypernodeID hn) const {
        return _fixed_vertices.isFixed(hn);
    }

    PartitionID fixedVertexBlock(const HypernodeID hn) const {
        return _fixed_vertices.fixedVertexBlock(hn);
    }

    void setMaxFixedVertexBlockWeight(const std::vector<HypernodeWeight>& max_block_weights) {
        _fixed_vertices.setMaxBlockWeight(max_block_weights);
    }

    const FixedVertexSupport<CompressedHypergraph>& fixedVertexSupport() const {
        return _fixed_vertices;
    }

    FixedVertexSupport<CompressedHypergraph> copyOfFixedVertexSupport() const {
        return _fixed_vertices.copy();
    }

    // ####################### Contract / Uncontract #######################

    // COMPRESSION_TODO
    CompressedHypergraph contract(parallel::scalable_vector<HypernodeID>& communities, bool deterministic = false);

    bool registerContraction(const HypernodeID, const HypernodeID) {
        throw UnsupportedOperationException(
            "registerContraction(u, v) is not supported in compressed hypergraph");
        return false;
    }

    size_t contract(const HypernodeID,
                    const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
        unused(max_node_weight);
        throw UnsupportedOperationException(
        "contract(v, max_node_weight) is not supported in compressed hypergraph");
        return 0;
    }

    void uncontract(const Batch&,
                    const UncontractionFunction& case_one_func = NOOP_BATCH_FUNC,
                    const UncontractionFunction& case_two_func = NOOP_BATCH_FUNC) {
        unused(case_one_func);
        unused(case_two_func);
        throw UnsupportedOperationException(
        "uncontract(batch) is not supported in compressed hypergraph");
    }

    VersionedBatchVector createBatchUncontractionHierarchy(const size_t) {
        throw UnsupportedOperationException(
            "createBatchUncontractionHierarchy(batch_size) is not supported in compressed hypergraph");
        return {};
    }

    // ####################### Remove / Restore Hyperedges #######################

    // To avoid rebuilding the entire _compressed_incident_nets, 0x00 padding bytes are left in place.
    void removeEdge(const HyperedgeID he) {
        ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");

        const auto HE = hyperedge(he);
        size_t he_pos = HE.firstEntry();
        const size_t he_end = HE.firstInvalidEntry();
        std::vector<HypernodeID> pins_of_he;
        const size_t esize = static_cast<size_t>(HE.size());
        pins_of_he.reserve(esize);
        {
            HypernodeID acc = 0;
            size_t emitted = 0;
            while (he_pos < he_end && emitted < esize) {
                const uint64_t gap = decode_varint_bounded(_compressed_incidence_array, he_pos, he_end);
                if (emitted > 0 && gap == 0) { continue; }
                acc += static_cast<HypernodeID>(gap);
                pins_of_he.push_back(acc);
                ++emitted;
            }
        }

        // For each pin, compute new bytes and overwrite in-place
        struct PendingUpdate { HypernodeID u; std::vector<uint8_t> bytes; size_t begin; };
        std::vector<PendingUpdate> pending; pending.reserve(pins_of_he.size());
        for (const HypernodeID u : pins_of_he) {
            if (!nodeIsEnabled(u)) continue;
            const auto HN = hypernode(u);
            const size_t header_begin = _hypernode_offsets[u];
            const size_t begin = HN.firstEntry();
            const size_t end = HN.firstInvalidEntry();
            // Decode current incidents
            const size_t expected = static_cast<size_t>(HN.degree());
            std::vector<HyperedgeID> inc; inc.reserve(expected);
            size_t cur = begin; HyperedgeID acc = 0; size_t emitted = 0;
            while (cur < end && emitted < expected) {
                const uint64_t gap = decode_varint_bounded(_compressed_incident_nets, cur, end);
                if (emitted > 0 && gap == 0) { continue; }
                acc += static_cast<HyperedgeID>(gap);
                if (acc != he && acc < _num_hyperedges) inc.push_back(acc);
                ++emitted;
            }
            // Re-encode
            std::vector<uint8_t> bytes;
            std::sort(inc.begin(), inc.end()); inc.erase(std::unique(inc.begin(), inc.end()), inc.end());
            encode_varint(static_cast<uint64_t>(inc.size()), bytes); // header
            HyperedgeID prev = 0; for (HyperedgeID x : inc) { ASSERT(x >= prev); encode_varint(static_cast<uint64_t>(x - prev), bytes); prev = x; }
            pending.push_back(PendingUpdate{u, std::move(bytes), header_begin});
        }

        // In-place overwrite; we intentionally leave any trailing unused bytes as-is
        for (auto& up : pending) {
            if (!up.bytes.empty()) {
                std::copy(up.bytes.begin(), up.bytes.end(), _compressed_incident_nets.begin() + up.begin);
            }
        }

        ++_num_removed_hyperedges;
        disableHyperedge(he);
    }
    
    // To avoid rebuilding the entire _compressed_incident_nets, 0x00 padding bytes are left in place.
    void removeLargeEdge(const HyperedgeID he) {
        ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");

        const auto HE = hyperedge(he);
        size_t he_pos = HE.firstEntry();
        const size_t he_end = HE.firstInvalidEntry();
        std::vector<HypernodeID> pins_of_he;
        const size_t esize = static_cast<size_t>(HE.size());
        pins_of_he.reserve(esize);
        {
            HypernodeID acc = 0;
            size_t emitted = 0;
            while (he_pos < he_end && emitted < esize) {
                const uint64_t gap = decode_varint_bounded(_compressed_incidence_array, he_pos, he_end);
                if (emitted > 0 && gap == 0) { continue; }
                acc += static_cast<HypernodeID>(gap);
                pins_of_he.push_back(acc);
                ++emitted;
            }
        }

        // Compute new bytes for each affected node and overwrite in-place
        struct PendingUpdateLE { HypernodeID u; std::vector<uint8_t> bytes; size_t begin; };
        std::vector<PendingUpdateLE> updates; updates.reserve(pins_of_he.size());
        for (const HypernodeID u : pins_of_he) {
            if (!nodeIsEnabled(u)) continue;
            const auto HN = hypernode(u);
            const size_t header_begin = _hypernode_offsets[u];
            const size_t begin = HN.firstEntry();
            const size_t end = HN.firstInvalidEntry();
            // Decode incidents
            const size_t expected = static_cast<size_t>(HN.degree());
            std::vector<HyperedgeID> inc; inc.reserve(expected);
            size_t cur = begin; HyperedgeID acc = 0; size_t emitted = 0;
            while (cur < end && emitted < expected) {
                const uint64_t gap = decode_varint_bounded(_compressed_incident_nets, cur, end);
                if (emitted > 0 && gap == 0) { continue; }
                acc += static_cast<HyperedgeID>(gap);
                if (acc != he && acc < _num_hyperedges) inc.push_back(acc);
                ++emitted;
            }
            std::vector<uint8_t> bytes;
            std::sort(inc.begin(), inc.end()); inc.erase(std::unique(inc.begin(), inc.end()), inc.end());
            encode_varint(static_cast<uint64_t>(inc.size()), bytes);
            HyperedgeID prev = 0; for (HyperedgeID x : inc) { ASSERT(x >= prev); encode_varint(static_cast<uint64_t>(x - prev), bytes); prev = x; }
            updates.push_back(PendingUpdateLE{u, std::move(bytes), header_begin});
        }
        for (auto& up : updates) {
            if (!up.bytes.empty()) {
                std::copy(up.bytes.begin(), up.bytes.end(), _compressed_incident_nets.begin() + up.begin);
            }
        }

        // Mirror static semantics: do not change removed count
        disableHyperedge(he);
    }

    void restoreLargeEdge(const HyperedgeID& he) {
        ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "is enabled");

        // Decode pins of the hyperedge
        const auto HE = hyperedge(he);
        size_t he_pos = HE.firstEntry();
        const size_t he_end = HE.firstInvalidEntry();
        std::vector<HypernodeID> pins_of_he;
        pins_of_he.reserve(HE.size());
        {
            HypernodeID acc = 0;
            size_t emitted = 0;
            const size_t expected = static_cast<size_t>(HE.size());
            while (he_pos < he_end && emitted < expected) {
                const uint64_t gap = decode_varint_bounded(_compressed_incidence_array, he_pos, he_end);
                if (emitted > 0 && gap == 0) { continue; }
                acc += static_cast<HypernodeID>(gap);
                pins_of_he.push_back(acc);
                ++emitted;
            }
        }

        // First pass: compute new bytes per affected node and check if in-place growth fits
        struct PendingUpdate { HypernodeID u; std::vector<uint8_t> bytes; size_t begin; size_t cap; };
        std::vector<PendingUpdate> pending; pending.reserve(pins_of_he.size());
        bool needs_repack = false;
        for (const HypernodeID u : pins_of_he) {
            if (!nodeIsEnabled(u)) continue;
            const auto HN = hypernode(u);
            const size_t header_begin = _hypernode_offsets[u];
            const size_t begin = HN.firstEntry();
            const size_t end = HN.firstInvalidEntry();
            // Decode existing incidents (bounded by stored degree and skipping zero-gap padding)
            const size_t expected = static_cast<size_t>(HN.degree());
            std::vector<HyperedgeID> inc; inc.reserve(expected + 1);
            size_t cur = begin; HyperedgeID acc = 0; size_t emitted = 0;
            while (cur < end && emitted < expected) {
                const uint64_t g = decode_varint_bounded(_compressed_incident_nets, cur, end);
                if (emitted > 0 && g == 0) { continue; }
                acc += static_cast<HyperedgeID>(g);
                if (acc < _num_hyperedges) inc.push_back(acc);
                ++emitted;
            }
            // Insert he and encode
            inc.push_back(he);
            std::sort(inc.begin(), inc.end()); inc.erase(std::unique(inc.begin(), inc.end()), inc.end());
            std::vector<uint8_t> bytes;
            encode_varint(static_cast<uint64_t>(inc.size()), bytes);
            HyperedgeID prev = 0; for (HyperedgeID x : inc) { ASSERT(x >= prev); encode_varint(static_cast<uint64_t>(x - prev), bytes); prev = x; }
            size_t next_begin = (u+1 < _num_hypernodes) ? _hypernode_offsets[u + 1] : _compressed_incident_nets.size();
            size_t cap = next_begin - header_begin; // bytes we can occupy without moving neighbors
            if (bytes.size() > cap) needs_repack = true;
            pending.push_back(PendingUpdate{u, std::move(bytes), header_begin, cap});
        }

        if (!needs_repack) {
            // In-place apply all updates (may extend vector for last node)
            for (auto& up : pending) {
                size_t needed_end = up.begin + up.bytes.size();
                if (needed_end > _compressed_incident_nets.size()) {
                    _compressed_incident_nets.resize(needed_end);
                }
                if (!up.bytes.empty()) {
                    std::copy(up.bytes.begin(), up.bytes.end(), _compressed_incident_nets.begin() + up.begin);
                }
            }
        } else {
            // Fallback: single repack using prepared bytes for affected nodes
            std::unordered_map<HypernodeID, std::vector<uint8_t>> overrides;
            overrides.reserve(pending.size());
            for (const auto& up : pending) { overrides.emplace(up.u, up.bytes); }
            // Snapshot current begins/ends
            std::vector<size_t> old_begin(_num_hypernodes), old_end(_num_hypernodes);
            for (HypernodeID u = 0; u < _num_hypernodes; ++u) { old_begin[u] = _hypernode_offsets[u]; old_end[u] = hypernode(u).firstInvalidEntry(); }
            // Helper to decode an untouched node from snapshot and encode again
            auto encode_from_snapshot = [&](HypernodeID u) {
                std::vector<HyperedgeID> inc; size_t cur = old_begin[u]; const size_t end = old_end[u];
                // header varint starts at old_begin[u], decode it to know how many to read
                size_t header_pos = cur;
                const size_t degree = static_cast<size_t>(decode_varint_bounded(_compressed_incident_nets, header_pos, end));
                cur = header_pos;
                HyperedgeID acc = 0; size_t emitted = 0;
                while (cur < end && emitted < degree) {
                    const uint64_t g = decode_varint_bounded(_compressed_incident_nets, cur, end);
                    if (emitted > 0 && g == 0) { continue; }
                    acc += static_cast<HyperedgeID>(g);
                    if (acc < _num_hyperedges) inc.push_back(acc);
                    ++emitted;
                }
                std::vector<uint8_t> out; if (!inc.empty()) { std::sort(inc.begin(), inc.end()); inc.erase(std::unique(inc.begin(), inc.end()), inc.end()); encode_varint(static_cast<uint64_t>(inc.size()), out); HyperedgeID prev = 0; for (HyperedgeID x : inc) { ASSERT(x >= prev); encode_varint(static_cast<uint64_t>(x - prev), out); prev = x; } }
                return out;
            };
            CompressedIncidentNets new_incident_nets; new_incident_nets.reserve(_compressed_incident_nets.size() + 16);
            size_t write_pos = 0;
            for (HypernodeID u = 0; u < _num_hypernodes; ++u) {
                if (!nodeIsEnabled(u)) continue;
                _hypernode_offsets[u] = write_pos;
                std::vector<uint8_t> bytes;
                auto it = overrides.find(u);
                if (it != overrides.end()) {
                    bytes = it->second;
                } else {
                    bytes = encode_from_snapshot(u);
                }
                new_incident_nets.insert(new_incident_nets.end(), bytes.begin(), bytes.end());
                write_pos += bytes.size();
            }
            _compressed_incident_nets.swap(new_incident_nets);
            // No explicit size metadata to update
        }

        // Re-enable the edge
        enableHyperedge(he);
    }

    parallel::scalable_vector<ParallelHyperedge> removeSinglePinAndParallelHyperedges() {
        throw UnsupportedOperationException(
            "removeSinglePinAndParallelHyperedges() is not supported in compressed hypergraph");
        return { };
    }

    void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>&) {
        throw UnsupportedOperationException(
            "restoreSinglePinAndParallelNets(hes_to_restore) is not supported in compressed hypergraph");
    }

    // ####################### Initialization / Reset Functions #######################


    // ! Reset internal community information
    void copyCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
        ASSERT(community_ids.size() == UI64(_num_hypernodes));
        doParallelForAllNodes([&](const HypernodeID& hn) {
            _community_ids[hn] = community_ids[hn];
        });
    }

    void setCommunityIDs(ds::Clustering&& communities) {
        ASSERT(communities.size() == initialNumNodes());
        _community_ids = std::move(communities);
    }

    // ! Copy compressed hypergraph in parallel
    CompressedHypergraph copy(parallel_tag_t) const;

    // ! Copy compressed hypergraph sequential
    CompressedHypergraph copy() const;

    void reset() { }

    void freeInternalData() {
        // Delete any temporary contraction buffer
        freeTmpContractionBuffer();
        // Release main containers
        _hypernode_offsets.clear();
        _hypernode_offsets.shrink_to_fit();
        _hypernode_enabled.clear();
        _hypernode_enabled.shrink_to_fit();
        _compressed_incident_nets.clear();
        _compressed_incident_nets.shrink_to_fit();
        _hyperedge_offsets.clear();
        _hyperedge_offsets.shrink_to_fit();
        _hyperedge_enabled.clear();
        _hyperedge_enabled.shrink_to_fit();
    _compressed_incidence_array.clear();
    _compressed_incidence_array.shrink_to_fit();
    _hypernode_weights.clear();
    _hyperedge_weights.clear();
        _community_ids.clear();
        _community_ids.shrink_to_fit();
        _fixed_vertices = FixedVertexSupport<CompressedHypergraph>();
        // Reset stats
        _num_hypernodes = 0;
        _num_removed_hypernodes = 0;
        _removed_degree_zero_hn_weight = 0;
        _num_hyperedges = 0;
        _num_removed_hyperedges = 0;
        _max_edge_size = 0;
        _num_pins = 0;
        _total_degree = 0;
        _total_weight = 0;
    }

    void freeTmpContractionBuffer() {
        if (_tmp_contraction_buffer) {
            delete(_tmp_contraction_buffer);
            _tmp_contraction_buffer = nullptr;
        }
    }

    void memoryConsumption(utils::MemoryTreeNode* parent) const;
    size_t memoryConsumptionKB() const;

    bool verifyIncidenceArrayAndIncidentNets() {
        throw UnsupportedOperationException(
            "verifyIncidenceArrayAndIncidentNets() not supported in compressed hypergraph");
        return false;
    }

private:
    friend class CompressedHypergraphFactory;
    template <typename Hypergraph,
          typename ConnectivityInformation>
    friend class PartitionedHypergraph;
    template<typename Hypergraph>
    friend class CommunitySupport;

    // ####################### Hypernode Information #######################

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode hypernode(const HypernodeID u) const {
        ASSERT(u < _num_hypernodes, "Hypernode" << u << " does not exist");
        return Hypernode(const_cast<CompressedHypergraph*>(this), u);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode hypernode(const HypernodeID u) {
        return Hypernode(this, u);
    }

   MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE IteratorRange<IncidentNetsIterator>
   incident_nets_of(const HypernodeID u, const size_t pos = 0) const {
       ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
       if (pos != 0) {
           throw UnsupportedOperationException(
               "incident_nets_of(u, pos) with pos > 0 is not supported in compressed hypergraph");
        }
       const size_t begin = hypernode(u).firstEntry();
       const size_t end = hypernode(u).firstInvalidEntry();
       size_t p = begin; const size_t deg = static_cast<size_t>(hypernode(u).degree());
       return IteratorRange<IncidentNetsIterator>(
           IncidentNetsIterator(_compressed_incident_nets, p, end, _num_hyperedges, &_hyperedge_enabled, deg),
           IncidentNetsIterator(_compressed_incident_nets, end));
   }

    // ####################### Hyperedge Information #######################

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge hyperedge(const HyperedgeID e) const {
        ASSERT(e < _num_hyperedges, "Hyperedge " << e << " does not exist");
        return Hyperedge(const_cast<CompressedHypergraph*>(this), e);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge hyperedge(const HyperedgeID e) {
        return Hyperedge(this, e);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID pinAt(const HyperedgeID e, const size_t local_pos) const {
        ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
        const auto he = hyperedge(e);
        ASSERT(local_pos < he.size());
        // Decode varints up to local_pos, skipping zero-gap padding bytes
        size_t cur = he.firstEntry();
        const size_t end = he.firstInvalidEntry();
        [[maybe_unused]] const size_t esize = static_cast<size_t>(he.size());
        ASSERT(local_pos < esize);
        HypernodeID acc = 0;
        size_t emitted = 0;
        while (cur < end && emitted <= local_pos) {
            const uint64_t gap = decode_varint_bounded(_compressed_incidence_array, cur, end);
            if (emitted > 0 && gap == 0) { continue; }
            acc += static_cast<HypernodeID>(gap);
            if (emitted == local_pos) return acc;
            ++emitted;
        }
        return acc;
    }



    void allocateTmpContractionBuffer() {
        if (!_tmp_contraction_buffer) {
            _tmp_contraction_buffer = new TmpContractionBuffer(
                _num_hypernodes, _num_hyperedges, _num_pins);
        }
    }

    // Member variables
    HypernodeID _num_hypernodes;
    HypernodeID _num_removed_hypernodes;
    HypernodeWeight _removed_degree_zero_hn_weight;
    HyperedgeID _num_hyperedges;
    HyperedgeID _num_removed_hyperedges;
    HypernodeID _max_edge_size;
    HypernodeID _num_pins;
    HypernodeID _total_degree;
    HypernodeWeight _total_weight;

    // Hypergraph structure with compression
    std::vector<size_t> _hypernode_offsets;   // CSR offsets into _compressed_incident_nets
    BitVector _hypernode_enabled;             // true = enabled, false = disabled
    CompressedIncidentNets _compressed_incident_nets;    // Compressed incident nets storage
    std::vector<size_t> _hyperedge_offsets;   // CSR offsets into _compressed_incidence_array
    BitVector _hyperedge_enabled;             // true = enabled, false = disabled
    CompressedIncidenceArray _compressed_incidence_array;  // Compressed pins storage
    // Lazy, two-level weight storage (default 1 when empty)
    mt_kahypar::ds::TwoLevelVector<WeightBaseT, HypernodeWeight> _hypernode_weights;
    mt_kahypar::ds::TwoLevelVector<WeightBaseT, HyperedgeWeight> _hyperedge_weights;

    // Communities and fixed vertices (uncompressed as requested)
    ds::Clustering _community_ids;
    FixedVertexSupport<CompressedHypergraph> _fixed_vertices;

    // Data reused throughout multilevel hierarchy
    TmpContractionBuffer* _tmp_contraction_buffer = nullptr;
};

} // namespace ds
} // namespace mt_kahypar