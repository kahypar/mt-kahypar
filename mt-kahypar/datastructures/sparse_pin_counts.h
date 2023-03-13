/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <tbb/concurrent_vector.h>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"


namespace mt_kahypar {
namespace ds {

/**
 * This is a sparse implementation of the pin count data structure.
 * Our original data structure for the pin count values takes O(k*|E|) space which would
 * be not feasible when k is large. The sparse implementation uses the observation
 * that the connectivity for most hyperedges is small in real-world hypergraphs.
 * The data structure stores for each hyperedge at most c <= k tuples of the form (block, pin_count),
 * where pin_count is the number of nodes part of the given block. If the connectivity of a
 * hyperedge is larger then c, then the data structure explicitly stores all pin count values
 * in an external pin count list. However, this external list is only initialized when the connectivity becomes larger then
 * c, which should happen rarely in practice. Thus, the data structures takes O(c * |E|) space when
 * overflows are negligible.
 *
 * The data structure supports concurrent read, but only one thread can modify the pin count values of
 * a hyperedge. Multiple writes to different hyperedges are supported.
 */
class SparsePinCounts {

  static constexpr bool debug = false;

  static constexpr size_t MAX_ENTRIES_PER_HYPEREDGE = 8; // = c

  struct PinCountHeader {
    // Stores the connectivity of a hyperedge
    PartitionID connectivity;
    // Flag that indicates whether or not pin counts are stored
    // in the external pin count list
    bool is_external;
  };

  // Stores the number of pins contained in a block
  struct PinCountEntry {
    PartitionID block;
    HypernodeID pin_count;
  };

 public:
  using Value = char;

  SparsePinCounts() :
    _num_hyperedges(0),
    _k(0),
    _entries_per_hyperedge(0),
    _size_of_pin_counts_per_he(0),
    _pin_count_in_part(),
    _pin_count_ptr(nullptr),
    _ext_pin_count_list(),
    _num_overflows(0) { }

  SparsePinCounts(const HyperedgeID num_hyperedges,
                  const PartitionID k,
                  const HypernodeID max_value,
                  const bool assign_parallel = true) :
    _num_hyperedges(0),
    _k(0),
    _entries_per_hyperedge(0),
    _size_of_pin_counts_per_he(0),
    _pin_count_in_part(),
    _pin_count_ptr(nullptr),
    _ext_pin_count_list(),
    _num_overflows(0) {
    initialize(num_hyperedges, k, max_value, assign_parallel);
  }

  SparsePinCounts(const SparsePinCounts&) = delete;
  SparsePinCounts & operator= (const SparsePinCounts &) = delete;

  SparsePinCounts(SparsePinCounts&& other) :
    _num_hyperedges(other._num_hyperedges),
    _k(other._k),
    _entries_per_hyperedge(other._entries_per_hyperedge),
    _pin_count_in_part(std::move(other._pin_count_in_part)),
    _pin_count_ptr(nullptr),
    _ext_pin_count_list(std::move(other._ext_pin_count_list)),
    _num_overflows(std::move(other._num_overflows)) {
    _pin_count_ptr = _pin_count_in_part.data();
  }

  SparsePinCounts & operator= (SparsePinCounts&& other) {
    _num_hyperedges = other._num_hyperedges;
    _k = other._k;
    _entries_per_hyperedge = other._entries_per_hyperedge;
    _size_of_pin_counts_per_he = other._size_of_pin_counts_per_he;
    _pin_count_in_part = std::move(other._pin_count_in_part);
    _pin_count_ptr = _pin_count_in_part.data();
    _ext_pin_count_list = std::move(other._ext_pin_count_list);
    _num_overflows = std::move(other._num_overflows);
    return *this;
  }

  // ! Initializes the data structure
  void initialize(const HyperedgeID num_hyperedges,
                  const PartitionID k,
                  const HypernodeID,
                  const bool assign_parallel = true) {
    _num_hyperedges = num_hyperedges;
    _k = k;
    _entries_per_hyperedge = std::min(
      static_cast<size_t>(k), MAX_ENTRIES_PER_HYPEREDGE);
    _size_of_pin_counts_per_he = sizeof(PinCountHeader) +
      sizeof(PinCountEntry) * _entries_per_hyperedge;
    _pin_count_in_part.resize("Refinement", "pin_count_in_part",
      _size_of_pin_counts_per_he * num_hyperedges, false, assign_parallel);
    _pin_count_ptr = _pin_count_in_part.data();
    _ext_pin_count_list.resize(_num_hyperedges);
    reset(assign_parallel);
  }

  void reset(const bool assign_parallel = true) {
    if ( assign_parallel ) {
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID he) {
        init_pin_count_of_hyperedge(he);
        _ext_pin_count_list[he].clear();
      });
    } else {
      for ( HyperedgeID he = 0; he < _num_hyperedges; ++he ) {
        init_pin_count_of_hyperedge(he);
        _ext_pin_count_list[he].clear();
      }
    }
  }

  // ! Returns the pin count of the hyperedge in the corresponding block
  inline HypernodeID pinCountInPart(const HyperedgeID he,
                                    const PartitionID id) const {
    const PinCountEntry* val = find_entry(he, id);
    return val ? val->pin_count : 0;
  }

  // ! Sets the pin count of the hyperedge in the corresponding block to value
  inline void setPinCountInPart(const HyperedgeID he,
                                const PartitionID id,
                                const HypernodeID value) {
    add_pin_count_entry(he, id, value);
  }

  // ! Increments the pin count of the hyperedge in the corresponding block
  inline HypernodeID incrementPinCountInPart(const HyperedgeID he,
                                             const PartitionID id) {
    PinCountEntry* val = find_entry(he, id);
    HypernodeID inc_pin_count = 0;
    if ( val ) {
      inc_pin_count = ++val->pin_count;
    } else {
      inc_pin_count = 1;
      add_pin_count_entry(he, id, inc_pin_count);
    }
    return inc_pin_count;
  }

  // ! Decrements the pin count of the hyperedge in the corresponding block
  inline HypernodeID decrementPinCountInPart(const HyperedgeID he,
                                             const PartitionID id) {
    PinCountEntry* val = find_entry(he, id);
    ASSERT(val);
    const HypernodeID dec_pin_count = --val->pin_count;
    if ( dec_pin_count == 0 ) {
      // Remove pin count entry
      // Note that only one thread can modify the pin count list of
      // a hyperedge at the same time. Therefore, this operation is thread-safe.
      PinCountHeader* head = header(he);
      --head->connectivity;
      if ( likely(!head->is_external) ) {
        PinCountEntry* back = entry(he, head->connectivity);
        *val = *back;
        back->block = kInvalidPartition;
        back->pin_count = 0;
      } else {
        // Note that in case the connectivity becomes smaller than c,
        // we do not fallback to the smaller pin count list bounded by c.
        size_t pos = 0;
        for ( ; pos < _ext_pin_count_list[he].size(); ++pos ) {
          if ( _ext_pin_count_list[he][pos].block == id ) {
            break;
          }
        }
        std::swap(_ext_pin_count_list[he][pos], _ext_pin_count_list[he][head->connectivity]);
        _ext_pin_count_list[he][head->connectivity].block = kInvalidPartition;
        _ext_pin_count_list[he][head->connectivity].pin_count = 0;
      }
    }
    return dec_pin_count;
  }

  // ! Returns the size in bytes of this data structure
  size_t size_in_bytes() const {
    // TODO: size of external list is missing
    return sizeof(char) * _pin_count_in_part.size();
  }

  static size_t num_elements(const HyperedgeID num_hyperedges,
                             const PartitionID k,
                             const HypernodeID) {
    const size_t entries_per_hyperedge = std::min(
      static_cast<size_t>(k), MAX_ENTRIES_PER_HYPEREDGE);
    const size_t size_of_pin_counts_per_he = sizeof(PinCountHeader) +
      sizeof(PinCountEntry) * entries_per_hyperedge;
    return size_of_pin_counts_per_he * num_hyperedges;
  }

 private:
  inline void init_pin_count_of_hyperedge(const HyperedgeID& he) {
    PinCountHeader* head = header(he);
    head->connectivity = 0;
    head->is_external = false;
    for ( size_t i = 0; i < _entries_per_hyperedge; ++i ) {
      PinCountEntry* pin_count = entry(he, i);
      pin_count->block = kInvalidPartition;
      pin_count->pin_count = 0;
    }
  }

  inline void add_pin_count_entry(const HyperedgeID he,
                                  const PartitionID id,
                                  const HypernodeID value) {
    // Assumes that the block with the given ID does not exist
    // and inserts it at the end of the pin count list
    // Note that only one thread can modify the pin count list of
    // a hyperedge at the same time. Therefore, this operation is thread-safe.
    PinCountHeader* head = header(he);
    if ( likely(!head->is_external) ) {
      const size_t connectivity = head->connectivity;
      if ( connectivity < _entries_per_hyperedge ) {
        // Still enough entries to add the pin count entry
        PinCountEntry* pin_count = entry(he, connectivity);
        pin_count->block = id;
        pin_count->pin_count = value;
      } else {
        // Connecitivity is now larger than c
        // => copy entries to external pin count list
        handle_overflow(he);
        add_pin_count_entry_to_external(he, id, value);
      }
    } else {
      add_pin_count_entry_to_external(he, id, value);
    }
    ++head->connectivity;
  }

  inline void handle_overflow(const HyperedgeID& he) {
    PinCountHeader* head = header(he);
    // Copy entries to external pin count list
    for ( size_t i = 0; i < _entries_per_hyperedge; ++i ) {
      _ext_pin_count_list[he].push_back(*entry(he, i));
    }
    ++_num_overflows;
    head->is_external = true;
  }

  inline void add_pin_count_entry_to_external(const HyperedgeID he,
                                              const PartitionID id,
                                              const HypernodeID value) {
    PinCountHeader* head = header(he);
    ASSERT(head->is_external);
    if ( static_cast<size_t>(head->connectivity) < _ext_pin_count_list[he].size() ) {
      // Reuse existing entry that was removed due to decrementing the pin count
      ASSERT(_ext_pin_count_list[he][head->connectivity].block == kInvalidPartition);
      _ext_pin_count_list[he][head->connectivity].block = id;
      _ext_pin_count_list[he][head->connectivity].pin_count = value;
    } else {
      _ext_pin_count_list[he].push_back(PinCountEntry { id, value });
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const PinCountEntry* find_entry(const HyperedgeID he, const PartitionID id) const {
    const PinCountHeader* head = header(he);
    if ( likely(!head->is_external) ) {
      // Due to concurrent writes, the connectivity can become larger than MAX_ENTRIES_PER_HYPEREDGE.
      const size_t connectivity =
        std::min(static_cast<size_t>(head->connectivity), MAX_ENTRIES_PER_HYPEREDGE);
      for ( size_t i = 0; i < connectivity; ++i ) {
        const PinCountEntry* value = entry(he, i);
        if ( value->block == id ) {
          return value;
        }
      }
    } else {
      const size_t num_entries = head->connectivity;
      for ( size_t i = 0; i < num_entries; ++i ) {
        const PinCountEntry& value = _ext_pin_count_list[he][i];
        if ( value.block == id ) {
          return &value;
        }
      }
    }
    return nullptr;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PinCountEntry* find_entry(const HyperedgeID he, const PartitionID id) {
    return const_cast<PinCountEntry*>(static_cast<const SparsePinCounts&>(*this).find_entry(he, id));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const PinCountHeader* header(const HyperedgeID he) const {
    ASSERT(he <= _num_hyperedges, "Hyperedge" << he << "does not exist");
    return reinterpret_cast<const PinCountHeader*>(_pin_count_ptr + he * _size_of_pin_counts_per_he);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PinCountHeader* header(const HyperedgeID he) {
    return const_cast<PinCountHeader*>(static_cast<const SparsePinCounts&>(*this).header(he));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const PinCountEntry* entry(const HyperedgeID he,
                                                                const size_t idx) const {
    ASSERT(he <= _num_hyperedges, "Hyperedge" << he << "does not exist");
    return reinterpret_cast<const PinCountEntry*>(_pin_count_ptr +
      he * _size_of_pin_counts_per_he + sizeof(PinCountHeader) + sizeof(PinCountEntry) * idx);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PinCountEntry* entry(const HyperedgeID he,
                                                          const size_t idx) {
    return const_cast<PinCountEntry*>(static_cast<const SparsePinCounts&>(*this).entry(he, idx));
  }

  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;

  // ! Number of blocks
  PartitionID _k;

  // ! Maximum number of pin count entries per hyperedge (= c)
  size_t _entries_per_hyperedge;

  // ! Size in bytes of the header struct and all pin count entries
  size_t _size_of_pin_counts_per_he;

  // ! Stores the pin count list bounded by c
  Array<char> _pin_count_in_part;
  char* _pin_count_ptr;

  // ! External pin count list that stores the pin count values when
  // ! the connectivity becomes larger than c.
  // ! Note that we have to use concurrent_vector since we allow concurrent
  // ! read while modyfing the entries.
  vec<tbb::concurrent_vector<PinCountEntry>> _ext_pin_count_list;

  parallel::AtomicWrapper<size_t> _num_overflows;
};
}  // namespace ds
}  // namespace mt_kahypar
