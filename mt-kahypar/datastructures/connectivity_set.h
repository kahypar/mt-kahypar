/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <atomic>
#include <type_traits>

#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/parallel/copyable_atomic.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/bit_magic.h"

namespace mt_kahypar {
namespace ds {

/**
 * Can store a set of partition ids (< k = number of blocks).
 *
 * Each connectivity set is associated with exactly one hyperedge. For each hyperedge it stores
 * the block ids which its corresponding pins belongs to. If a pin of the hyperedge is
 * moved to a different block the connectivity can increase or decrease.
 * This class provides a lock-free and thread-safe way to handle connectivity changes in a multi-
 * threaded setting. The block ids are stored in a bitset. If the connectivity increased or decreased
 * one can add or remove the corresponding block id from the set with add(...) or remove(...).
 * Also, the class provides a thread-safe iterator over a snapshot of the connectivity set.
 *
 * IMPORTANT
 * ---------
 * In some situations the number of set ones in the bitset can differ from value returned by size().
 * Think of the situation, where a vertex move decreased the connectivity of block x. In parallel,
 * an other vertex moves to block x and increased the connectivity again. Since, connectivity of block
 * x first decreased and than increased, the corresponding bit for block x is set before updating the
 * connectivity set. However, since both moves happens simultanously it can happen that the connectivity
 * set is first updated with the increase operation and afterwards with the decrease operation.
 * To handle that case, we perform a XOR operation on the bit position for block x ( bit for block x XOR 1 ).
 * In the concrete example this means, that the increase operation for block x would set the corresponding
 * bit to zero, but connectivity is temporary inceased. However, in such situations there must be always an
 * decrease operation running in parallel, which will decrease the connectivity again and set the bit
 * for block x to 1 again.
 * In fact of that, the connectivity set reflects the real connectivity of a hyperedge if all parallel
 * operations are finished. In a parallel setting, this set can be used as indicator for the real connectivity
 * set with the assumptions that conflicts happens very rarely.
 */
class ConnectivitySet {

  static constexpr bool debug = false;

  using PartitionID = int32_t;
  using ConnectivityAtomic = parallel::CopyableAtomic<PartitionID>;
  using BitsetAtomic = parallel::CopyableAtomic<uint8_t>;


 public:
  /**
   * Thread-safe iterator over the set of block ids of the connectivity set.
   * To do so, the iterator makes a copy of the 8-bit word (which has at least one set bit)
   * over which it will iterate next. During that time changes on that part of the bitset will
   * not be reflected in the iterator. However, any other modification on the bitset will
   * be also immediately visible in the iterator.
   */
  class ConnectivitySetIterator :
    public std::iterator<std::forward_iterator_tag,  // iterator_category
                         PartitionID,                // value_type
                         std::ptrdiff_t,             // difference_type
                         const PartitionID*,         // pointer
                         PartitionID> {              // reference

   public:

    ConnectivitySetIterator() = default;

    ConnectivitySetIterator(const ConnectivitySetIterator& other) = default;
    ConnectivitySetIterator& operator= (const ConnectivitySetIterator& other) = default;

    ConnectivitySetIterator(ConnectivitySetIterator&& other) = default;
    ConnectivitySetIterator& operator= (ConnectivitySetIterator&& other) = default;

    ~ConnectivitySetIterator() = default;


    ConnectivitySetIterator(const parallel::scalable_vector<BitsetAtomic>& bitset) :
      _bitset(bitset),
      _current_id(-1),
      _current_bitset_id(0),
      _current_bitset(_bitset[0]) {
      next();
    }

    ConnectivitySetIterator(const parallel::scalable_vector<BitsetAtomic>& bitset,
                            const PartitionID current_bitset_id,
                            const uint8_t current_bitset) :
      _bitset(bitset),
      _current_id(-1),
      _current_bitset_id(current_bitset_id),
      _current_bitset(current_bitset) { }

    // ! Returns the id of the element the iterator currently points to.
    PartitionID operator* () const {
      ASSERT(_current_id < utils::count[_current_bitset]);
      return 8 * _current_bitset_id + utils::select[_current_bitset][_current_id + 1];
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    ConnectivitySetIterator& operator++ () {
      next();
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    ConnectivitySetIterator operator++ (int) {
      ConnectivitySetIterator copy = *this;
      operator++ ();
      return copy;
    }

    // ! Convenience function for range-based for-loops
    friend ConnectivitySetIterator en(const std::pair<ConnectivitySetIterator,
                                                      ConnectivitySetIterator>& iter_pair);
    // ! Convenience function for range-based for-loops
    friend ConnectivitySetIterator begin(const std::pair<ConnectivitySetIterator,
                                                         ConnectivitySetIterator>& iter_pair);

    bool operator!= (const ConnectivitySetIterator& rhs) {
      return _current_bitset_id != rhs._current_bitset_id ||
             _current_bitset != rhs._current_bitset;
    }

    bool operator== (const ConnectivitySetIterator& rhs) {
      return _current_bitset_id == rhs._current_bitset_id &&
             _current_bitset == rhs._current_bitset;
    }

   private:
    void next() {
      ++_current_id;
      while ( _current_id + 1 > utils::count[_current_bitset] ) {
        _current_id = 0;
        ++_current_bitset_id;
        if ( _current_bitset_id < (PartitionID) _bitset.size() ) {
          _current_bitset = _bitset[_current_bitset_id];
        } else {
          _current_bitset = 0;
          _current_bitset_id = _bitset.size();
          break;
        }
      }
    }

    const parallel::scalable_vector<BitsetAtomic>& _bitset;
    PartitionID _current_id;
    PartitionID _current_bitset_id;
    uint8_t _current_bitset;
  };

  ConnectivitySet(const PartitionID k) :
    _connectivity(0),
    _bitset() {
    size_t num_entries = k / 8 + ( k % 8 > 0 ? 1 : 0 );
    _bitset.assign(num_entries, BitsetAtomic(0));
  }

  ConnectivitySet(const ConnectivitySet&) = delete;
  ConnectivitySet& operator= (const ConnectivitySet&) = delete;

  ConnectivitySet(ConnectivitySet&& other) = default;
  ConnectivitySet& operator= (ConnectivitySet&&) = default;

  ~ConnectivitySet() = default;

  std::pair<ConnectivitySetIterator, ConnectivitySetIterator> connectivitySet() const {
    return std::make_pair( ConnectivitySetIterator(_bitset),
                           ConnectivitySetIterator(_bitset, _bitset.size(), 0) );
  }

  bool contains(const PartitionID id) const {
    ASSERT(id >= 0 && id <  ( (PartitionID) ( 8 * _bitset.size() ) ) );
    size_t entry = id / 8;
    size_t offset = ( id % 8 ) + 1;
    return _bitset[entry] & utils::bitmask[offset];
  }

  void add(const PartitionID id) {
    ASSERT(id >= 0 && id <  ( (PartitionID) ( 8 * _bitset.size() ) ) );
    ++_connectivity;
    XOR(id);
  }

  void remove(const PartitionID id) {
    ASSERT(id >= 0 && id <  ( (PartitionID) ( 8 * _bitset.size() ) ) );
    --_connectivity;
    XOR(id);
  }

  /*!
   * Clears the connectivity set
   * Note, this function is not thread-safe and should be not called
   * in a multi-threaded setting.
   */
  void clear() {
    _connectivity = 0;
    for ( size_t i = 0; i < _bitset.size(); ++i ) {
      _bitset[i] &= utils::bitmask[0];
    }
  }

  PartitionID size() const {
    return _connectivity;
  }

 private:

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void XOR(const PartitionID id) {
    ASSERT(id >= 0 && id <  ( (PartitionID) ( 8 * _bitset.size() ) ) );
    size_t entry = id / 8;
    size_t offset = ( id % 8 ) + 1;
    _bitset[entry] ^= utils::bitmask[offset];
  }

  ConnectivityAtomic _connectivity;
  parallel::scalable_vector<BitsetAtomic> _bitset;
};

} // namespace ds
} // namespace mt_kahypar