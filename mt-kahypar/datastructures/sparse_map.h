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
/*
 * Sparse map based on sparse set representation of
 * Briggs, Preston, and Linda Torczon. "An efficient representation for sparse sets."
 * ACM Letters on Programming Languages and Systems (LOPLAS) 2.1-4 (1993): 59-69.
 */

#pragma once

#include <cmath>

#include "mt-kahypar/macros.h"
#include "kahypar/datastructure/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

/*
 * Sparse map based on sparse set representation of
 * Briggs, Preston, and Linda Torczon. "An efficient representation for sparse sets."
 * ACM Letters on Programming Languages and Systems (LOPLAS) 2.1-4 (1993): 59-69.
 */

/*!
 * Sparse map implementation that uses a fixed size.
 * In contrast to the implementation in KaHyPar (see kahypar/datastructure/sparse_map.h),
 * which uses as size the cardinality of the key universe, hash collisions have to be handled
 * explicitly. Hash collisions are resolved with linear probing.
 * Advantage of the implementation is that it uses significantly less space than the
 * version in KaHyPar and should be therefore more cache-efficient.
 * Note, there is no fallback strategy if all slots of the sparse map are occupied by an
 * element. Please make sure that no more than MAP_SIZE elements are inserted into the
 * sparse map. Otherwise, the behavior is undefined.
 */
template <typename Key = Mandatory,
          typename Value = Mandatory>
class FixedSizeSparseMap {

  struct MapElement {
    Key key;
    Value value;
  };

  struct SparseElement {
    MapElement* element;
    size_t timestamp;
  };

 public:
  static constexpr size_t MAP_SIZE = 16384;

  static_assert(is_power_of_two(MAP_SIZE), "Size of map is not a power of two!");

  explicit FixedSizeSparseMap(const Value initial_value) :
    _initial_value(initial_value),
    _data(std::make_unique<uint8_t[]>(
      MAP_SIZE * sizeof(MapElement) + MAP_SIZE * sizeof(SparseElement))),
    _size(0),
    _timestamp(1),
    _sparse(reinterpret_cast<SparseElement*>(_data.get())),
    _dense(reinterpret_cast<MapElement*>(_data.get() +  + sizeof(SparseElement) * MAP_SIZE)) {
    memset(_data.get(), 0, MAP_SIZE * (sizeof(MapElement) + sizeof(SparseElement)));
  }

  FixedSizeSparseMap(const FixedSizeSparseMap&) = delete;
  FixedSizeSparseMap& operator= (const FixedSizeSparseMap& other) = delete;

  ~FixedSizeSparseMap() = default;

  size_t size() const {
    return _size;
  }

  const MapElement* begin() const {
    return _dense;
  }

  const MapElement* end() const {
    return _dense + _size;
  }

  MapElement* begin() {
    return _dense;
  }

  MapElement* end() {
    return _dense + _size;
  }

  bool contains(const Key key) const {
    SparseElement* s = find(key);
    return containsValidElement(key, s);
  }

  Value& operator[] (const Key key) {
    SparseElement* s = find(key);
    if ( containsValidElement(key, s) ) {
      ASSERT(s->element);
      return s->element->value;
    } else {
      return addElement(key, _initial_value, s)->value;
    }
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key));
    return find(key)->element->value;
  }

  void clear() {
    _size = 0;
    ++_timestamp;
  }

 private:
  inline SparseElement* find(const Key key) const {
    size_t hash = key & ( MAP_SIZE - 1 );
    while ( _sparse[hash].timestamp == _timestamp ) {
      ASSERT(_sparse[hash].element);
      if ( _sparse[hash].element->key == key ) {
        return &_sparse[hash];
      }
      hash = (hash + 1) & ( MAP_SIZE - 1 );
    }
    return &_sparse[hash];
  }

  inline bool containsValidElement(const Key key,
                                   const SparseElement* s) const {
    unused(key);
    ASSERT(s);
    const bool is_contained = s->timestamp == _timestamp;
    ASSERT(!is_contained || s->element->key == key);
    return is_contained;
  }

  inline MapElement* addElement(const Key key,
                                const Value value,
                                SparseElement* s) {
    ASSERT(find(key) == s);
    _dense[_size] = MapElement { key, value };
    *s = SparseElement { &_dense[_size++], _timestamp };
    return s->element;
  }

  const Value _initial_value;
  std::unique_ptr<uint8_t[]> _data;

  size_t _size;
  size_t _timestamp;
  SparseElement* _sparse;
  MapElement* _dense;
};

} // namespace ds
} // namespace mt_kahypar