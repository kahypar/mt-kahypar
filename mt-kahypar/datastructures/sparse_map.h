/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "mt-kahypar/macros.h"
#include "kahypar/datastructure/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename Key = Mandatory,
          typename Value = Mandatory>
class CacheEfficientSparseMap {

  static constexpr size_t CACHE_EFFICIENT_MAP_SIZE = 1024;

  struct MapElement {
    Key key;
    Value value;
  };

  struct SparseElement {
    MapElement* element;
    size_t timestamp;
  };

 public:
  // TODO(heuer): Only use cache-efficient part if max_size > CACHE_EFFICIENT_MAP_SIZE
  explicit CacheEfficientSparseMap(const Key max_size, const Value initial_value) :
    _max_size(max_size),
    _initial_value(initial_value),
    _data(std::make_unique<uint8_t[]>(
      CACHE_EFFICIENT_MAP_SIZE * sizeof(MapElement) +
      CACHE_EFFICIENT_MAP_SIZE * sizeof(SparseElement))),
    _size(0),
    _timestamp(1),
    _sparse(nullptr),
    _dense(nullptr) {

    _sparse = reinterpret_cast<SparseElement*>(_data.get());
    _dense = reinterpret_cast<MapElement*>(
      _data.get() + sizeof(SparseElement) * CACHE_EFFICIENT_MAP_SIZE);
    for (size_t i = 0; i < CACHE_EFFICIENT_MAP_SIZE; ++i) {
      _dense[i] = MapElement { std::numeric_limits<Key>::max(), initial_value };
      _sparse[i] = SparseElement { &_dense[i], 0 };
    }
  }

  CacheEfficientSparseMap(const CacheEfficientSparseMap&) = delete;
  CacheEfficientSparseMap& operator= (const CacheEfficientSparseMap& other) = delete;

  ~CacheEfficientSparseMap() = default;

  size_t size() const {
    return _size;
  }

  bool contains(const Key key) const {
    MapElement* element = find(key);
    return element != nullptr;
  }

  void add(const Key key, const Value value) {
    if ( _size < CACHE_EFFICIENT_MAP_SIZE - 2 ) {
      addElement(key, value);
    }
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

  Value& operator[] (const Key key) {
    MapElement* element = find(key);
    if ( element != nullptr ) {
      return element->value;
    } else if ( _size < CACHE_EFFICIENT_MAP_SIZE - 2 ) {
      return addElement(key, _initial_value)->value;
    }
    return _dense[0].value;
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key));
    return find(key)->value;
  }

  void clear() {
    _size = 0;
    ++_timestamp;
  }

 private:
  inline MapElement* find(const Key key) const {
    size_t hash = key % CACHE_EFFICIENT_MAP_SIZE;
    while ( _sparse[hash].timestamp == _timestamp ) {
      if ( _sparse[hash].element->key == key ) {
        return _sparse[hash].element;
      }
      hash = (hash + 1) % CACHE_EFFICIENT_MAP_SIZE;
    }
    return nullptr;
  }

  inline MapElement* addElement(const Key key, const Value value) {
    ASSERT(_size < CACHE_EFFICIENT_MAP_SIZE - 2);
    size_t hash = key % CACHE_EFFICIENT_MAP_SIZE;
    while ( _sparse[hash].timestamp == _timestamp ) {
      hash = (hash + 1) % CACHE_EFFICIENT_MAP_SIZE;
    }

    _dense[_size] = MapElement { key, value };
    _sparse[hash] = SparseElement { &_dense[_size++], _timestamp };
    return _sparse[hash].element;
  }

  const size_t _max_size;
  const Value _initial_value;
  std::unique_ptr<uint8_t[]> _data;

  size_t _size;
  size_t _timestamp;
  SparseElement* _sparse;
  MapElement* _dense;
};

} // namespace ds
} // namespace mt_kahypar