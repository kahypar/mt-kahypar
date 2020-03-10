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

#include <cmath>

#include "mt-kahypar/macros.h"
#include "kahypar/datastructure/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename Key = Mandatory,
          typename Value = Mandatory>
class CacheEfficientSparseMap {

  static constexpr size_t CACHE_EFFICIENT_MAP_SIZE = 8192;
  static constexpr size_t REHASH_THRESHOLD = 4096;

  static_assert(REHASH_THRESHOLD < CACHE_EFFICIENT_MAP_SIZE,
    "Rehash threshold must be smaller than map size");

  struct MapElement {
    Key key;
    Value value;
  };

  struct SparseElement {
    MapElement* element;
    size_t timestamp;
  };

 public:
  explicit CacheEfficientSparseMap(const Key max_size, const Value initial_value) :
    _max_size(std::pow(2.0, std::ceil(std::log2(static_cast<double>(max_size))))),
    _initial_value(initial_value),
    _data(std::make_unique<uint8_t[]>(
      _max_size * sizeof(MapElement) +
      _max_size * sizeof(SparseElement))),
    _size(0),
    _capacity(std::min(static_cast<size_t>(_max_size), CACHE_EFFICIENT_MAP_SIZE)),
    _rehash_threshold(static_cast<size_t>(_max_size) < CACHE_EFFICIENT_MAP_SIZE ? _max_size : REHASH_THRESHOLD),
    _timestamp(1),
    _sparse(nullptr),
    _dense(nullptr),
    _rehash(0),
    _ops(0) {
    _dense = reinterpret_cast<MapElement*>(_data.get());
    _sparse = reinterpret_cast<SparseElement*>(
      _data.get() + sizeof(MapElement) * _capacity);
    memset(_data.get(), 0, _max_size * (sizeof(MapElement) + sizeof(SparseElement)));
  }

  CacheEfficientSparseMap(const CacheEfficientSparseMap&) = delete;
  CacheEfficientSparseMap& operator= (const CacheEfficientSparseMap& other) = delete;

  ~CacheEfficientSparseMap() {
    LOG << V(_rehash) << V(_ops);
  }

  size_t size() const {
    return _size;
  }

  bool contains(const Key key) const {
    SparseElement* s = find(key);
    return containsValidElement(key, s);
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
    SparseElement* s = find(key);
    ++_ops;
    if ( containsValidElement(key, s) ) {
      ASSERT(s->element);
      return s->element->value;
    } else if ( _size < _rehash_threshold ) {
      return addElement(key, _initial_value, s)->value;
    } else {
      rehash();
      ASSERT(key < _max_size);
      s = &_sparse[key];
      ASSERT(s->timestamp < _timestamp);
      return addElement(key, _initial_value, s)->value;
    }
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key));
    return find(key)->element->value;
  }

  void clear() {
    _size = 0;
    _capacity = std::min(static_cast<size_t>(_max_size), CACHE_EFFICIENT_MAP_SIZE);
    if ( _max_size >= CACHE_EFFICIENT_MAP_SIZE && _rehash_threshold == _max_size ) {
      _dense = reinterpret_cast<MapElement*>(_data.get());
      _sparse = reinterpret_cast<SparseElement*>(
        _data.get() + sizeof(MapElement) * _capacity);
      memset(_sparse, 0, _capacity * sizeof(SparseElement));
    }
    _rehash_threshold = _max_size < CACHE_EFFICIENT_MAP_SIZE ? _max_size : REHASH_THRESHOLD;
    ++_timestamp;
  }

 private:
  inline SparseElement* find(const Key key) const {
    size_t hash = key & ( _capacity - 1 );
    while ( _sparse[hash].timestamp == _timestamp ) {
      ASSERT(_sparse[hash].element);
      if ( _sparse[hash].element->key == key ) {
        return &_sparse[hash];
      }
      hash = (hash + 1) & ( _capacity - 1 );
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

  void rehash() {
    ++_rehash;
    ASSERT(_capacity < _max_size);
    _capacity = _max_size;
    _rehash_threshold = _max_size;
    _dense = reinterpret_cast<MapElement*>(_data.get());
    _sparse = reinterpret_cast<SparseElement*>(
      _data.get() + sizeof(MapElement) * _capacity);

    for ( MapElement& element : *this ) {
      const Key& key = element.key;
      ASSERT(key < _max_size);
      ASSERT(_sparse[key].timestamp < _timestamp);
      _sparse[key] = SparseElement { &element, _timestamp };
    }
  }

  const size_t _max_size;
  const Value _initial_value;
  std::unique_ptr<uint8_t[]> _data;

  size_t _size;
  size_t _capacity;
  size_t _rehash_threshold;
  size_t _timestamp;
  SparseElement* _sparse;
  MapElement* _dense;

  size_t _rehash;
  size_t _ops;
};

} // namespace ds
} // namespace mt_kahypar