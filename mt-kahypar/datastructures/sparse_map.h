/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <cmath>

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"

namespace mt_kahypar {
namespace ds {

/*
 * Sparse map based on sparse set representation of
 * Briggs, Preston, and Linda Torczon. "An efficient representation for sparse sets."
 * ACM Letters on Programming Languages and Systems (LOPLAS) 2.1-4 (1993): 59-69.
 */

template <typename Key = Mandatory,
          typename Value = Mandatory,
          typename Derived = Mandatory>
class SparseMapBase {
 protected:
  struct MapElement {
    Key key;
    Value value;
  };

 public:
  SparseMapBase(const SparseMapBase&) = delete;
  SparseMapBase& operator= (const SparseMapBase&) = delete;

  SparseMapBase& operator= (SparseMapBase&&) = delete;

  size_t size() const {
    return _size;
  }

  void setMaxSize(const size_t max_size) {
    ASSERT(_sparse);
    _dense = reinterpret_cast<MapElement*>(_sparse + max_size);
  }

  bool contains(const Key key) const {
    return static_cast<const Derived*>(this)->containsImpl(key);
  }

  void add(const Key key, const Value value) {
    static_cast<Derived*>(this)->addImpl(key, value);
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

  void clear() {
    static_cast<Derived*>(this)->clearImpl();
  }

  Value& operator[] (const Key key) {
    const size_t index = _sparse[key];
    if (!contains(key)) {
      _dense[_size] = MapElement { key, Value() };
      _sparse[key] = _size++;
      return _dense[_size - 1].value;
    }
    return _dense[index].value;
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key), V(key));
    return _dense[_sparse[key]].value;
  }

  void freeInternalData() {
    _size = 0;
    _data = nullptr;
    _sparse = nullptr;
    _dense = nullptr;
  }

 protected:
  explicit SparseMapBase(const size_t max_size) :
    _size(0),
    _data(nullptr),
    _sparse(nullptr),
    _dense(nullptr) {
    allocate_data(max_size);
  }

  ~SparseMapBase() = default;

  SparseMapBase(SparseMapBase&& other) :
    _size(other._size),
    _data(std::move(other._data)),
    _sparse(std::move(other._sparse)),
    _dense(std::move(other._dense)) {
    other._size = 0;
    other._data = nullptr;
    other._sparse = nullptr;
    other._dense = nullptr;
  }

  void allocate_data(const size_t max_size) {
    ASSERT(!_data && !_sparse);
    const size_t num_elements = (max_size * sizeof(MapElement) + max_size * sizeof(size_t)) / sizeof(size_t);
    char* data = parallel::MemoryPool::instance().request_unused_mem_chunk(num_elements, sizeof(size_t));
    if ( data ) {
      _sparse = reinterpret_cast<size_t*>(data);
    } else {
      _data = parallel::make_unique<size_t>(num_elements);
      _sparse = reinterpret_cast<size_t*>(_data.get());
    }
    _dense = reinterpret_cast<MapElement*>(_sparse + max_size);
  }

  size_t _size;
  parallel::tbb_unique_ptr<size_t> _data;
  size_t* _sparse;
  MapElement* _dense;
};


template <typename Key = Mandatory,
          typename Value = Mandatory>
class SparseMap final : public SparseMapBase<Key, Value, SparseMap<Key, Value> >{
  using Base = SparseMapBase<Key, Value, SparseMap<Key, Value> >;
  friend Base;

 public:
  explicit SparseMap(const Key max_size) :
    Base(max_size) { }

  SparseMap(const SparseMap&) = delete;
  SparseMap& operator= (const SparseMap& other) = delete;

  SparseMap(SparseMap&& other) :
    Base(std::move(other)) { }

  SparseMap& operator= (SparseMap&& other) {
    _data = std::move(other._data);
    _sparse = std::move(other._sparse);
    _size = 0;
    _dense = std::move(other._dense);
    other._size = 0;
    other._data = nullptr;
    other._sparse = nullptr;
    other._dense = nullptr;
    return *this;
  }

  ~SparseMap() = default;

  void remove(const Key key) {
    const size_t index = _sparse[key];
    if (index < _size && _dense[index].key == key) {
      std::swap(_dense[index], _dense[_size - 1]);
      _sparse[_dense[index].key] = index;
      --_size;
    }
  }

 private:
  bool containsImpl(const Key key) const {
    const size_t index = _sparse[key];
    return index < _size && _dense[index].key == key;
  }


  void addImpl(const Key key, const Value value) {
    const size_t index = _sparse[key];
    if (index >= _size || _dense[index].key != key) {
      _dense[_size] = { key, value };
      _sparse[key] = _size++;
    }
  }

  void clearImpl() {
    _size = 0;
  }

  using Base::_data;
  using Base::_sparse;
  using Base::_dense;
  using Base::_size;
};

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

  static constexpr size_t MAP_SIZE = 32768; // Size of sparse map is approx. 1 MB

  static_assert(MAP_SIZE && ((MAP_SIZE & (MAP_SIZE - 1)) == 0UL), "Size of map is not a power of two!");

  explicit FixedSizeSparseMap(const Value initial_value) :
    _map_size(0),
    _initial_value(initial_value),
    _data(nullptr),
    _size(0),
    _timestamp(1),
    _sparse(nullptr),
    _dense(nullptr) {
    allocate(MAP_SIZE);
  }

  explicit FixedSizeSparseMap(const size_t max_size,
                              const Value initial_value) :
    _map_size(0),
    _initial_value(initial_value),
    _data(nullptr),
    _size(0),
    _timestamp(1),
    _sparse(nullptr),
    _dense(nullptr) {
    allocate(max_size);
  }

  FixedSizeSparseMap(const FixedSizeSparseMap&) = delete;
  FixedSizeSparseMap& operator= (const FixedSizeSparseMap& other) = delete;

  FixedSizeSparseMap(FixedSizeSparseMap&& other) :
    _map_size(other._map_size),
    _initial_value(other._initial_value),
    _data(std::move(other._data)),
    _size(other._size),
    _timestamp(other._timestamp),
    _sparse(std::move(other._sparse)),
    _dense(std::move(other._dense)) {
    other._data = nullptr;
    other._sparse = nullptr;
    other._dense = nullptr;
  }

  ~FixedSizeSparseMap() = default;

  size_t capacity() const {
    return _map_size;
  }

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

  void setMaxSize(const size_t max_size) {
    if ( max_size > _map_size ) {
      freeInternalData();
      allocate(max_size);
    }
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

  void freeInternalData() {
    _size = 0;
    _timestamp = 0;
    _data = nullptr;
    _sparse = nullptr;
    _dense = nullptr;
  }

 private:
  inline SparseElement* find(const Key key) const {
    ASSERT(_size < _map_size);
    size_t hash = key & ( _map_size - 1 );
    while ( _sparse[hash].timestamp == _timestamp ) {
      ASSERT(_sparse[hash].element);
      if ( _sparse[hash].element->key == key ) {
        return &_sparse[hash];
      }
      hash = (hash + 1) & ( _map_size - 1 );
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

  void allocate(const size_t size) {
    if ( _data == nullptr ) {
      _map_size = align_to_next_power_of_two(size);
      _data = std::make_unique<uint8_t[]>(
        _map_size * sizeof(MapElement) + _map_size * sizeof(SparseElement));
      _size = 0;
      _timestamp = 1;
      _sparse = reinterpret_cast<SparseElement*>(_data.get());
      _dense = reinterpret_cast<MapElement*>(_data.get() +  + sizeof(SparseElement) * _map_size);
      memset(_data.get(), 0, _map_size * (sizeof(MapElement) + sizeof(SparseElement)));
    }

  }

  size_t align_to_next_power_of_two(const size_t size) const {
    return std::pow(2.0, std::ceil(std::log2(static_cast<double>(size))));
  }

  size_t _map_size;
  const Value _initial_value;
  std::unique_ptr<uint8_t[]> _data;

  size_t _size;
  size_t _timestamp;
  SparseElement* _sparse;
  MapElement* _dense;
};

template <typename Key = Mandatory,
          typename Value = Mandatory>
class DynamicSparseMap {

  struct MapElement {
    Key key;
    Value value;
  };

  struct SparseElement {
    MapElement* element;
    size_t timestamp;
  };

 public:
  explicit DynamicSparseMap() :
    _capacity(0),
    _initial_value(),
    _data(nullptr),
    _size(0),
    _timestamp(1),
    _sparse(nullptr),
    _dense(nullptr) {
  }

  DynamicSparseMap(const DynamicSparseMap&) = delete;
  DynamicSparseMap& operator= (const DynamicSparseMap& other) = delete;

  DynamicSparseMap(DynamicSparseMap&& other) :
    _capacity(other._capacity),
    _initial_value(other._initial_value),
    _data(std::move(other._data)),
    _size(other._size),
    _timestamp(other._timestamp),
    _sparse(std::move(other._sparse)),
    _dense(std::move(other._dense)) {
    other._data = nullptr;
    other._sparse = nullptr;
    other._dense = nullptr;
  }

  ~DynamicSparseMap() = default;

  size_t capacity() const {
    return _capacity;
  }

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

  void initialize(size_t capacity) {
    allocate(align_to_next_power_of_two(capacity));
  }

  bool contains(const Key key) const {
    SparseElement* s = find(key, _sparse, _capacity);
    return containsValidElement(key, s);
  }

  Value& operator[] (const Key key) {
    SparseElement* s = find(key, _sparse, _capacity);
    if ( containsValidElement(key, s) ) {
      ASSERT(s->element);
      return s->element->value;
    } else {
      if ( _size + 1 > _capacity / 5UL ) {
        grow();
        s = find(key, _sparse, _capacity);
      }
      return addElement(key, _initial_value, s, _dense, _size)->value;
    }
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key));
    return find(key, _sparse, _capacity)->element->value;
  }

  const Value* get_if_contained(const Key key) const {
    SparseElement* s = find(key, _sparse, _capacity);
    if ( containsValidElement(key, s) ) {
      return &s->element->value;
    } else {
      return nullptr;
    }
  }

  void clear() {
    _size = 0;
    ++_timestamp;
  }

  void freeInternalData() {
    _size = 0;
    _timestamp = 0;
    _data = nullptr;
    _sparse = nullptr;
    _dense = nullptr;
  }

  size_t size_in_bytes() const {
    return _capacity * (sizeof(SparseElement) + sizeof(MapElement));
  }

 private:
  inline SparseElement* find(const Key key,
                             SparseElement* sparse,
                             const size_t capacity) const {
    size_t hash = key & ( capacity - 1 );
    while ( sparse[hash].timestamp == _timestamp ) {
      ASSERT(sparse[hash].element);
      if ( sparse[hash].element->key == key ) {
        return &sparse[hash];
      }
      hash = (hash + 1) & ( capacity - 1 );
    }
    return &sparse[hash];
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
                                SparseElement* s,
                                MapElement* dense,
                                size_t& size) {
    dense[size] = MapElement { key, value };
    *s = SparseElement { &dense[size++], _timestamp };
    return s->element;
  }

  void allocate(const size_t size) {
    if ( _data == nullptr ) {
      _capacity = align_to_next_power_of_two(size);
      _data = std::make_unique<uint8_t[]>(
        _capacity * sizeof(MapElement) + _capacity * sizeof(SparseElement));
      _size = 0;
      _timestamp = 1;
      _sparse = reinterpret_cast<SparseElement*>(_data.get());
      _dense = reinterpret_cast<MapElement*>(_data.get() + sizeof(SparseElement) * _capacity);
      memset(_data.get(), 0, _capacity * (sizeof(MapElement) + sizeof(SparseElement)));
    }
  }

  void grow() {
    const size_t capacity = 2UL * _capacity;
    std::unique_ptr<uint8_t[]> data = std::make_unique<uint8_t[]>(
      capacity * sizeof(MapElement) + capacity * sizeof(SparseElement));
    SparseElement* sparse = reinterpret_cast<SparseElement*>(data.get());
    MapElement* dense = reinterpret_cast<MapElement*>(data.get() + sizeof(SparseElement) * capacity);
    memset(data.get(), 0, capacity * (sizeof(MapElement) + sizeof(SparseElement)));

    rehash(sparse, dense, capacity);

    _data = std::move(data);
    _sparse = sparse;
    _dense = dense;
    _capacity = capacity;
  }

  void rehash(SparseElement* sparse, MapElement* dense, const size_t capacity) {
    size_t size = 0;
    for ( const MapElement& element : *this ) {
      SparseElement* slot = find(element.key, sparse, capacity);
      addElement(element.key, element.value, slot, dense, size);
    }
    ASSERT(size == _size);
  }

  constexpr size_t align_to_next_power_of_two(const size_t size) const {
    return std::pow(2.0, std::ceil(std::log2(static_cast<double>(size))));
  }

  size_t _capacity;
  const Value _initial_value;
  std::unique_ptr<uint8_t[]> _data;

  size_t _size;
  size_t _timestamp;
  SparseElement* _sparse;
  MapElement* _dense;
};

struct EmptyStruct { };

template<typename Key> using DynamicSparseSet = DynamicSparseMap<Key, EmptyStruct>;

template <typename Key = Mandatory,
          typename Value = Mandatory>
class DynamicFlatMap {

  struct MapElement {
    Key key;
    Value value;
    uint32_t timestamp;
  };

 public:
  explicit DynamicFlatMap() :
    _data(nullptr),
    _size(0),
    _capacity(0),
    _timestamp(1) {
  }

  DynamicFlatMap(const DynamicFlatMap&) = delete;
  DynamicFlatMap& operator= (const DynamicFlatMap& other) = delete;

  DynamicFlatMap(DynamicFlatMap&& other) :
    _data(std::move(other._data)),
    _size(other._size),
    _capacity(other._capacity),
    _timestamp(other._timestamp) {
    other._data = nullptr;
    other._size = 0;
    other._capacity = 0;
  }

  ~DynamicFlatMap() = default;

  size_t capacity() const {
    return _capacity;
  }

  size_t size() const {
    return _size;
  }

  void initialize(size_t capacity) {
    allocate(align_to_next_power_of_two(capacity));
  }

  bool contains(const Key key) const {
    MapElement* s = find(key);
    return containsValidElement(key, s);
  }

  Value& operator[] (const Key key) {
    MapElement* s = find(key);
    if (containsValidElement(key, s)) {
      return s->value;
    } else {
      if (_size + 1 > _capacity / 5UL) {
        grow();
        s = find(key);
      }
      *s = MapElement { key, Value(), _timestamp };
      _size++;
      return s->value;
    }
  }

  const Value& get(const Key key) const {
    ASSERT(contains(key));
    return find(key)->value;
  }

  const Value* get_if_contained(const Key key) const {
    MapElement* s = find(key);
    if (containsValidElement(key, s)) {
      return &s->value;
    } else {
      return nullptr;
    }
  }

  void clear() {
    _size = 0;
    ++_timestamp;
  }

  void freeInternalData() {
    _size = 0;
    _timestamp = 0;
    _data = nullptr;
  }

  size_t size_in_bytes() const {
    return _capacity * sizeof(MapElement);
  }

 private:
  inline MapElement* find(const Key key) const {
    size_t hash = key & (_capacity - 1);
    while (_data[hash].timestamp == _timestamp ) {
      if (_data[hash].key == key) {
        return &_data[hash];
      }
      hash = (hash + 1) & (_capacity - 1);
    }
    return &_data[hash];
  }

  inline bool containsValidElement(const Key key,
                                   const MapElement* s) const {
    unused(key);
    ASSERT(s);
    const bool is_contained = s->timestamp == _timestamp;
    ASSERT(!is_contained || s->key == key);
    return is_contained;
  }

  void allocate(const size_t capacity) {
    _capacity = capacity;
    _data = std::make_unique<MapElement[]>(_capacity);
    memset(_data.get(), 0, _capacity * sizeof(MapElement));
  }

  void grow() {
    const size_t capacity = 2UL * _capacity;
    std::unique_ptr<MapElement[]> old_data = std::move(_data);
    const MapElement* old_data_begin = old_data.get();
    const MapElement* old_data_end = old_data.get() + _capacity;
    allocate(capacity);

    rehash(old_data_begin, old_data_end);
  }

  void rehash(const MapElement* old_data_begin, const MapElement* old_data_end) {
    for (const MapElement* element = old_data_begin; element < old_data_end; ++element) {
      if (element->timestamp == _timestamp) {
        MapElement* slot = find(element->key);
        *slot = MapElement { element->key, element->value, _timestamp };
      }
    }
  }

  constexpr size_t align_to_next_power_of_two(const size_t size) const {
    return std::pow(2.0, std::ceil(std::log2(static_cast<double>(size))));
  }

  std::unique_ptr<MapElement[]> _data;
  size_t _size;
  size_t _capacity;
  uint32_t _timestamp;
};

} // namespace ds
} // namespace mt_kahypar