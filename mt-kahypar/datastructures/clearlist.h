/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <vector>

#include "boost/dynamic_bitset.hpp"

namespace mt_kahypar {
namespace ds {
template <class T>
class ClearListSet {
 public:
  ClearListSet(size_t n) :
    contained_keys(),
    set(n) { }

 public:
  void insert(const T e) {
    set.set(e);
    contained_keys.push_back(e);
  }

  void remove(const T e) {
    set.reset(e);
  }

  bool contains(const T e) const {
    return set[e];
  }

  const std::vector<T>& keys() const { return contained_keys; }
  std::vector<T>& keys() { return contained_keys; }

  template <typename F>
  void forKeys(F f) {
    for (T& e : contained_keys)
      if (contains(e))
        f(e);
  }

  void clear() {
    for (T& e : keys())
      remove(e);
    contained_keys.clear();
  }

 private:
  std::vector<T> contained_keys;
  boost::dynamic_bitset<> set;
};


template <class K, class V>
class ClearListMap {
 public:
  ClearListMap(size_t size, V defaultValue) :
    defaultValue(defaultValue),
    contained_keys(),
    map(size, defaultValue) { }

  void insert(const K k, const V v) {
    contained_keys.push_back(k);
    map[k] = v;
  }

  void add(const K k, const V v) {
    if (!contains(k))
      insert(k, v);
    else
      map[k] += v;
  }

  void remove(const K k) {
    map[k] = defaultValue;
  }

  bool contains(const K k) const {
    return map[k] != defaultValue;
  }

  const V & operator[] (const K k) const { return map[k]; }

  V & operator[] (const K k) { return map[k]; }

  const std::vector<K>& keys() const { return contained_keys; }
  std::vector<K>& keys() { return contained_keys; }

  template <typename F>
  void forKeyValuePairs(F f) {
    for (K& k :keys())
      if (contains(k))
        f(k, map[k]);
  }

  void clear() {
    for (K& k : keys())
      remove(k);
    contained_keys.clear();
    // assert( std::all_of( map.begin(), map.end(), [&](const V& v) { return v == defaultValue; } ) );
  }

  void merge(ClearListMap<K, V>& other) {
    for (K& k : other.keys())
      add(k, other[k]);
    other.clear();
  }

 private:
  V defaultValue;
  std::vector<K> contained_keys;
  std::vector<V> map;
};
}  // namespace ds
}  // namespace mt_kahypar
