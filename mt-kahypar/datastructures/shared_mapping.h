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
#pragma once

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <cmath>

#include "tbb/concurrent_vector.h"

#include "kahypar/macros.h"
#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename Identifier = Mandatory,
          typename Key = Mandatory,
          typename Value = Mandatory>
class SharedMapping {

  static constexpr Identifier INVALID_IDENTIFIER = std::numeric_limits<Identifier>::max();

  struct Element {
    Element() :
      _id(INVALID_IDENTIFIER),
      _value() { }

    Element(const Identifier id, const Value value) :
      _id(id),
      _value(value) { }

    bool acquire(const Identifier id, const Value value) {
      Identifier expected = INVALID_IDENTIFIER;
      Identifier desired = id;
      if ( __atomic_compare_exchange(&_id, &expected, &desired,
          false, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED) ) {
        _value = value;
        return true;
      }
      return false;
    }

    void release(const Identifier id) {
      unused(id);
      ASSERT(id == _id);
      __atomic_store_n(&_id, INVALID_IDENTIFIER, __ATOMIC_RELAXED);
    }

    Identifier _id;
    Value _value;
  };

  struct MapElement {
    MapElement() :
      primary_elem(),
      size(0),
      secondary_elems() { }

    Element primary_elem;
    size_t size;
    tbb::concurrent_vector<Element> secondary_elems;
  };

 public:
  explicit SharedMapping(const Key max_size) :
    _elements(max_size)  { }

  void add(const Identifier id, const Key key, const Value value) {
    ASSERT(!contains(id, key));
    ASSERT(static_cast<size_t>(key) < _elements.size());
    MapElement& map_elem = _elements[key];
    if ( unlikely(!map_elem.primary_elem.acquire(id, value)) ) {
      bool success = false;
      for ( Element& elem : map_elem.secondary_elems ) {
        if ( elem.acquire(id, value) ) {
          success = true;
          break;
        }
      }
      if ( !success ) {
        map_elem.secondary_elems.emplace_back(id, value);
      }
      __atomic_fetch_add(&map_elem.size, 1, __ATOMIC_RELAXED);
    }
  }

  void remove(const Identifier id, const Key key) {
    ASSERT(contains(id, key));
    ASSERT(static_cast<size_t>(key) < _elements.size());
    MapElement& map_elem = _elements[key];
    if ( __atomic_load_n(&map_elem.primary_elem._id , __ATOMIC_RELAXED) == id ) {
      map_elem.primary_elem.release(id);
    } else if ( __atomic_load_n(&map_elem.size, __ATOMIC_RELAXED) > 0 ) {
      for ( Element& elem : map_elem.secondary_elems ) {
        if ( __atomic_load_n(&elem._id, __ATOMIC_RELAXED) == id ) {
          elem.release(id);
          __atomic_fetch_sub(&map_elem.size, 1, __ATOMIC_RELAXED);
          break;
        }
      }
    }
  }

  const Value* get_if_contained(const Identifier id, const Key key) const {
    ASSERT(static_cast<size_t>(key) < _elements.size());
    const MapElement& map_elem = _elements[key];
    if ( __atomic_load_n(&map_elem.primary_elem._id , __ATOMIC_RELAXED) == id ) {
      return &map_elem.primary_elem._value;
    } else if ( __atomic_load_n(&map_elem.size, __ATOMIC_RELAXED) > 0 ) {
      for ( const Element& elem : map_elem.secondary_elems ) {
        if ( __atomic_load_n(&elem._id, __ATOMIC_RELAXED) == id ) {
          return &elem._value;
        }
      }
    }
    return nullptr;
  }

  bool contains(const Identifier id, const Key key) const {
    ASSERT(static_cast<size_t>(key) < _elements.size());
    const MapElement& map_elem = _elements[key];
    if ( __atomic_load_n(&map_elem.primary_elem._id , __ATOMIC_RELAXED) == id ) {
      return true;
    } else if ( __atomic_load_n(&map_elem.size, __ATOMIC_RELAXED) > 0 ) {
      for ( const Element& elem : map_elem.secondary_elems ) {
        if ( __atomic_load_n(&elem._id, __ATOMIC_RELAXED) == id ) {
          return true;
        }
      }
    }
    return false;
  }

  void clear() {
    tbb::parallel_for(0UL, _elements.size(), [&](const size_t i) {
      _elements[i].primary_elem._id = INVALID_IDENTIFIER;
      _elements[i].size = 0;
      _elements[i].secondary_elems.clear();
    });
  }

 private:
  vec<MapElement> _elements;
};

} // namespace ds
} // namespace mt_kahypar