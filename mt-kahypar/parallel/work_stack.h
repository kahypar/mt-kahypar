/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "../definitions.h"

namespace mt_kahypar {


template<typename T>
class WorkStack {
public:

  WorkStack(size_t maxNumElements) :
          size(0),
          elements(maxNumElements, T())
  {

  }

  void push_back(const T& el) {
    const size_t old_size = size.fetch_add(1, std::memory_order_acq_rel);
    assert(old_size < elements.size());
    elements[old_size] = el;
  }

  bool try_pop(T& dest) {
    const size_t old_size = size.fetch_sub(1, std::memory_order_acq_rel);
    if (old_size > 0 && old_size < capacity()) {
      dest = elements[old_size - 1];
      return true;
    }
    return false;
  }

  bool empty() const {
    const size_t s = unsafe_size();
    return s == 0 || s >= capacity();
  }

  size_t unsafe_size() const {
    return size.load(std::memory_order_acq_rel);
  }

  size_t capacity() const {
    return elements.size();
  }

  vec<T>& get_underlying_container() {
    return elements;
  }

  void clear() {
    size.store(0);
  }

  void shrink_to_fit() {
    elements.resize(size);
    elements.shrink_to_fit();
  }

  void shuffle() {
    utils::Randomize::instance().parallelShuffleVector(elements, 0, unsafe_size());
  }

  CAtomic<size_t> size;
  vec<T> elements;
};


}