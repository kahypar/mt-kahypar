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

#include <vector>
#include <cstdint>
#include <functional>
#include <algorithm>
#include <cassert>

namespace mt_kahypar {
namespace ds {

using PosT = uint32_t;
static constexpr PosT invalid_position = std::numeric_limits<PosT>::max();

template<typename KeyT, typename IdT, typename Comparator = std::less<KeyT>, uint32_t arity = 4>
class Heap {

public:
  static_assert(arity > 1);


  explicit Heap(std::vector<IdT>& external_positions) : comp(), heap(), positions(external_positions) {
    heap.reserve(1024);
  }

  IdT top() const {
    return heap[0].id;
  }

  KeyT topKey() const {
    return heap[0].key;
  }

  void deleteTop() {
    assert(!empty());
    positions[heap[0].id] = invalid_position;
    positions[heap.back().id] = 0;
    heap[0] = heap.back();
    heap.pop_back();
    siftDown(0);
  }

  void insert(const KeyT k, const IdT e) {
    assert(!contains(e));
    const PosT pos = size();
    positions[e] = pos;
    heap.push_back({k, e});
    siftUp(pos);
  }

  void remove(const IdT e) {
    assert(!empty() && contains(e));
    PosT pos = positions[e];
    const KeyT removedKey = heap[pos].key, lastKey = heap.back().key;
    positions[e] = invalid_position;
    positions[heap.back().id] = pos;
    heap[pos] = heap.back();
    heap.pop_back();
    if (comp(removedKey, lastKey)) {
      siftUp(pos);
    } else {
      siftDown(pos);
    }
  }

  // assumes semantics of comp = std::less, i.e. we have a MaxHeap and increaseKey moves the element up the tree. if comp = std::greater, increaseKey will still move the element up the tree
  void increaseKey(const IdT e, const KeyT newKey) {
    assert(contains(e));
    const PosT pos = positions[e];
    assert(comp(heap[pos].key, newKey));
    heap[pos].key = newKey;
    siftUp(pos);
  }

  // assumes semantics of comp = std::less, i.e. we have a MaxHeap and increaseKey moves the element up the tree
  void decreaseKey(const IdT e, const KeyT newKey) {
    assert(contains(e));
    const PosT pos = positions[e];
    assert(comp(newKey, heap[pos].key));
    heap[pos].key = newKey;
    siftDown(pos);
  }

  void adjustKey(const IdT e, const KeyT newKey) {
    assert(contains(e));
    const PosT pos = positions[e];
    if (comp(heap[pos].key, newKey)) {
      increaseKey(e, newKey);
    } else {
      decreaseKey(e, newKey);
    }
  }

  bool contains(const IdT e) const {
    assert(fits(e));
    return positions[e] != invalid_position && heap[positions[e]].id == e;
  }

  PosT size() const {
    return static_cast<PosT>(heap.size());
  }

  bool empty() const {
    return size() == 0;
  }

  KeyT key(const PosT pos) const {
    return heap[pos].key;
  }

  IdT id(const PosT pos) const {
    return heap[pos].id;
  }


protected:

  bool fits(const IdT id) const {
    return id < positions.size();
  }

  PosT parent(const PosT pos) const {
    return (pos - 1) / arity;
  }

  PosT firstChild(const PosT pos) const {
    return pos*arity + 1;
  }

  void siftUp(PosT pos) {
    const KeyT k = heap[pos].key;
    const IdT id = heap[pos].id;

    PosT parent_pos = parent(pos);
    while (pos > 0 && comp(heap[parent_pos].key, k)) {    // eliminate pos > 0 check by a sentinel at position zero?
      positions[ heap[parent_pos].id ] = pos;
      heap[pos] = heap[parent_pos];
      pos = parent_pos;
      parent_pos = parent(pos);
    }
    positions[id] = pos;
    heap[pos].id = id;
    heap[pos].key = key;
  }

  void siftDown(PosT pos) {
    const KeyT k = heap[pos].key;
    const IdT id = heap[pos].id;
    const PosT initial_pos = pos;

    PosT first = firstChild(pos);
    while (first < size() && first != pos) {
      PosT largestChild;

      if constexpr (arity > 2) {
        largestChild = first;
        KeyT largestChildKey = heap[largestChild].key;

        // find child with largest key for MaxHeap / smallest key for MinHeap
        const PosT firstInvalid = std::min(size(), firstChild(pos + 1));
        for (PosT c = first + 1; c < firstInvalid; ++c) {
          if ( comp(largestChildKey, heap[c].key) ) {
            largestChildKey = heap[c].key;
            largestChild = c;
          }
        }

        if (comp(largestChildKey, k) || largestChildKey == k) {
          break;
        }

      } else {
        assert(arity == 2);

        const PosT second = std::min(first + 1, size() - 1);    // TODO this branch is not cool. maybe make the while loop condition secondChild(pos) < size() ?
        const KeyT k1 = heap[first].key, k2 = heap[second].key;
        const bool c2IsLarger = comp(k2, k1);

        if (comp(k, c2IsLarger ? k2 : k1)) {
          break;
        }

        largestChild = c2IsLarger ? second : first;
      }

      positions[ heap[largestChild].id ] = pos;
      heap[pos] = heap[largestChild];
      pos = largestChild;
      first = firstChild(pos);
    }

    if (pos != initial_pos) {
      positions[id] = pos;
      heap[pos].key = key;
      heap[pos].id = id;
    }
  }

  struct HeapElement {
    KeyT key;
    IdT id;
  };

  Comparator comp;                // comp(heap[parent(pos)].key, heap[pos].key) returns true if the element at pos should move upward --> comp = std::less for MaxHeaps
                                  // similarly comp(heap[child(pos)].key, heap[pos].key) returns false if the element at pos should move downward
  std::vector<HeapElement> heap;
  std::vector<PosT>& positions;    // this is ref because BlockQueues need to share handles
};


// used to initialize handles in ExclusiveHandleHeap before handing a ref to Heap
template<typename IdT>
struct HandlesPBase {
  explicit HandlesPBase(size_t n) : handles(n, invalid_position) { }
  std::vector<IdT> handles;
};

template<typename KeyT, typename IdT, typename Comparator, uint32_t arity>
class ExclusiveHandleHeap : protected HandlesPBase<IdT>, public Heap<KeyT, IdT, Comparator, arity> {
public:
  explicit ExclusiveHandleHeap(size_t nHandles) : HandlesPBase<IdT>(nHandles), Heap<KeyT, IdT, Comparator, arity>(this->handles) { }
};


// using BlockQueue = ExclusiveHandleHeap<PartitionID, GainT, std::less<GainT>>;
// using VertexQueue = MaxHeap<HypernodeID, GainT>;   // these need external handles

// we will probably not use these. and instead just do a scan over k entries. should be fast enough
// using DQueue = MaxHeap<PartitionID, HyperedgeWeight>;


}
}
