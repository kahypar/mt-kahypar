/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <atomic>
#include <type_traits>
#include <limits>
#include <cassert>

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/bit_ops.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"

#include "mt-kahypar/macros.h"
#include "hypergraph_common.h"

namespace mt_kahypar {
namespace ds {


/**
 *      The connectivity set of a hyperedge is the set of parts of the partition, that it has pins in.
 *      For each hyperedge we maintain its connectivity set in a packed format (std::vector<uint64_t>)
 *      and implement the necessary bitset functionality ourselves, i.e. add, remove, contains, clear, iteration.
 *      That is because we want atomic updates to support safe parallel modification of the partition.
 *      Adding/removing a part are both implemented as a toggle of the corresponding bit, in case an add and
 *      a remove operation are interweaved. However, this means the user must ensure that no two threads simultaneously try
 *      to add a part. One correct way is to keep an atomic count of pins for each hyperedge and part. Then only the thread
 *      raising the counter from zero to one performs the add, and only the thread decreasing the counter from one to zero
 *      performs the removal.
 */
class ConnectivitySets {
public:

  static constexpr bool debug = false;
  using UnsafeBlock = uint64_t;

  ConnectivitySets() :
    _k(0),
    _num_hyperedges(0),
    _num_blocks_per_hyperedge(0),
    _bits() { }

  ConnectivitySets(const HyperedgeID num_hyperedges,
                   const PartitionID k,
                   const bool assign_parallel = true) :
    _k(k),
    _num_hyperedges(num_hyperedges),
    _num_blocks_per_hyperedge(k / BITS_PER_BLOCK + (k % BITS_PER_BLOCK != 0)),
    _bits() {
      if ( num_hyperedges > 0 ) {
        _bits.resize("Refinement", "connectivity_set",
          static_cast<size_t>(num_hyperedges) * _num_blocks_per_hyperedge, true, assign_parallel);
      }
    }

  void add(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
  }

  void remove(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
  }

  bool contains(const HyperedgeID he, const PartitionID p) const {
    const size_t div = p / BITS_PER_BLOCK;
    const size_t rem = p % BITS_PER_BLOCK;
    const size_t pos = static_cast<size_t>(he) * _num_blocks_per_hyperedge + div;
    return _bits[pos].load(std::memory_order_relaxed) & (UnsafeBlock(1) << rem);
  }

  // not threadsafe
  void clear(const HyperedgeID he) {
    const size_t start = static_cast<size_t>(he) * _num_blocks_per_hyperedge;
    const size_t end = ( static_cast<size_t>(he) + 1 ) * _num_blocks_per_hyperedge;
    for (size_t i = start; i < end; ++i) {
      _bits[i].store(0, std::memory_order_relaxed);
    }
  }

  void reset() {
    for (size_t i = 0; i < _bits.size(); ++i) {
      _bits[i].store(0, std::memory_order_relaxed);
    }
  }

  PartitionID connectivity(const HyperedgeID he) const {
    PartitionID conn = 0;
    const size_t start = static_cast<size_t>(he) * _num_blocks_per_hyperedge;
    const size_t end = ( static_cast<size_t>(he) + 1 ) * _num_blocks_per_hyperedge;
    for (size_t i = start; i < end; ++i) {
      conn += utils::popcount_64(_bits[i].load(std::memory_order_relaxed));
    }
    return conn;
  }

  void freeInternalData() {
    parallel::free(_bits);
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    parent->addChild("Connectivity Bit Vector", sizeof(Block) * _bits.size());
  }

  static size_t num_elements(const HyperedgeID num_hyperedges,
                             const PartitionID k) {
    return static_cast<size_t>(num_hyperedges) * (k / BITS_PER_BLOCK + (k % BITS_PER_BLOCK != 0));
  }

private:
  static constexpr int BITS_PER_BLOCK = std::numeric_limits<UnsafeBlock>::digits;
  using Block = parallel::IntegralAtomicWrapper<UnsafeBlock>;
  using BlockIterator = Array<Block>::const_iterator;


	PartitionID _k;
	HyperedgeID _num_hyperedges;
	PartitionID _num_blocks_per_hyperedge;
	Array<Block> _bits;

	void toggle(const HyperedgeID he, const PartitionID p) {
	  assert(p < _k);
	  assert(he < _num_hyperedges);
    const size_t div = p / BITS_PER_BLOCK, rem = p % BITS_PER_BLOCK;
    const size_t idx = static_cast<size_t>(he) * _num_blocks_per_hyperedge + div;
    assert(idx < _bits.size());
	  _bits[idx].fetch_xor(UnsafeBlock(1) << rem, std::memory_order_relaxed);
	}

public:

  class Iterator : public std::iterator<std::forward_iterator_tag, PartitionID, std::ptrdiff_t, const PartitionID*, PartitionID> {
  public:
    Iterator(BlockIterator first, PartitionID part, PartitionID k) : currentPartition(part), _k(k), firstBlockIt(first) {
      findNextBit();
    }

    PartitionID operator*() const {
      return currentPartition;
    }

    Iterator& operator++() {
      findNextBit();
      return *this;
    }

    Iterator operator++(int ) {
      const Iterator res = *this;
      findNextBit();
      return res;
    }

    bool operator==(const Iterator& o) const {
      return currentPartition == o.currentPartition && firstBlockIt == o.firstBlockIt;
    }

    bool operator!=(const Iterator& o) const {
      return !operator==(o);
    }

  private:
    PartitionID currentPartition;
    PartitionID _k;
    BlockIterator firstBlockIt;

    void findNextBit() {
      ++currentPartition;
      UnsafeBlock b = currentBlock()->load(std::memory_order_release);
      while (b >> (currentPartition % BITS_PER_BLOCK) == 0 && currentPartition < _k) {
        currentPartition += (BITS_PER_BLOCK - (currentPartition % BITS_PER_BLOCK));   // skip rest of block
        b = currentBlock()->load(std::memory_order_release);
      }
      if (currentPartition < _k) {
        currentPartition += utils::lowest_set_bit_64(b >> (currentPartition % BITS_PER_BLOCK));
      } else {
        currentPartition = _k;
      }
    }

    BlockIterator currentBlock() const {
      return firstBlockIt + currentPartition / BITS_PER_BLOCK;
    }

  };

	Iterator hyperedgeBegin(const HyperedgeID he) const {
	  return Iterator(_bits.cbegin() + static_cast<size_t>(he) * _num_blocks_per_hyperedge, -1, _k);
	}

	Iterator hyperedgeEnd(const HyperedgeID he) const {
	  return Iterator(_bits.cbegin() + static_cast<size_t>(he) * _num_blocks_per_hyperedge, _k-1, _k);
	}


  IteratorRange<Iterator> connectivitySet(const HyperedgeID he) const {
    return IteratorRange<Iterator>(hyperedgeBegin(he), hyperedgeEnd(he));
  }

};



}  // namespace ds
}  // namespace mt_kahypar
