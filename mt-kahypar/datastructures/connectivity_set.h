/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>
#include <type_traits>
#include <limits>
#include <cassert>

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
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

  ConnectivitySets(const HyperedgeID numEdges, const PartitionID k) : k(k),
                                                                      numEdges(numEdges),
                                                                      numBlocksPerHyperedge(k / BITS_PER_BLOCK + (k % BITS_PER_BLOCK != 0)),
                                                                      bits(numEdges * numBlocksPerHyperedge)
  {
  }

  void add(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
  }

  void remove(const HyperedgeID he, const PartitionID p) {
    toggle(he, p);
  }

  bool contains(const HyperedgeID he, const PartitionID p) const {
    const PartitionID div = p / BITS_PER_BLOCK;
    const PartitionID rem = p % BITS_PER_BLOCK;
    return bits[he * numBlocksPerHyperedge + div].load(std::memory_order_relaxed) & (UnsafeBlock(1) << rem);
  }

  // not threadsafe
  void clear(const HyperedgeID he) {
    for (size_t i = he * numBlocksPerHyperedge; i < (he + 1) * numBlocksPerHyperedge; ++i) {
      bits[i].store(0, std::memory_order_relaxed);
    }
  }

  void reset() {
    for (size_t i = 0; i < bits.size(); ++i) {
      bits[i].store(0, std::memory_order_relaxed);
    }
  }

  PartitionID connectivity(const HyperedgeID he) const {
    PartitionID conn = 0;
    for (size_t i = he * numBlocksPerHyperedge; i < (he + 1) * numBlocksPerHyperedge; ++i) {
      conn += utils::popcount_64(bits[i].load(std::memory_order_relaxed));
    }
    return conn;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    parent->addChild("Connectivity Bit Vector", sizeof(Block) * bits.size());
  }

private:
  using UnsafeBlock = uint64_t;
  static constexpr int BITS_PER_BLOCK = std::numeric_limits<UnsafeBlock>::digits;
  using Block = parallel::IntegralAtomicWrapper<UnsafeBlock>;
  using BlockIterator = parallel::scalable_vector<Block>::const_iterator;


	PartitionID k;
	HyperedgeID numEdges;
	PartitionID numBlocksPerHyperedge;
	parallel::scalable_vector<Block> bits;

	void toggle(const HyperedgeID he, const PartitionID p) {
	  assert(p < k);
	  assert(he < numEdges);
    const size_t div = p / BITS_PER_BLOCK, rem = p % BITS_PER_BLOCK;
    const size_t idx = he * numBlocksPerHyperedge + div;
    assert(idx < bits.size());
	  bits[idx].fetch_xor(UnsafeBlock(1) << rem, std::memory_order_relaxed);
	}

public:

  class Iterator : public std::iterator<std::forward_iterator_tag, PartitionID, std::ptrdiff_t, const PartitionID*, PartitionID> {
  public:
    Iterator(BlockIterator first, PartitionID part, PartitionID k) : currentPartition(part), k(k), firstBlockIt(first) {
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
    PartitionID k;
    BlockIterator firstBlockIt;

    void findNextBit() {
      ++currentPartition;
      UnsafeBlock b = currentBlock()->load(std::memory_order_release);
      while (b >> (currentPartition % BITS_PER_BLOCK) == 0 && currentPartition < k) {
        currentPartition += (BITS_PER_BLOCK - (currentPartition % BITS_PER_BLOCK));   // skip rest of block
        b = currentBlock()->load(std::memory_order_release);
      }
      if (currentPartition < k) {
        currentPartition += utils::lowest_set_bit_64(b >> (currentPartition % BITS_PER_BLOCK));
      } else {
        currentPartition = k;
      }
    }

    BlockIterator currentBlock() const {
      return firstBlockIt + currentPartition / BITS_PER_BLOCK;
    }

  };

	Iterator hyperedgeBegin(const HyperedgeID he) const {
	  return Iterator(bits.begin() + he * numBlocksPerHyperedge, -1, k);
	}

	Iterator hyperedgeEnd(const HyperedgeID he) const {
	  return Iterator(bits.begin() + he * numBlocksPerHyperedge, k-1, k);
	}


  IteratorRange<Iterator> connectivitySet(const HyperedgeID he) const {
    return IteratorRange<Iterator>(hyperedgeBegin(he), hyperedgeEnd(he));
  }

};



}  // namespace ds
}  // namespace mt_kahypar
